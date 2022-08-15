import argparse
from ctypes import alignment
import dataclasses
import operator
import json
import re
import struct
import os
from typing import Dict, List, Tuple, TypedDict, Union, Callable
import operator
import sys
from Bio.PDB.MMCIFParser import FastMMCIFParser, MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.Align import MultipleSeqAlignment
from Bio import pairwise2
import itertools
from asyncio import run
import itertools
from dotenv import load_dotenv
from ribetl.ciftools.bsite_mixed import BindingSite, open_structure
import numpy as np
flatten = itertools.chain.from_iterable
n1      = np.array



class SeqMatch:

	def __init__(self,
		sourceseq:str,
		targetseq:str, 
		source_residues:List[int]) -> None:
		"""A container for origin and target sequences when matching the resiudes of a ligand binding site
		to another protein's sequence through BioSeq's Align
		 """

		#* Computed indices of the ligand-facing in the source sequence.
		self.src     :str      = sourceseq
		self.src_ids:List[int] = source_residues

		#* Indices of predicted residues in target sequence. To be filled.
		self.tgt     :str      = targetseq
		self.tgt_ids:List[int] = []
		
		_            = pairwise2.align.globalxx(self.src,self.tgt, one_alignment_only=True)
		self.src_aln = _[0].seqA
		self.tgt_aln = _[0].seqB

		self.aligned_ids = []

		for src_resid in self.src_ids:
			self.aligned_ids.append(self.forwards_match(self.src_aln,src_resid))

		self.aligned_ids = list(filter(lambda x: x != None, self.aligned_ids ))

		for aln_resid in self.aligned_ids:
			if self.tgt_aln[aln_resid] == '-':
				continue
			self.tgt_ids.append(self.backwards_match(self.tgt_aln,aln_resid))

	def backwards_match(self, alntgt:str, resid:int):
		"""Returns the target-sequence index of a residue in the (aligned) target sequence"""
		if resid > len(alntgt):
			exit(IndexError(f"Passed residue with invalid index ({resid}) to back-match to target.Seqlen:{len(alntgt)}"))
		counter_proper = 0
		for i,char in enumerate(alntgt):
			if i == resid:
				return counter_proper
			if char =='-':
				continue
			else: 
				counter_proper  +=1

	def forwards_match(self,alnsrc:str, resid:int):
		"""Returns the index of a source-sequence residue in the aligned source sequence."""
		count_proper = 0
		for alignment_indx,char in enumerate( alnsrc ):
			if count_proper == resid:
				return alignment_indx
			if char =='-':
				continue
			else: 
				count_proper  +=1

	@staticmethod
	def hl_subseq(sequence:str, subsequence:str, index:int=None):
		"""Highlight subsequence"""
		CRED = '\033[91m'
		CEND = '\033[0m'
		_ = [ ]
		if index != None:
			return sequence[:index-1] + CRED + sequence[index] + CEND +sequence[index+1:]
		for item in re.split(re.compile(f'({subsequence})'),sequence):
			if item == subsequence:
				_.append(CRED + item + CEND)
			else:
				_.append(item)
		return ''.join(_)

	@staticmethod
	def hl_ixs(sequence:str,  ixs:List[int]):
		"""Highlight indices"""
		CRED = '\033[91m'
		CEND = '\033[0m'
		_ = ''
		for i,v in enumerate(sequence):
			if i in ixs: _ += CRED + v +CEND
			else: 	 	 _ += v
		return _


#! Include sequence into the ligand profiles
#! exclude the ligand itself from the transposition
#! match matchable chains from target to prediction

#? Protocol: 
#? for every chain in origin, grab same nomenclature in tgt
#? track ids back, apply forth
#? apply ids 

#3j7z  ERY ---> to 5hl7, 3j9z, 6otl
#6AZ1  PAR ---> 5t2a l.donovani

# * Ecoli structs :  3j7z, 7k00, 6q97, 5j30
# ! yeast : 6z6n, 5mrc, 3j6b,6xir, 4u4n

#? PAR:
# 6az1 , 5tcu, 5iqr, 5el7,4wsd,4l71

#? KIR:
# 5afi, 4v8q, 4v5s,4v5g, 

# source_struct = str(sys.argv[1] ).upper()
# target_struct = str(sys.argv[2] ).upper()
# ligand        = str(sys.argv[3]).upper()

def open_bsite(
	  poly_lig_flag    : str,
	  rcsb_id          : str,
	  auth_asym_id     : Union[str, bool]=False,
	  ligand_chemicalId: Union[str, bool]=False
	)->BindingSite or FileNotFoundError     : 

	if poly_lig_flag == 'polymer':
		with open(os.path.join(os.getenv('STATIC_ROOT' ),rcsb_id,f'POLYMER_{auth_asym_id}.json'), 'rb') as infile:
			data = json.load(infile)
		return BindingSite(data)

	else:
		assert(poly_lig_flag=='ligand')
		with open(os.path.join(os.getenv('STATIC_ROOT' ),rcsb_id,f'LIGAND_{ligand_chemicalId}.json'), 'rb') as infile:
			data = json.load(infile)
		return BindingSite(data)
	
def init_transpose_ligand(
	source_struct: str,
	target_struct: str,
	target_profile:dict,
	binding_site : BindingSite
	)->dict            : 

	origin_chains = {
	}
	target_chains = {
	}

	#? For every chain in a ligand file, if it has nomenclature, append its residues, strand and sequence
	for chain in binding_site.data:

		if len(binding_site[chain]['nomenclature'] ) <1:
			continue

		else:

			resids :List[int] = [
				resid for resid in [*map(lambda x : x['residue_id'], binding_site[chain]['residues'])]
			]

			origin_chains[binding_site.data[chain]['nomenclature'][0]] = {
				# 'strand'      : chain,
				'seq'         : binding_site.data[chain]['sequence'],
				'auth_asym_id': binding_site.data[chain]['auth_asym_id'],
				'ids'         : resids
			}

	target_polymers = [*target_profile['rnas'],*target_profile['proteins']]

	for nom in origin_chains:

		# goal is to look up the chain in the target struct __by nomenclature__
		matches =  [*filter(lambda tgt_poly: nom in tgt_poly['nomenclature'], target_polymers)]
		if len( matches )  < 1:
			continue

		seq          = matches[0]['entity_poly_seq_one_letter_code'] 
		auth_asym_id = matches[0]['auth_asym_id']

		target_chains[nom] ={
			'seq'         : seq,
			'auth_asym_id': auth_asym_id
		}


	prediction ={}

	for name in origin_chains:

		# If no chain with a given nomenclature found in target --> skip it.
		if name not in target_chains:
			continue

		src_ids = origin_chains[name]['ids']

		src     = origin_chains[name]['seq']
		tgt     = target_chains[name]['seq']

		sq      = SeqMatch(src,tgt,src_ids)

		src_aln = sq.src_aln
		tgt_aln = sq.tgt_aln

		aln_ids = sq.aligned_ids
		tgt_ids = sq.tgt_ids


		prediction[name] = {
			"source":{
				"src"         : src,
				"src_ids"     : src_ids,
				"auth_asym_id": origin_chains[name]['auth_asym_id']
			},
			"target":{
				"tgt"         : tgt,
				"tgt_ids"     : tgt_ids,
				'auth_asym_id': target_chains[name]['auth_asym_id']
			},
			"alignment" :{
				"aln_ids": aln_ids,
				"src_aln": src_aln,
				"tgt_aln": tgt_aln,
			},
		}

		# print(f"Chain {name}")
		# print("Source length:", len(sq.src))
		# print("Target length:", len(sq.tgt))
		# print("Aligned ids" , src_ids)
		# print("To ------->" , sq        .aligned_ids  )
		# print("To targets :", sq        .tgt_ids      )


		# print("ORG   : \t",SeqMatch.hl_ixs(sq.src    , sq.src_ids    ),"\n")
		print("ORG AL: \t",SeqMatch.hl_ixs(sq.src_aln, sq.aligned_ids),"\n")
		print("TGT AL: \t",SeqMatch.hl_ixs(sq.tgt_aln, sq.aligned_ids),"\n")
		# print("TGT   : \t",SeqMatch.hl_ixs(sq.tgt    , sq.tgt_ids    ),"\n")

	return prediction


if __name__ =="__main__":


	load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
	STATIC_ROOT = os.environ.get('STATIC_ROOT')
	def root_self(rootname: str = ''):
		"""Returns the rootpath for the project if it's unique in the current folder tree."""
		root = os.path.abspath(__file__)[:os.path.abspath(
		__file__).find(rootname)+len(rootname)]
		sys.path.append(root)

	root_self('ribetl')

	prs = argparse.ArgumentParser()

	prs.add_argument('-src', '--source_structure', type=str, required=True)
	prs.add_argument('-tgt', '--target_structure', type=str, required=True)
	prs.add_argument('-poly', '--polymer', type=str)
	prs.add_argument('-lig' , '--ligand' , type=str)

	args = prs.parse_args()

	SRC_STRUCT = args.source_structure.upper()
	TGT_STRUCT = args.target_structure.upper()

	LIG   = args.ligand if args.ligand is not None else False
	POLY  = args.polymer if args.polymer is not None else False
	if LIG:
		bsite_path = os.path.join(STATIC_ROOT, SRC_STRUCT, f"LIGAND_{LIG}.json")
	elif POLY:
		bsite_path = os.path.join(STATIC_ROOT, SRC_STRUCT, f"POLYMER_{POLY}.json")
	else:
		raise argparse.ArgumentTypeError("Provide either a ligand chem. id or a ligand-like polymer's auth_asym_id.")

	_target_handle_json = open_structure(TGT_STRUCT,'json')
	bsite      = open_bsite(
		"ligand" if LIG else "polymer",
		SRC_STRUCT,
		auth_asym_id      = POLY,
		ligand_chemicalId = LIG)

	prediction = init_transpose_ligand(SRC_STRUCT,TGT_STRUCT,_target_handle_json,bsite)
	fname = f'PREDICTION_{LIG if LIG else POLY}_{SRC_STRUCT}_{TGT_STRUCT}.json'
	with open(os.path.join(os.getenv('STATIC_ROOT' ),TGT_STRUCT,fname), 'w') as outfile:
		json.dump(prediction,outfile)
		print("\033[093mSucessfully saved prediction {}\033[0m".format(fname))