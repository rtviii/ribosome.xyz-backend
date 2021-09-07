from ctypes import alignment
import dataclasses
import operator
import json
from pprint import pprint
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
import numpy as np
import ciftools.binding_site as bsite
flatten = itertools.chain.from_iterable
n1      = np.array


class SeqMatch:

	def __init__(self,sourceseq:str,targetseq:str, source_residues:List[int]) -> None:
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
			print(len(self.tgt_aln))
			print(self.tgt_aln)
			print(self.aligned_ids)
			print("looking for res", aln_resid)
			if self.tgt_aln[aln_resid] == '-':
				print("Omitting aligned residue ", aln_resid)
				continue
			self.tgt_ids.append(self.backwards_match(self.tgt_aln,aln_resid))

	def backwards_match(self, alntgt:str, resid:int):
		"""Returns the target-sequence  index of a residue in the (aligned) target sequence"""
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
#* match matchable chains from target to prediction

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



def transpose_ligand(source_struct, target_struct,chemid):

	source_struct = source_struct.upper()
	target_struct = target_struct.upper()
	chemid        = chemid.upper()



	with open(os.path.join(os.getenv('STATIC_ROOT'), source_struct, f'LIGAND_{chemid}.json'), 'rb') as infile:
		data = json.load(infile)
	bs = bsite.BindingSite(data)


	origin_chains = {
	}

	target_chains = {
	}

	#* For every chain in a ligand file, if it has nomenclature, append its residues, strand and sequence
	for chain in bs.data:
		if len(bs.data[chain]['nomenclature'] ) <1:
			continue
		else:
			resids :List[int] = [
				resid for  resid in [*map(lambda x : x['residue_id'], bs.data[chain]['residues'])]
			]
			origin_chains[bs.data[chain]['nomenclature'][0]] = {
				'strand': chain,
				'seq'   : bs.data[chain]['sequence'],
				'ids'   : resids
			}

	for nom in origin_chains:
		name_matches = []
		cypher       = f"""match (n:RibosomeStructure {{rcsb_id:"{target_struct}"}})-[]-(c)-[]-(r {{class_id:"{nom}"}}) return c.entity_poly_seq_one_letter_code, c.entity_poly_strand_id, c.asym_ids"""
		response     = bsite._neoget(cypher)
		if len( response )  < 1:
			print(f"No chain-class matches for {nom} in {target_struct} in the database.")
			continue
		else:
			match = response[0]

		seq                 = match[0]
		strand              = match[1]
		asymid              = match[2][0]

		target_chains[nom] ={
			'seq'   : seq,
			'strand': strand,
			'asymid': asymid,
		}
	#! """Only the chains with nomenclature matches in source and origin make their way into the prediction file """



	prediction ={}


	for name in origin_chains:
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
				"src"    : src,
				"src_ids": src_ids,
				"strand" : origin_chains[name]['strand']
			},
			"target":{
				"tgt"    : tgt,
				"tgt_ids": tgt_ids,
				'strand' : target_chains[name]['strand']
			},
			"alignment" :{
				"aln_ids": aln_ids,
				"src_aln": src_aln,
				"tgt_aln": tgt_aln,
			},
		}

		print(f"Chain {name}")
		print("Source length:", len(sq.src))
		print("Target length:", len(sq.tgt))
		# print("Aligned ids" , src_ids)
		# print("To ------->" , sq        .aligned_ids  )
		# print("To targets :", sq        .tgt_ids      )


		print("ORG   : \t",SeqMatch.hl_ixs(sq.src    , sq.src_ids    ),"\n")
		# print("ORG AL: \t",SeqMatch.hl_ixs(sq.src_aln, sq.aligned_ids),"\n")
		# print("TGT AL: \t",SeqMatch.hl_ixs(sq.tgt_aln, sq.aligned_ids),"\n")
		print("TGT   : \t",SeqMatch.hl_ixs(sq.tgt    , sq.tgt_ids    ),"\n")

	fname = f'PREDICTION_{chemid}_{source_struct}_{target_struct}.json'

	with open(f'/home/rxz/dev/ribetl/static/{target_struct}/{fname}', 'w') as outfile:
		json.dump(prediction,outfile)
	print("Sucessfully saved prediction {}".format(fname))