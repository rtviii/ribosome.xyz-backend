import argparse
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
from dotenv import load_dotenv
from .bsite_mixed import BindingSite, open_structure
from .transpose_ligand import  SeqMatch 


""" The goal is to have a module that, 
given two files (assume both same typep, both exist for now)
and a residue range:
1. retrieves them
2. seq-aligns them
3. maps the range on the alignment
4. backtracks the new ranges to original chains
5. clips the backtracked ranges out of the chains in pymol
6. superimpose the individual snippets
"""



def ranged_super(
	polymer_class: str,
	src_struct   : str,
	tgt_struct   : str,
	rng        : Tuple[int,int],

)->Tuple[str, Tuple[int,int], str, Tuple[int,int]]:

	"""Return a bundle of path + mapped range for a source and a target structure
	for a given polymer class. Feed this into pymol 
	to chop up on the ranges and superimpose resultant snippets.""" 

	# assert(rng[1]-rng[0]>20)

	rstart    ,rend = rng



	json_src = open_structure(src_struct.upper(),'json')
	json_tgt = open_structure(tgt_struct.upper(),'json')


	tgt_chainid, src_chainid = [ None, None ]
	tgt_seq    , src_seq     = [ None, None ]


	for chain in [ *json_src['proteins'], *json_src['rnas'] ]:
		if polymer_class in chain['nomenclature']:
			src_chainid = chain['auth_asym_id']
			src_seq     = chain['entity_poly_seq_one_letter_code']

	for chain in [ *json_tgt['proteins'], *json_tgt['rnas'] ]:
		if polymer_class in chain['nomenclature']:
			tgt_chainid = chain['auth_asym_id']
			tgt_seq     = chain['entity_poly_seq_one_letter_code']

	if None in [tgt_chainid, src_chainid, tgt_seq, src_seq]: print(f"Couldn't find a polymer of matching class {polymer_class}"); exit(1)


	ixs= [*range(rstart,rend)]
	sm = SeqMatch(src_seq,tgt_seq, ixs)


	target_range = ( sm.tgt_ids[0],sm.tgt_ids[-1] )
	source_range = ( sm.src_ids[0],sm.src_ids[-1] )

	src_chain_path = os.path.join(os.environ.get( "STATIC_ROOT" ), src_struct.upper(), "CHAINS", "{}_STRAND_{}.cif".format(src_struct.upper(),src_chainid))
	tgt_chain_path = os.path.join(os.environ.get( "STATIC_ROOT" ), tgt_struct.upper(), "CHAINS", "{}_STRAND_{}.cif".format(tgt_struct.upper(),tgt_chainid))

	return (src_chain_path, source_range, tgt_chain_path, target_range )




if __name__ =="__main__":
	load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
	STATIC_ROOT = os.environ.get('STATIC_ROOT')
	[print(x) for x in ranged_super('uL24', '3j7z','5afi', ( 20,50 ))]

# 	prs = argparse.ArgumentParser()


# 	prs.add_argument('-s', '--source_struct' , type=str, required=True)
# 	prs.add_argument('-t', '--target_struct' , type=str, required=True)
# 	prs.add_argument('-c', '--poly_class'    , type=str, required=True)
# 	prs.add_argument('-r', '--residue_range' , type=str, required=True)

# 	args = prs.parse_args()

# 	src_struct      =            args.source_struct.upper()
# 	tgt_struct      =            args.target_struct.upper()
# 	poly_class      =            args.poly_class
# 	rstart    ,rend = [* map(int,args.residue_range        .split("-")) ]

# 	cif_src  = open_structure(src_struct.upper(),'cif')
# 	cif_tgt  = open_structure(tgt_struct.upper(),'cif')

# 	json_src = open_structure(src_struct.upper(),'json')
# 	json_tgt = open_structure(tgt_struct.upper(),'json')


# 	tgt_chainid, src_chainid = [ None, None ]
# 	tgt_seq    , src_seq     = [ None, None ]


# 	for chain in [ *json_src['proteins'], *json_src['rnas'] ]:
# 		if poly_class in chain['nomenclature']:
# 			src_chainid = chain['auth_asym_id']
# 			src_seq     = chain['entity_poly_seq_one_letter_code']

# 	for chain in [ *json_tgt['proteins'], *json_tgt['rnas'] ]:
# 		if poly_class in chain['nomenclature']:
# 			tgt_chainid = chain['auth_asym_id']
# 			tgt_seq     = chain['entity_poly_seq_one_letter_code']

# 	if None in [tgt_chainid, src_chainid, tgt_seq, src_seq]: print(f"Couldn't find a polymer of matching class {poly_class}"); exit(1)


# 	ixs= [*range(rstart,rend)]
# 	sm = SeqMatch(src_seq,tgt_seq, ixs)


# 	# print(SeqMatch.hl_ixs(sm.src_aln, ixs))
# 	# print(SeqMatch.hl_ixs(sm.tgt_aln, ixs))
# 	# print("\n <----> \n")
# 	# print(SeqMatch.hl_ixs(sm.src, sm.src_ids))
# 	# print(SeqMatch.hl_ixs(sm.tgt, sm.tgt_ids))
# 	# print("\n")

# 	target_range = [ sm.tgt_ids[0],sm.tgt_ids[-1] ]
# 	source_range = [ sm.src_ids[0],sm.src_ids[-1] ]

# 	src_chain_path = os.path.join(STATIC_ROOT, src_struct.upper(), "CHAINS", "{}_STRAND_{}.cif".format(src_struct.upper(),src_chainid))
# 	tgt_chain_path = os.path.join(STATIC_ROOT, tgt_struct.upper(), "CHAINS", "{}_STRAND_{}.cif".format(tgt_struct.upper(),tgt_chainid))


# 	print(src_chain_path)
# 	print(source_range)
# 	print(tgt_chain_path)
# 	print(target_range)

# 	# tgt_struct  = open_structure(tgt_struct,'cif')