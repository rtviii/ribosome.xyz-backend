from ast import expr_context
import builtins
from cmath import log
import pdb
from pydoc import resolve
import struct
from neo4j_connector.views import get_struct
from ribetl.ciftools.neoget import _neoget
import dataclasses
import json
from pprint import pprint
from neo4j import GraphDatabase, Result
from dotenv import load_dotenv
import os
from typing import Dict, List, Tuple, TypedDict, Union, Callable
import operator
import sys
import pandas as pd
from Bio.PDB.MMCIFParser import FastMMCIFParser, MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
import argparse
import itertools
from dataclasses import dataclass,field
from asyncio import run
import itertools
import numpy as np
flatten = itertools.chain.from_iterable
n1      = np.array


load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get('STATIC_ROOT')

#? By electrostatic charge
AMINO_ACIDS = {
"ALA":0,
'ARG':1,
'ASN':0,
'ASP':-1,
'CYS':0,
'GLN':0,
'GLU':-1,
'GLY':0,
'HIS':0,
'ILE':0,
'LEU':0,
'LYS':1,
'MET':0,
'PHE':0,
'PRO':0,
'SER':0,
'THR':0,
'TRP':0,
'TYR':0,
'VAL':0,
'SEC':0,
'PYL':0}

NUCLEOTIDES = ['A','T','C','G','U']

def struct_path(pdbid:str, pftype: str ):

    if pftype == 'cif':
        return os.path.join(STATIC_ROOT, pdbid.upper(),f"{pdbid.upper()}.cif")

    if pftype == 'json':
        return os.path.join(STATIC_ROOT, pdbid.upper(),f"{pdbid.upper()}.json")

def open_structure(pdbid:str, cifpath: str = None)->Structure:
    pdbid = pdbid.upper()
    if cifpath == None: cifpath = os.path.join(STATIC_ROOT, pdbid, f'{pdbid}.cif')
    return FastMMCIFParser(QUIET=True).get_structure(pdbid,cifpath)

def open_ligand(pdbid:str, ligid:str, ligpath: str = None):
    pdbid = pdbid.upper()
    ligid = ligid.upper()

    if ligpath == None: ligpath = os.path.join(STATIC_ROOT, pdbid, f'LIGAND_{ligid}.json')
    with open(ligpath, 'rb') as infile: 
            data = json.load(infile)
    return  data

@dataclass(unsafe_hash=True, order=True)
class ResidueLite     : 

      residue_name     : str               = field(              hash=True ,compare=False)
      residue_id       : int               = field(              hash=True ,compare=True)
      parent_strand_id : str               = field(              hash=True ,compare=False)
      
      @staticmethod
      def res2reslite(r:Residue):
          biopy_id_tuple = r.get_full_id()
          parent_chain   = biopy_id_tuple[2]
          resname        = r.resname
          resid          = r.id[1]
          return ResidueLite(resname,resid,parent_chain)

@dataclass
class BindingSiteChain: 
      sequence        : str
      nomenclature    : List[str]
      asym_ids        : List[str]
      residues        : List[ResidueLite]

class BindingSite:
    def __init__(self,data:Dict[str,BindingSiteChain]) -> None:
        self.data : Dict[str,BindingSiteChain] = data
    def to_json(self,pathtofile:str)->None:
        with open(pathtofile, 'w') as outf: 
            serialized = {}
            for x in self.data.items():
                serialized.update({x[0]: dataclasses.asdict(x[1])})  
            json.dump(serialized,outf)
            print(f"\033[91mSaved  {pathtofile} successfuly.\033[0m")

    def to_csv(self,pathtofile:str)->None:

        k = [
        "chainname",
        "nomenclature",
        "residue_id",
        "residue_name"
        ]

        serialized = dict.fromkeys(k,[])
        pprint(serialized)

def getLigandResIds(ligchemid:str, struct: Structure)->List[Residue]:
    """Returns a list of dictionaries specifying each _ligand_ of type @ligchemid as a biopython-residue inside a given @struct."""
    """*ligchemids are of type https://www.rcsb.org/ligand/IDS"""

    ligandResidues: List[Residue] = list(filter(lambda x: x.get_resname() == ligchemid, list( struct.get_residues() )))

    return ligandResidues

async def matchStrandToClass(pdbid:str, strand_id:str)->List[str]:
    """Request Ban nomenclature classes from the db given a protein's entity_poly_strand_id."""
    with open(struct_path(pdbid, 'json'),'r') as infile:
        profile=  json.load(infile)
    for p in profile['proteins']:
        if p['entity_poly_strand_id'] == strand_id:
            return p['nomenclature']
    for r in profile['rnas']:
        if r['entity_poly_strand_id'] == strand_id:
            return r['nomenclature']
    return []

def get_liglike_ids(pdbid:str)->List[str]:
    pdbid       = pdbid.upper()
    with open(struct_path(pdbid, 'json'), 'r') as file:
        profile = json.load(file)
    liglike_strand_ids = []
    for r in profile['rnas']:
        if r['ligand_like'] == True:
            liglike_strand_ids.append(r['entity_poly_strand_id'])
    for p in profile['proteins']:
        if p['ligand_like'] == True:
            liglike_strand_ids.append(p['entity_poly_strand_id'])
    return liglike_strand_ids

def get_liglike_nbrs(
      polymer: dict,
      ligand_residues: List[Residue],
      struct         : Structure,
    )-> BindingSite  : 

    """KDTree search the neighbors of a given list of residues(which constitue a ligand) 
    and return unique having tagged them with a ban identifier proteins within 5 angstrom of these residues. """

    
    pdbid = struct.get_id().upper()

    with open(struct_path(pdbid,'json'),'r') as strfile:
        profile       = json.load(strfile)
        poly_entities = [*profile['proteins'], *profile['rnas']]

    print(f"Parsing ligand-like polymer.")

    ns           = NeighborSearch(list(struct.get_atoms()))
    nbr_residues = []

    # print(f"Ligand has len(ligand_residues) residue(s).")
    #? Searching Phase
    # a ligand consists of residues
    for lig_res in ligand_residues:
        for atom in lig_res.child_list:
                nbr_residues.extend(ns.search(atom.get_coord(), 10,level='R'))

    #? Filtering phase
    #Convert residues to the dataclass, filter non-unique
    nbr_residues = list(set([* map(ResidueLite.res2reslite, nbr_residues) ]))

    #Filter the ligand itself, water and other special residues 
    nbr_residues = list(filter(lambda resl:resl.residue_name  in [*AMINO_ACIDS.keys(),  *NUCLEOTIDES], nbr_residues))
    nbr_dict     = {}
    chain_names  = list(set(map(lambda _:  _.parent_strand_id, nbr_residues)))

    for c in chain_names:

        for poly_entity in poly_entities:
            if ',' in poly_entity['entity_poly_strand_id']:
                if c in poly_entity['entity_poly_strand_id'].split(','):
                    nomenclature = poly_entity['nomenclature'                   ]
                    strands      = poly_entity['entity_poly_strand_id'          ]
                    asymid       = poly_entity['asym_ids'                       ]
                    seq          = poly_entity['entity_poly_seq_one_letter_code']
            elif c == poly_entity['entity_poly_strand_id']:
                    nomenclature = poly_entity['nomenclature'                   ]
                    strands      = poly_entity['entity_poly_strand_id'          ]
                    asymid       = poly_entity['asym_ids'                       ]
                    seq          = poly_entity['entity_poly_seq_one_letter_code']

        # #!RESP may arrive as a single 3-tuple of strand, asymid and sequence or might be an array of 3-tuples of chains whose strand ids contained the queried name. Then filter.
        # cypher=f"""match (n:RibosomeStructure{{rcsb_id:"{struct.get_id().upper()}"}})-[]-(x) where x.entity_poly_strand_id contains "{c}"
        #     return x.entity_poly_strand_id, x.asym_ids,x.entity_poly_seq_one_letter_code"""
        # resp = _neoget(cypher)
        # nomenclature  = list(flatten(run(matchStrandToClass(pdbid,c))))

        nbr_dict[c]= BindingSiteChain(
            seq,
            nomenclature,
            asymid ,
            sorted([residue for residue in nbr_residues if residue.parent_strand_id == c], key=operator.attrgetter('residue_id')))

    vpprint(nbr_dict)
    return BindingSite(nbr_dict)

def dropions(s:str): return False if "ion" in s[1].lower() else  True

def parse_and_save_ligand(ligid:str, rcsbid:str):

    print(rcsbid.upper())
    print(ligid)
    STATIC_ROOT = os.environ.get('STATIC_ROOT')
    print(STATIC_ROOT)

    outfile_json = os.path.join(STATIC_ROOT,rcsbid.upper(), f'LIGAND_{ligid}.json')


    print("\033[92m GOT PATH {} \033[0m".format(outfile_json))

    struct                   = openStructutre  (rcsbid                   )
    residues : List[Residue] = getLigandResIds (ligid  , struct          )
    bs:BindingSite           = get_ligand_nbrs (ligid ,residues , struct )

    bs.to_json(outfile_json)


if __name__ =="__main__" :

    parser      = argparse.ArgumentParser(                                             )
    parser .add_argument ('-l','--ligand'    , type  = str                )
    parser .add_argument ('-s','--structure' , type  = str ,required =True)
    parser .add_argument ('-V','--verbose'   , action='store_true'        )
    parser .add_argument ('-f','--force'     , action='store_true'        )
    args = parser.parse_args()

    VERBOSE     = args.verbose
    PDBID       = args.structure.upper()
    LIGID       = args.ligand

    print("\t\t\t\033[92m * \033[0m")

    if LIGID == None:
        struct_ligands = n1([ *filter(dropions, get_lig_ids_struct(PDBID) ) ],dtype=object)

        print(f"\tParsing ligands for structure {PDBID}: ")
        print("Found:",end="")
        pprint(struct_ligands)
        if len( struct_ligands ) < 1: print("None found. Exited."); exit(1)
        for l in struct_ligands:
            parse_and_save_ligand(l[0].upper(), PDBID)
    else:
        parse_and_save_ligand(LIGID.upper(),PDBID)

    print("\t\t\t\033[92m - \033[0m")

