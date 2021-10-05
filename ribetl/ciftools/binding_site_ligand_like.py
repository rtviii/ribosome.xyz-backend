from ast import expr_context
import builtins
from cmath import log
from pathlib import Path
import pdb
from pydoc import resolve
import struct
# from ribetl.ciftools.neoget import _neoget
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

def get_poly_nbrs(
        asym_id      : str,
        residues     : List[Residue],
        struct       : Structure,
      )-> BindingSite: 

    """KDTree search the neighbors of a given list of residues(which constitue a ligand) 
    and return unique having tagged them with a ban identifier proteins within 5 angstrom of these residues. """

    pdbid        = struct.get_id().upper()

    ns           = NeighborSearch(list(struct.get_atoms()))
    nbr_residues = []

    # print(f"Ligand has len(ligand_residues) residue(s).")
    #? Searching Phase
    # a ligand consists of residues
    for poly_res in residues:
        for atom in poly_res.child_list:
                nbr_residues.extend(ns.search(atom.get_coord(), 10,level='R'))

    #? Filtering phase
    #Convert residues to the dataclass, filter non-unique
    nbr_residues = list(set([* map(ResidueLite.res2reslite, nbr_residues) ]))

    #Filter the ligand itself, water and other special residues 
    #i.e. "check that each residue is either an amino-acid or a nucleotide. Everything else is an edge case.
    # nbr_residues = list(filter(lambda res: res.residue_name  in [*AMINO_ACIDS.keys(), *NUCLEOTIDES], nbr_residues))
    nbr_residues = list(filter(lambda res: res.parent_strand_id != asym_id and res.residue_name  in [*AMINO_ACIDS.keys(), *NUCLEOTIDES], nbr_residues))

    nbr_dict     = {}
    chain_names  = list(set(map(lambda _:  _.parent_strand_id, nbr_residues)))

    with open(struct_path(pdbid,'json'),'rb') as strfile:
        profile       = json.load(strfile)
        poly_entities = [*profile['proteins'], *profile['rnas']]

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

        nbr_dict[c]= BindingSiteChain(
            seq,
            nomenclature,
            asymid ,
            sorted([residue for residue in nbr_residues if residue.parent_strand_id == c], key=operator.attrgetter('residue_id')))

    return BindingSite(nbr_dict)



def get_polymer_residues(rcsb_id:str, auth_asym_id:str, struct: Structure)->List[Residue]:
    c:Chain = structure[0][auth_asym_id]
    return [*c.get_residues()]

@dataclass(unsafe_hash=True, order=True)
class PolymerRefs              : 
    parent_rcsb_id          : str = field( hash=True ,compare=False)
    auth_asym_ids           : List[str] = field( hash=True ,compare=True )
    rcsb_pdbx_description   : str = field( hash=True ,compare=False)
    entity_poly_seq_length  : int = field( hash=True ,compare=False)
    entity_poly_polymer_type: str = field( hash=True ,compare=False)

def parse_polymer_nbhd(poly:PolymerRefs, structure:Structure):

    STATIC_ROOT  = os.environ.get('STATIC_ROOT')
    for asym_id in poly.auth_asym_ids:
        outfile_json = os.path.join(STATIC_ROOT,poly.parent_rcsb_id.upper(), f'LIGAND_LIKE_{asym_id}.json')
        residues : List[Residue] = get_polymer_residues(poly.parent_rcsb_id, asym_id, structure)
        bs:BindingSite           = get_poly_nbrs(asym_id,residues , structure )
        bs.to_json(outfile_json)


# ! some data classes are due. Nope. don't complicate it.
# Grab chain ---> residues ---> apply neighbor search on  sparsely. done.



def get_liglike_polymers(struct:str)->List[PolymerRefs]:
    """Given an rcsb id, open the profile for the corresponding structure 
    and return references to all polymers marked ligand-like"""

    with open(struct_path(struct,'json'),'r') as infile:
        profile = json.load(infile)
    liglike = []
    for i in [*profile['rnas'], *profile['proteins']]:
        if i['ligand_like'] == True:
            liglike = [*liglike, 
            PolymerRefs(
                i['parent_rcsb_id'],
                i['auth_asym_ids'],
                i['rcsb_pdbx_description'],
                i['entity_poly_seq_length'],
                i['entity_poly_polymer_type'],
            )]

    return liglike

      
if __name__ =="__main__" :
    
    def root_self(rootname:str=''):
        """Returns the rootpath for the project if it's unique in the current folder tree."""
        root=os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
        sys.path.append(root)
    root_self('ribetl')

    from ciftools.neoget import _neoget
    parser      = argparse.ArgumentParser(                                             )
    parser .add_argument ('-s','--structure' , type  = str ,required =True)
    args = parser.parse_args()
    PDBID       = args.structure.upper()

    print("\t\t\t\033[92m * ")
    structure     = open_structure      (PDBID)
    liglike_polys = get_liglike_polymers(PDBID) 
    print(f"\tParsing ligand-like polymers for structure {PDBID}: ")
    if len( liglike_polys ) < 1: print("None found. Exited."); exit(1)
    for polyref in liglike_polys:
        parse_polymer_nbhd(polyref, structure)

    print("\t\t\t * \033[0m")

