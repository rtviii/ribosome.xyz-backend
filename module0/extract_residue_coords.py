from pymol import cmd
import numpy as np
import math
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import MMCIF2Dict
from Bio.PDB.Structure import Structure
from Bio.PDB import FastMMCIFParser
from Bio.PDB import StructureBuilder
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Polypeptide import PPBuilder,Polypeptide
from argparse import ArgumentParser
from tunnel import tunnel_atoms


def getAlphaCarbonResidue(res:Residue):
    atoms = [ *res.get_atoms() ]
    alphaCarbons =  [ *filter(lambda atom: atom.get_fullname() == "CA", atoms) ]
    #If residue doesnt have an alpha carbon, return any other 
    return alphaCarbons[0] if len(alphaCarbons)> 0 else atoms[0]


def main(): 
    struct: Structure = FastMMCIFParser(QUIET=True).get_structure(structure_id='5NWY',filename='5nwy.cif')
    chain_s  : Chain     = struct[0].child_dict['s']    
    residues = [ *chain_s.get_residues() ]
    for residue in residues:
        j = getAlphaCarbonResidue(residue)
        print(f" Atom {j.get_full_id()} at {j.get_coord()}")

main()


