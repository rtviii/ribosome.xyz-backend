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

def cli():
    parser = ArgumentParser(description='Residue-wide homology analysis. Primary focus is the exit tunnel')
    parser.add_argument('--verbose', help='To enable atomwise logging. Not useful in any way.')
    parser.add_argument('-p','--path',help="Path to file(Relative?)")
    parser.add_argument('-t','--target-protein', help="Protein to operate inside the chosen strucutre")
    parser.add_argument('-r','--radius', help='clustering radius', type=float)
    return parser.parse_args()

def euc_dist(n1, n2):
    # a scipy function ideally
    return math.sqrt((n1[0]-n2[0] )**2 +(n1[1]-n2[1] )**2 + (n1[2]-n2[2])**2 )

def main(): 
#   uL4   : Chain = struct[0].child_dict['LC']
    struct: Structure = FastMMCIFParser(QUIET=True).get_structure(structure_id='4UG0',filename='4ug0.cif')
    uL22  : Chain     = struct[0].child_dict['LP']
    tunnelcoords = [ np.asfarray(atom) for atom in tunnel_atoms]
    adjresidues  = get_adjacent_residues(uL22, tunnelcoords,10)
    [ print(res) for res in adjresidues]

    
    
def get_adjacent_residues(chain:Chain,hull, radius:int):
    # @param hull -- the coordinate set for the tunnel as a convex hull
    # @param radius -- angstrom radius to treshold adjacency
    adjacent_residues =set()
    for resn in [ *chain.get_residues()]:
        for atom in hull:
            if euc_dist(getAlphaCarbonResidue(resn).get_coord(), atom) < radius:
                adjacent_residues.add(resn)
    
    return adjacent_residues


def getAlphaCarbonResidue(res:Residue):
    atoms = [ *res.get_atoms() ]
    alphaCarbons =  [ *filter(lambda atom: atom.get_fullname() == "CA", atoms) ]
    #If residue doesnt have an alpha carbon, return any other 
    return alphaCarbons[0] if len(alphaCarbons)> 0 else atoms[0]


def getResidueById(chainObj:Chain, id:int):
    chains = chainObj.child_dict
    return next(chains[key] for key in chains if id in key)


def driver():   
    args = cli()
    main()

if __name__=="__main__":
    print(f"""⋯⋯⋅⋅⋅⋅┈⋅⋄⋅⋅⋅⋄⋯⋅⋅⋅⋅┈⋅⋄⋅⋅⋅⋄⋅⋅⋅⋅┈⋅⋄⋅⋅⋅⋄⋯⋯⋅⋅⋅⋅┈⋅
    A template module for ribosome.xyz.
    Running {__file__} as standalone module.
    ⋯⋯⋅⋅⋅⋅┈⋅⋄⋅⋅⋅⋄⋯⋅⋅⋅⋅┈⋅⋄⋅⋅⋅⋄⋅⋅⋅⋅┈⋅⋄⋅⋅⋅⋄⋯⋯⋅⋅⋅⋅┈⋅
    \n\n
    """)
    driver()