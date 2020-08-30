from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import MMCIF2Dict
from Bio.PDB.Structure import Structure
from Bio.PDB import FastMMCIFParser
from Bio.PDB import StructureBuilder
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Polypeptide import PPBuilder,Polypeptide
from Bio.PDB import PDBIO



pdbid     = '3J9M'
filepath  = f'./{pdbid}.cif'
cifparser = FastMMCIFParser(QUIET=True)
try:
    with open(filepath) as infile:
        structure:Structure = cifparser.get_structure(pdbid, infile)[0]
except:
        print("Failed to open {}".format(filepath))

mychain:Chain = structure['D']

writer = PDBIO()
writer.set_structure(mychain)
print(f"Got chain {mychain}")
# with open('pdbchain.pdb', "w+") as out:
writer.save('newchain.pdb')


