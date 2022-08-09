import requests
from Bio.PDB import MMCIF2Dict
from Bio.PDB.MMCIFParser import FastMMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain



def get_struct(pdb_id):
    url = "https://files.rcsb.org/download/" + pdb_id + ".pdb"
    r = requests.get(url)
    return r.text

RCSB_ID  = "3J7Z"
URL      = "https://api.ribosome.xyz/static_files/download_structure/?struct_id=3j7z"
response = requests.get(URL)
with open("{}.cif".format(RCSB_ID), "wb") as outfile:
    outfile.write(response.content)
    cifstruct = FastMMCIFParser(QUIET=True).get_structure(RCSB_ID, outfile)
    print(cifstruct)




