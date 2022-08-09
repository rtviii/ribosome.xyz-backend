from io import FileIO
from re import template
import requests
from Bio.PDB import MMCIF2Dict
from Bio.PDB.MMCIFParser import FastMMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
import gemmi
import tempfile

RCSB_ID  = "3J7Z"
URL      = "https://api.ribosome.xyz/static_files/download_structure/?struct_id=4v7s"
response = requests.get(URL)

# with open("{}.cif".format(RCSB_ID), "wb") as outfile:
#     x = outfile.write(response.content)
#     cifstruct = FastMMCIFParser(QUIET=True).get_structure(RCSB_ID, x)
#     print(cifstruct)

with tempfile.NamedTemporaryFile() as outfile:
    outfile.write(response.content)
    outfile.flush()
    filename = outfile.name
    # cifstruct = FastMMCIFParser(QUIET=True).get_structure(RCSB_ID, filename)
    # y = MMCIF2Dict.MMCIF2Dict(filename)
    # doc   = gemmi.cif.read_file(filename)

    cif_block = gemmi.cif.read(filename)[0]
    structure = gemmi.make_structure_from_block(cif_block)
    for entity in structure.entities:
        print(dir(entity))
        print(entity.name)
        print(entity.polymer_type)
        print(entity.full_sequence)
        
    # for chain in structure[0]:
    #     print(chain.one_letter_code)





