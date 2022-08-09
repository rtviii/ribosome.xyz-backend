from io import FileIO
from re import template
import requests
from Bio.PDB import MMCIF2Dict
from Bio.PDB.MMCIFParser import FastMMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.Align import MultipleSeqAlignment
from Bio.PDB.Chain import Chain
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
import gemmi
import tempfile

# RCSB_ID  = "3J7Z"
# URL      = "https://api.ribosome.xyz/static_files/download_structure/?struct_id=4v7s"
# response = requests.get(URL)

# with open("{}.cif".format(RCSB_ID), "wb") as outfile:
#     outfile.write(response.content)

cif = gemmi.cif.read_file('3J7Z.cif')
# doc =cif.document()
block = cif.sole_block()

# print(block.find_value('_ribosome_nomenclature.entity_poly.pdbx_strand_id.AA'))
# for ( key,val )  in zip( block.find_loop('_ribosome_nomenclature.entity_poly.pdbx_strand_id'), block.find_loop('_ribosome_nomenclature.polymer_class')):
#     print(key,val)
# for ( chain_id,one_letter_code )  in zip( 
#                                          block.find_loop('_entity_poly.pdbx_strand_id'),
#                                          block.find_loop('_entity_poly.pdbx_seq_one_letter_code')
#                                         ):
#     print(chain_id,one_letter_code)



s1= "AAAAASDF"
s2= "TAAAAASDF"
s3= "TAAGGAASDF"
print(SeqRecord))
# ms= MultipleSeqAlignment([s1,s2,s3])
# print(ms)

# for i in block:
#     print(block)



# with tempfile.NamedTemporaryFile() as outfile:
#     outfile.write(response.content)
#     outfile.flush()
#     filename = outfile.name
#     # cifstruct = FastMMCIFParser(QUIET=True).get_structure(RCSB_ID, filename)
#     # y = MMCIF2Dict.MMCIF2Dict(filename)
#     # doc   = gemmi.cif.read_file(filename)

#     cif_block = gemmi.cif.read(filename)[0]
#     structure = gemmi.make_structure_from_block(cif_block)
#     # for entity in structure.entities:
#     #     print(dir(entity))
#     #     print(entity.name)
#     #     print(entity.polymer_type)
#     #     print(entity.full_sequence)
        
#     for chain in structure:
#         print(chain)





