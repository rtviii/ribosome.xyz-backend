from itertools import chain
import os
import sys
import dotenv
from gemmi import cif
import gemmi
import ribetl.ciftools.bsite_mixed as bsite

PDBID       = sys.argv[1].upper()
dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get('STATIC_ROOT')

cifpath     = os.path.join(STATIC_ROOT,PDBID, PDBID+".cif")
cifmodified = os.path.join(STATIC_ROOT,PDBID, PDBID+"_modified.cif")

structprofile           = bsite.open_structure(PDBID, 'json')

def make_nom_dict(profile)->dict:
	nomdict = {}
	for i in [*profile['rnas'], *profile['proteins']]:
		for aa in  i['auth_asym_ids']:
			nomdict[aa]  = i['nomenclature']
	return nomdict

doc   = cif.read_file(cifpath)
block = doc.sole_block()
loop  = block.init_loop('_ribosome_nomenclature.', ['entity_poly.pdbx_strand_id', 'polymer_class'])
for i in make_nom_dict(structprofile).items():
	print(i)
	loop.add_row([i[0], cif.quote("unclassified" if len( i[1] )==0 else i[1][0])])

doc.write_file(cifmodified)
print("\033[91m Wrote {} \033[0m".format(cifmodified))

# strand_id = "0"
# nom = "uL22"
# structure   = gemmi.read_structure(cifpath)
# chainA = structure[0]['0']
# print(chainA)
# print(chainA.name)
# chainA.name ='B_uL22'

# s = gemmi.Structure()

# # deletes chains other than B
# for model in structure:
#     for chain in model:
#         if chain.name != 'B':
#             del chain[:]

# s.name = "model_chain_{}_uL22".format(strand_id)
# s.make_mmcif_document().write_file("somechain.cif")