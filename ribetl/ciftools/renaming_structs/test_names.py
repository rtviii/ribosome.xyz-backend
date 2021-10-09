from Bio.PDB.Entity import Entity
from Bio.PDB.Structure import Structure
from Bio.PDB.MMCIFParser import FastMMCIFParser, MMCIFParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.mmcifio import MMCIFIO
import os, sys
import ribetl.ciftools.dict_ops as dict_ops
import ribetl.ciftools.bsite_mixed as bsite
import os, sys
import json
import dotenv

PDBID       = sys.argv[1].upper()
dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get('STATIC_ROOT')

structprofile           = bsite.open_structure(PDBID, 'json')
struct_cif   :Structure = bsite.open_structure(PDBID, 'cif' )

def make_nom_dict(profile):
	nomdict = {}
	for i in [*profile['rnas'], *profile['proteins']]:
		for aa in  i['auth_asym_ids']:
			nomdict[aa]  = i['nomenclature']

	return nomdict




pdbid           = str(sys.argv[1]).upper()
struct_filename = str(sys.argv[2])
parser          = FastMMCIFParser              (QUIET=True                )
io              = MMCIFIO                      (                          )
struct= parser.get_structure('somestruct', '4UG0_STRAND_Lh.cif')

print (struct[0].child_dict['Lh'])

# for chain in struct[0].child_list:

# 	chain:Chain

# 	id = chain.get_id()
# 	if id in chain_classes:
# 		struct[0][id].id = f"{chain_classes[id]}"
		# print("got nomcalss",nomclass_result[0][0])
	# else:
	# 	struct[0][strand_id].id = f"{strand_id}[-]"
	# 	print("got empty")

# io.set_structure(struct)
# result_name = f'{pdbid}_renamed.cif'
# io.save(result_name)
# print("Saved renamed structure to {result_name}")