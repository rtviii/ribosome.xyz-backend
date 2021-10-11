from Bio.PDB.Structure import Structure
from Bio.PDB.MMCIFParser import FastMMCIFParser, MMCIFParser
from Bio.PDB.Chain import Chain
from Bio.PDB.mmcifio import MMCIFIO
import ribetl.ciftools.bsite_mixed as bsite
import os, sys
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

def flatten(x):
    if isinstance(x, list):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]



parser          = FastMMCIFParser              (QUIET=True                )
io              = MMCIFIO                      (                          )
cifpath = os.path.join(STATIC_ROOT, PDBID, "{}.cif".format(PDBID))

CEND = f'\033[93m got sturct {PDBID}\033[0m'
nomdict = make_nom_dict(structprofile)
# Keep track of whether certain names have already appeared in the structure to avoid collisions
nomenclature_count = {}
for chain in struct_cif[0].child_list:
	chain:Chain
	id = chain.get_id()
	if id in nomdict:
		if len(nomdict[id]) >0:
			nom_class = nomdict[id][0] 
			if nom_class in nomenclature_count:
				assigned= f'{nom_class}_{nomenclature_count[nom_class]}'
				nomenclature_count[nom_class]+=1
			else:
				assigned= f'{nom_class}'
				nomenclature_count.update({nom_class:1})
			struct_cif[0][id].id = assigned

io.set_structure(struct_cif)
io.save(cifpath)
print(f"Saved renamed structure to {cifpath}")