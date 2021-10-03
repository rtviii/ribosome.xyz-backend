import pdb
from pprint import pprint
from turtle import update
from Bio.PDB.Entity import Entity
from Bio.PDB.Structure import Structure
from Bio.PDB.MMCIFParser import FastMMCIFParser, MMCIFParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.mmcifio import MMCIFIO
import os, sys
import json
from asyncio import run
from typing import TypedDict

import dotenv




pdbid       = sys.argv[1].upper()


dotenv.load_dotenv(dotenv_path='./../../../rxz_backend/.env')
STATIC_ROOT = os.environ.get('STATIC_ROOT')

jsonprofile = os.path.join(STATIC_ROOT,pdbid,f'{pdbid.upper()}.json')
cifprofile  = os.path.join(STATIC_ROOT,pdbid,f'{pdbid.upper()}.cif')



# Grab the json profile, iterate throuhg rnas and proteins extracting chain_ids and their nomenclature
def chain_nomenclature_dict(path:str):
	with open(path,'r') as infile:
		profile = json.load(infile)
	chain_nom_dict = {}
	
	for rp in profile['proteins']:
		chain_nom_dict.update({rp['entity_poly_strand_id']:rp['nomenclature']})

	for rna in profile['rnas']:
		chain_nom_dict.update({rna['entity_poly_strand_id']:rna['nomenclature']})
	return chain_nom_dict

def flatten(x):
    if isinstance(x, list):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]



parser          = FastMMCIFParser              (QUIET=True                )
io              = MMCIFIO                      (                          )
struct:Structure= parser         .get_structure(pdbid, cifprofile)
CEND = f'\033[93m got sturct {pdbid}\033[0m'
nomdict = chain_nomenclature_dict(jsonprofile)
# Keep track of whether certain names have already appeared in the structure to avoid collisions
nomenclature_count = {}
for chain in struct[0].child_list:
	chain:Chain
	id = chain.get_id()
	if id in nomdict:
		if len(nomdict[id] ) >0:
			nom_class = nomdict[id][0] 
			if nom_class in nomenclature_count:
				assigned= f'{nom_class}_{nomenclature_count[nom_class]}'
				nomenclature_count[nom_class]+=1
			else:
				assigned= f'{nom_class}'
				nomenclature_count.update({nom_class:1})
			struct[0][id].id = assigned

io.set_structure(struct)
io.save(cifprofile)
print(f"Saved renamed structure to {cifprofile}")