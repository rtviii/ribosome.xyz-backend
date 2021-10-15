# paromomycin : 50
# ery : 13
# kirromycin: 7
# viomycin : 12
# Apidaecin : 5
# listerin : 5
# neomycin :5 

from ast import arg, parse
import os,sys
from pprint import pprint
import json
import dotenv
from pymol import cmd
import argparse
import glob
from ciftools.bsite_mixed import BindingSite



# 5JU8 5JTE 4WFN 4V7X 4V7U 3J7z 1YI2 6S0X 6ND6      --- ERY ---> 7aqc
# 5JU8 5JTE 4WFN 4V7X 4V7U 3J7z 1YI2 6S0z 6S0X 6ND6 --- PAR ---> 4lfz
# 6XZB 6XZA 6XZ7 6OF1                               --- DI0 ---> 3j9w 

def get_noms(a:dict):
	_ = []
	for i in a:

		for name in a[i]['nomenclature']:
			if name in _:
				print("DUPLICATE")
			else:
				_ = [*_, [name,len(a[i]['residues'])]  ]
	return _
sys.path.append('/home/rxz/dev/riboxyzbackend/')
dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get("STATIC_ROOT")
nomid       = sys.argv[1]
bsites      = []
for x in glob.glob(STATIC_ROOT + f"/*/POLYMER_*.json"):
	try:
		with open(x, 'rb') as infile:
			data = json.load(infile)	
			bsites.append(BindingSite(data))
	except:
		...
for y in glob.glob(STATIC_ROOT + f"/*/LIGAND_*.json"):
	try:
		with open(y, 'rb') as infile:
			data = json.load(infile)	
			bsites.append(BindingSite(data))
	except:
		...



# have to go through the profiles given that liglike polymer names are not preserved in the filenames
# just start with ligandlike:true --> grab the corresponding profile, look for a given neighbor in nomenclatures
# construct a dict of how many ligands have it match and and at which residues
# plot line by line with strips of ligands of the same kind


# for bsite in 

# pprint(PAR_bsites[0].data)
# pprint(get_noms(PAR_bsites[0].data))

