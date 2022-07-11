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
from ribetl.ciftools.bsite_mixed import BindingSite



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
STATIC_ROOT =os.environ.get("STATIC_ROOT")

chemid= sys.argv[1].upper()
PAR_bsites = []
for x in glob.glob(STATIC_ROOT + f"/*/LIGAND_{chemid}*.json"):
	with open(x, 'rb') as infile:
		data = json.load(infile)	
		PAR_bsites.append(BindingSite(data))

# pprint(PAR_bsites[0].data)
# pprint(get_noms(PAR_bsites[0].data))

count ={
}
for bsite in PAR_bsites:
	for name,matchlen in get_noms(bsite.data):
		print(name,matchlen)
		if name not in count:
			count[name] = ( 1, [ matchlen ] )
		else:
			count[name] = (  count[name][0]+1 , [*count[name][1] ,matchlen])

print(count)
print(f"{chemid} interface members:")
for key,value in sorted(count.items(), key=lambda item: item[1][0]):
    print("{}:         {}".format(key, value))
print(f"Total number of binding sites for {chemid}:", len(PAR_bsites))

