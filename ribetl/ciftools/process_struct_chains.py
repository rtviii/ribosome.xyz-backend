from pprint import pprint
import struct
from termios import CIBAUD
from Bio.PDB.Chain import Chain
import os
import sys
from Bio.PDB.Structure import Structure
import itertools
import itertools
import dotenv
import numpy as np
import render_binding_site as bsite
from Bio.PDB import MMCIF2Dict, MMCIFIO
import dict_ops as dops
flatten = itertools.chain.from_iterable
n1      = np.array

if __name__ == "__main__":
	dotenv.load_dotenv(dotenv_path='/home/rxz/dev/ribetl/.env')

PDBID  = sys.argv[1].upper()
io     = MMCIFIO()
struct = dops.get_dict('../3J7Z.cif')

if not os.path.exists(os.path.join(os.getenv('STATIC_ROOT'),PDBID,'CHAINS')):
	os.mkdir(os.path.join(os.getenv('STATIC_ROOT'),PDBID,'CHAINS'))

struct = bsite.openStructutre(PDBID)




for chain in struct[0].child_list:
	chain:Chain

	io.set_structure(chain)
	destination = os.path.join(os.getenv( 'STATIC_ROOT' ),PDBID,'CHAINS', '{}_STRAND_{}.cif'.format(PDBID,chain.get_id()) )
	io.save(destination)

	cdict = dops.get_dict(destination)
	# print(cdict.keys())
	cdict['data_']=f'3J7Z_{chain.get_id()}'
	print(cdict['data_'])
	io.set_dict(cdict)
	io.save(destination)

	# id = chain.get_id()
	# if id in cmap.keys():
	# 	if cmap[id] in struct[0].child_dict:
	# 		continue
	# 	struct[0][id].id = f"{cmap[id]}"
	# io          = MMCIFIO()
	# io.set_structure(chain)
	# io.save(destination)





