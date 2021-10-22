from pprint import pprint
import struct
from Bio.PDB.Chain import Chain
import os
import sys
from Bio.PDB.Structure import Structure
import itertools
import itertools
import dotenv
import numpy as np
from Bio.PDB import MMCIF2Dict, MMCIFIO
flatten = itertools.chain.from_iterable
import ribetl.ciftools.bsite_mixed as bsite
import gemmi 
n1      = np.array

dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get('STATIC_ROOT')

def get_dict(path:str,)->dict:
	return MMCIF2Dict.MMCIF2Dict(path)

def make_nom_dict(profile):
	nomdict = {}
	for i in [*profile['rnas'], *profile['proteins']]:
		nomdict[i['auth_asym_id']]  = i['nomenclature']
	return nomdict



def inject_dict( pdbid:str):

	cifpath       = os.path.join(STATIC_ROOT,pdbid, pdbid+".cif")
	cifmodified   = os.path.join(STATIC_ROOT,pdbid, pdbid+"_modified.cif")
	structprofile = bsite.open_structure(pdbid, 'json')

	doc   = gemmi.cif.read_file(cifpath)
	block = doc.sole_block()
	loop  = block.init_loop('_ribosome_nomenclature.', ['entity_poly.pdbx_strand_id', 'polymer_class'])

	nomd = make_nom_dict(structprofile)
	for i in make_nom_dict(structprofile).items():
		loop.add_row([i[0], gemmi.cif.quote("unclassified" if len( i[1] )==0 else i[1][0])])

	doc.write_file(cifmodified)
	print("\033[91m Wrote {} \033[0m".format(cifmodified))


def process_chains( pdbid:str):

	io          = MMCIFIO()
	STATIC_ROOT = os.environ.get('STATIC_ROOT')

	if not os.path.exists(os.path.join(STATIC_ROOT,pdbid,'CHAINS')):
		os.mkdir(os.path.join(STATIC_ROOT,pdbid,'CHAINS'))

	structprofile           = bsite.open_structure(pdbid, 'json')
	struct_cif   :Structure = bsite.open_structure(pdbid, 'cif' )

	model = gemmi.read_structure(os.path.join(STATIC_ROOT, pdbid, f'{pdbid}.cif'))[0]
	nomd  = make_nom_dict(structprofile)

	for chain in struct_cif[0].child_list:
		io.set_structure(chain)
		destination = os.path.join(STATIC_ROOT,pdbid,'CHAINS', '{}_STRAND_{}.cif'.format(pdbid,chain.get_id()) ) 
		io.save(destination)
		cdict = get_dict(destination)
		if len(nomd[chain.get_id()]) < 1:
			cdict['data_']=f'{pdbid}_{chain.get_id()}'
		else:
			cdict['data_']=f'{pdbid}_{chain.get_id()}_{nomd[chain.get_id()][0]}'
		io.set_dict(cdict)
		print("Saved ", destination)
		io.save(destination)


pdbid = sys.argv[1]	

inject_dict(pdbid.upper())
process_chains(pdbid.upper())