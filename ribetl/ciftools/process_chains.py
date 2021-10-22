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

def get_dict(path:str,)->dict:
	return MMCIF2Dict.MMCIF2Dict(path)
def make_nom_dict(profile):
	nomdict = {}
	for i in [*profile['rnas'], *profile['proteins']]:
		for aa in  i['auth_asym_ids']:
			nomdict[aa]  = i['nomenclature']
	return nomdict


if __name__ == "__main__":

	dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
	STATIC_ROOT = os.environ.get("STATIC_ROOT")
	PDBID       = sys.argv[1].upper()
	io          = MMCIFIO()

	if not os.path.exists(os.path.join(os.environ.get('STATIC_ROOT'),PDBID,'CHAINS')):
		os.mkdir(os.path.join(os.environ.get('STATIC_ROOT'),PDBID,'CHAINS'))


	structprofile           = bsite.open_structure(PDBID, 'json')
	struct_cif   :Structure = bsite.open_structure(PDBID, 'cif' )

	model = gemmi.read_structure(os.path.join(STATIC_ROOT, PDBID, f'{PDBID}.cif'))[0]
	nomd  = make_nom_dict(structprofile)

	for chain in struct_cif[0].child_list:
		io.set_structure(chain)
		destination = os.path.join(os.environ.get('STATIC_ROOT'),PDBID,'CHAINS', '{}_STRAND_{}.cif'.format(PDBID,chain.get_id()) ) if len(nomd[chain.get_id()]) <1 else  os.path.join(os.environ.get( 'STATIC_ROOT' ),PDBID,'CHAINS', '{}_STRAND_{}_{}.cif'.format(PDBID,chain.get_id(), nomd[chain.get_id()][0]) )
		io.save(destination)
		cdict = get_dict(destination)

		if len(nomd[chain.get_id()]) < 1:
			cdict['data_']=f'{PDBID}_{chain.get_id()}'
		else:
			cdict['data_']=f'{PDBID}_{chain.get_id()}_{nomd[chain.get_id()][0]}'
		io.set_dict(cdict)
		print("Saved ", destination)
		io.save(destination)