import json
from pprint import pprint
from Bio.PDB.Chain import Chain
import os
from Bio.PDB.Structure import Structure
import itertools
import itertools
import dotenv
import numpy as np
from Bio.PDB import MMCIF2Dict, MMCIFIO, FastMMCIFParser
try: 
    from rxz_backend.settings import DOTENV_PATH 
except:
    ...

flatten = itertools.chain.from_iterable
import gemmi 
n1      = np.array
import argparse


def get_dict(path:str,)->dict:
	return MMCIF2Dict.MMCIF2Dict(path)

def make_nom_dict(profile):
	nomdict = {}
	for i in [*profile['rnas'], *profile['proteins']]:
		nomdict[i['auth_asym_id']]  = i['nomenclature']
	return nomdict

def struct_path(pdbid: str, pftype: str):
    if pftype == 'cif':
        return os.path.join(os.environ.get("STATIC_ROOT"), pdbid.upper(), f"{pdbid.upper()}.cif")
    if pftype == 'json':
        return os.path.join(os.environ.get("STATIC_ROOT"), pdbid.upper(), f"{pdbid.upper()}.json")

def open_structure(pdbid: str, pftype:str):
    pdbid = pdbid.upper()
    if pftype == 'cif':
        cifpath = struct_path(pdbid,'cif')
        try:
            return FastMMCIFParser(QUIET=True).get_structure(pdbid, cifpath)
        except Exception as e:
            return f"\033[93m Parser Error in structure {pdbid} \033[0m : {e}"

    if pftype == 'json':
        with open(struct_path(pdbid, 'json'), 'rb') as _:
            return json.load(_)

def inject_dict( pdbid:str):

	cifpath       = os.path.join(STATIC_ROOT,pdbid, pdbid+".cif")
	cifmodified   = os.path.join(STATIC_ROOT,pdbid, pdbid+"_modified.cif")
	structprofile = open_structure(pdbid, 'json')

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

	structprofile           = open_structure(pdbid, 'json')
	struct_cif   :Structure = open_structure(pdbid, 'cif' )

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



parser = argparse. ArgumentParser(description='Split structure into constituent polymers and inject new nomencalture into the .cif file')
parser. add_argument ("-s"    , "--structure", type= str , help="RCSB ID of structure to process")
parser. add_argument ("-env"    , "--dotenv_path", type= str , help="Fallback dotenv path. Needed to locate the static files folder")
 
args  = parser.parse_args()
pdbid = args.structure

dotenv.load_dotenv(dotenv_path=args.dotenv_path or DOTENV_PATH)
STATIC_ROOT = os.environ.get('STATIC_ROOT')

if not pdbid:
    print("Provide structure ID with -s arg")
    exit(0)
else:
	inject_dict(pdbid.upper())
	process_chains(pdbid.upper())