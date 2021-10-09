from ast import arg, parse
import os,sys
import pprint
import json

import dotenv
from pymol import cmd
import argparse

sys.path.append('/home/rxz/dev/riboxyzbackend/')

# 5JU8 5JTE 4WFN 4V7X 4V7U 3J7z 1YI2 6S0X 6ND6      --- ERY ---> 7aqc
# 5JU8 5JTE 4WFN 4V7X 4V7U 3J7z 1YI2 6S0z 6S0X 6ND6 --- PAR ---> 4lfz
# 6XZB 6XZA 6XZ7 6OF1                               --- DI0 ---> 3j9w 

dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT =os.environ.get("STATIC_ROOT")


print("STATIC ROOT:", STATIC_ROOT)

@cmd.extend
def see_lig_pred(LIG,SRC,TGT):
	cmd.delete('all')
	LIG, SRC, TGT = [_.upper() for _ in [LIG, SRC,TGT]]

	cmd.load('/home/rxz/dev/riboxyzbackend/ribetl/static/{}/{}.cif'.format(TGT,TGT))

	predfile = '/home/rxz/dev/riboxyzbackend/ribetl/static/{}/PREDICTION_{}_{}_{}.json'.format(TGT,LIG,SRC,TGT)
	with  open(predfile, 'rb') as infile:
		data = json.load(infile)
	cmd.color('gray40','all')
	for chain in data:
		tgt_strand = data[ chain ]['target']['strand']
		tgt_resids =  data[ chain ]['target']['tgt_ids']

		for resid in tgt_resids:
			cmd.color('cyan', f'c. {tgt_strand} and resi {resid}')
			cmd.show('sticks', f'c. {tgt_strand} and resi {resid}')

@cmd.extend
def see_poly_pred(POLY,SRC,TGT):
	cmd.delete('all')
	POLY, SRC, TGT = [_.upper() for _ in [POLY, SRC,TGT]]

	cif_path        = os.path.join(STATIC_ROOT,TGT,'{}.cif'.format(TGT))
	prediction_path = os.path.join(STATIC_ROOT,TGT,'PREDICTION_{}_{}_{}.json'.format(POLY,SRC,TGT))

	cmd.load(cif_path)
	if not os.path.isfile(prediction_path):
		os.system("python3 /home/rxz/dev/riboxyzbackend/ribetl/ciftools/transpose_ligand.py --poly {} -src {} -tgt {}".format(POLY,SRC,TGT))  

	with  open(prediction_path, 'rb') as infile:
		data = json.load(infile)

	cmd.color('gray40','all')
	for chain in data:
		tgt_strand = data[ chain ]['target']['strand']
		tgt_resids = data[ chain ]['target']['tgt_ids']

		for resid in tgt_resids:
			cmd.color('cyan', f'c. {tgt_strand} and resi {resid}')
			cmd.show('sticks', f'c. {tgt_strand} and resi {resid}')

@cmd.extend
def see_ligand(LIG,SRC):
	cmd.delete('all')
	LIG, SRC = [_.upper() for _ in [LIG, SRC]]

	cmd.load('/home/rxz/dev/riboxyzbackend/ribetl/static/{}/{}.cif'.format(SRC,SRC))
	with  open('/home/rxz/dev/riboxyzbackend/ribetl/static/{}/LIGAND_{}.json'.format(SRC,LIG), 'rb') as infile:
		data = json.load(infile)
	cmd.color('gray40','all')
	cmd.color('cyan',f'resn {LIG}')

	for chain in data:
		'nomenclature'
		resids = [ _['residue_id'] for _ in data[chain]['residues']]
		for i in resids:
			cmd.color('green',f'c. {chain} and resi {i}')

@cmd.extend
def see_poly(POLY,SRC):
	cmd.delete('all')
	POLY, SRC = [_.upper() for _ in [POLY, SRC]]
	cmd.load('/home/rxz/dev/riboxyzbackend/ribetl/static/{}/{}.cif'.format(SRC,SRC))
	with  open('/home/rxz/dev/riboxyzbackend/ribetl/static/{}/POLYMER_{}.json'.format(SRC,POLY), 'rb') as infile:
		data = json.load(infile)
		pprint.pprint(data)
	cmd.color('gray40','all')
	cmd.color('cyan',f'c. {POLY}')

	for chain in data:
		'nomenclature'
		resids = [ _['residue_id'] for _ in data[chain]['residues']]
		for i in resids:
			cmd.color('green',f'c. {chain} and resi {i}')

		



		


