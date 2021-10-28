from ast import arg, parse
import os,sys
import pprint
import json
import dotenv
from pymol import cmd
import argparse
sys.path.append('/home/rxz/dev/riboxyzbackend')
import ribetl.ciftools.transpose_ligand as transpose_ligand
import ribetl.ciftools.bsite_mixed as bsite_mixed
sys.path.append('../../')


# 5JU8 5JTE 4WFN 4V7X 4V7U 3J7z 1YI2 6S0X 6ND6      --- ERY ---> 7aqc
# 5JU8 5JTE 4WFN 4V7X 4V7U 3J7z 1YI2 6S0z 6S0X 6ND6 --- PAR ---> 4lfz
# 6XZB 6XZA 6XZ7 6OF1                               --- DI0 ---> 3j9w 

dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env'
dotenv.load_dotenv(dotenv_path=dotenv_path)
STATIC_ROOT = os.environ.get("STATIC_ROOT")

print("Environment injected from:", dotenv_path)
print("STATIC ROOT:", STATIC_ROOT)

@cmd.extend
def see_lig_pred(LIG,SRC,TGT):
	cmd.delete('all')

	LIG, SRC, TGT = [_.upper() for _ in [LIG, SRC,TGT]]

	target_path     = os.path.join(STATIC_ROOT, TGT,f"{TGT}.cif")
	prediction_path = os.path.join(STATIC_ROOT, TGT,f"PREDICTION_{LIG}_{SRC}_{TGT}.json")
	cmd.load(target_path)

	if not os.path.isfile(prediction_path):

		with open(os.path.join(STATIC_ROOT,TGT,"{}.json".format(TGT)), 'rb') as infile:
			target_handle = json.load(infile)
		with open(os.path.join(STATIC_ROOT,SRC,"LIGAND_{}.json".format(LIG)), 'rb') as bsitef:
			bsite = json.load(bsitef)
			print("bsite is ", bsite)
		prediction= transpose_ligand.init_transpose_ligand(SRC,TGT,target_handle,bsite_mixed.BindingSite(bsite))
		fname      = f'PREDICTION_{LIG}_{SRC}_{TGT}.json'
		with open(os.path.join(os.getenv('STATIC_ROOT' ),TGT,fname), 'w') as outfile:
			json.dump(prediction,outfile)
			print("\033[093mSucessfully saved prediction {}\033[0m".format(fname))
		# os.system("python3 /home/rxz/dev/riboxyzbackend/ribetl/ciftools/transpose_ligand.py --ligand {} -src {} -tgt {}".format(LIG,SRC,TGT))  

	with  open(prediction_path, 'rb') as infile:
		data = json.load(infile)


	cmd.color('gray40','all')
	for chain in data:
		tgt_auth_asym_id = data[ chain ]['target']['auth_asym_id']
		tgt_resids =  data[ chain ]['target']['tgt_ids']

		for resid in tgt_resids:
			cmd.color('cyan', f'c. {tgt_auth_asym_id} and resi {resid}')
			cmd.show('sticks', f'c. {tgt_auth_asym_id} and resi {resid}')

@cmd.extend
def see_poly_pred(POLY,SRC,TGT):
	cmd.delete('all')

	target_path     = os.path.join(STATIC_ROOT, TGT,f"{TGT}.cif")
	prediction_path = os.path.join(STATIC_ROOT,TGT,'PREDICTION_{}_{}_{}.json'.format(POLY,SRC,TGT))

	cmd.load(target_path)
	if not os.path.isfile(prediction_path):
		os.system("python3 /home/rxz/dev/riboxyzbackend/ribetl/ciftools/transpose_ligand.py --poly {} -src {} -tgt {}".format(POLY,SRC,TGT))  

	with  open(prediction_path, 'rb') as infile:
		data = json.load(infile)


	cmd.color('gray40','all')
	for chain in data:
		tgt_auth_asym_id = data[ chain ]['target']['tgt_auth_asym_id']
		tgt_resids = data[ chain ]['target']['tgt_ids']

		for resid in tgt_resids:
			cmd.color('cyan', f'c. {tgt_auth_asym_id} and resi {resid}')
			cmd.show('sticks', f'c. {tgt_auth_asym_id} and resi {resid}')

@cmd.extend
def see_ligand(LIG,SRC): 
	cmd.delete('all')
	LIG, SRC = [_.upper() for _ in [LIG, SRC]]
	source_path     = os.path.join(STATIC_ROOT, SRC,f"{SRC}.cif")
	lig_path     = os.path.join(STATIC_ROOT, SRC,f"LIGAND_{LIG}.json")

	cmd.load(source_path)
	with  open(lig_path, 'rb') as infile:
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

	source_path = os.path.join(STATIC_ROOT, SRC,f"{SRC}.cif")
	poly_path   = os.path.join(STATIC_ROOT, SRC,f"POLYMER_{POLY}.json")

	cmd.load(source_path)

	with  open(poly_path, 'rb') as infile:
		data = json.load(infile)
		pprint.pprint(data)
	cmd.color('gray40','all')
	cmd.color('cyan',f'c. {POLY}')

	for chain in data:
		'nomenclature'
		resids = [ _['residue_id'] for _ in data[chain]['residues']]
		for i in resids:
			cmd.color('green',f'c. {chain} and resi {i}')

		



		


