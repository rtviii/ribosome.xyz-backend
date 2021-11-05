from ast import arg, parse
import os,sys
import pprint
import json
import dotenv
from pymol import cmd
import argparse
sys.path.append('/home/rxz/dev/riboxyzbackend')
sys.path.append('../../')

import ribetl.ciftools.transpose_ligand as transpose_ligand
from ribetl.ciftools.super import ranged_super, ranged_super_class
import ribetl.ciftools.bsite_mixed as bsite_mixed

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

		
@cmd.extend
def rsuper(src_struct,tgt_struct,rstart,rend, poly_class,):
	cmd.delete('all')

	[ c1 ,r1, c2,r2 ] = ranged_super_class(src_struct, tgt_struct,( int(rstart), int(rend )), poly_class)

	n1 = c1.split("/")[-1].split('.')[0]
	n2 = c2.split("/")[-1].split('.')[0]

	src_snip = "{}_{}".format(src_struct, poly_class)
	tgt_snip = "{}_{}".format(tgt_struct, poly_class)

	cmd.load(c1)
	cmd.select("resi {}-{} and m. {} ".format(*r1,n1))
	cmd.create(src_snip,"sele")
	cmd.color('cyan',src_snip)
	cmd.delete(n1)

	cmd.load(c2)
	cmd.select("resi {}-{} and m. {} ".format(*r2,n2))
	cmd.create(tgt_snip,"sele")
	cmd.color('wheat',tgt_snip)
	cmd.delete(n2)

	# print("Aligned with CEALIGN.")
	# cmd.cealign("{}_{}".format(src_struct, poly_class),"{}_{}".format(tgt_struct, poly_class))

	print("Aligned with SUPER.")
	cmd.super("{}_{}".format(src_struct, poly_class),"{}_{}".format(tgt_struct, poly_class))

	# print("Aligned with ALIGN.")
	# cmd.align(src_snip,tgt_snip)


