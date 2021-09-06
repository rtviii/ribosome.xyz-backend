from ast import arg, parse
import os,sys
import json
from pymol import cmd
import argparse



# 5JU8 5JTE 4WFN 4V7X 4V7U 3J7z 1YI2 6S0X 6ND6      --- ERY ---> 7aqc
# 5JU8 5JTE 4WFN 4V7X 4V7U 3J7z 1YI2 6S0z 6S0X 6ND6 --- PAR ---> 4lfz
# 6XZB 6XZA 6XZ7 6OF1                               --- DI0 ---> 3j9w 



@cmd.extend
def see_prediction(LIG,SRC,TGT):
	cmd.delete('all')
	LIG, SRC, TGT = [_.upper() for _ in [LIG, SRC,TGT]]

	cmd.load('/home/rxz/dev/ribetl/static/{}/{}.cif'.format(TGT,TGT))
	predfile = '/home/rxz/dev/ribetl/static/{}/PREDICTION_{}_{}_{}.json'.format(TGT,LIG,SRC,TGT)

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
def see_ligand(LIG,SRC):
	cmd.delete('all')
	LIG, SRC = [_.upper() for _ in [LIG, SRC]]
	cmd.load('/home/rxz/dev/ribetl/static/{}/{}.cif'.format(SRC,SRC))
	with  open('/home/rxz/dev/ribetl/static/{}/LIGAND_{}.json'.format(SRC,LIG), 'rb') as infile:
		data = json.load(infile)
	cmd.color('gray40','all')
	cmd.color('cyan',f'resn {LIG}')

	for chain in data:
		'nomenclature'
		resids = [ _['residue_id'] for _ in data[chain]['residues']]
		for i in resids:
			cmd.color('green',f'c. {chain} and resi {i}')

		


