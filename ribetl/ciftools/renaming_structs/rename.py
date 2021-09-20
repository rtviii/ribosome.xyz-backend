from pprint import pprint
from Bio.PDB.Entity import Entity
from Bio.PDB.Structure import Structure
from Bio.PDB.MMCIFParser import FastMMCIFParser, MMCIFParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.mmcifio import MMCIFIO
import os, sys
from asyncio import run
from neo4jreq import get_classes, get_rna
from typing import TypedDict


def flatten(x):
    if isinstance(x, list):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

pdbid         = str( sys.argv[1] ).upper()
chain_classes = dict([ *map(lambda x: [x[1],x[0]],get_classes(pdbid)) ])
rnas          = flatten( get_rna(pdbid) )

pprint(chain_classes)

pprint(rnas)



# @d rcsb_description to regex through
def rna_class(d:str)->str:
	x = d.lower()
	...

	#! TODO: Conditions for each class right now <<<<<< hacky
	#? Incorporate RNA Nomenclature into the schema
	...
	return 'other'

class Rna(TypedDict):
	pdbx_description:str
	strand_id       :str

# for x in rnas:
# 	x:Rna
# 	print(x)
# 	sort_rna_description(x['pdbx_description'])


parser          = FastMMCIFParser              (QUIET=True                )
io              = MMCIFIO                      (                          )
pdbid           = pdbid           .upper       (                          )
struct:Structure= parser         .get_structure(f'{pdbid}', f'{pdbid}.cif')

for chain in struct[0].child_list:

	chain:Chain

	id = chain.get_id()
	if id in chain_classes:
		struct[0][id].id = f"{chain_classes[id]}"

io.set_structure(struct)
result_name = f'{pdbid}_renamed.cif'
io.save(result_name)
print(f"Saved renamed structure to {result_name}")