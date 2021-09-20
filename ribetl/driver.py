import enum
from dotenv import load_dotenv
import os
import ciftools.dict_ops as dops
import itertools 
import Bio.PDB.mmcifio as mmcifio
from ciftools.neoget import get_nom_cmap
import itertools
import logging

import numpy as np
flatten = itertools.chain.from_iterable
n1      = np.array
STATIC_ROOT = os  .getenv   ('STATIC_ROOT')

io = mmcifio.MMCIFIO()


#! I guess we'll have to document the types of exceptions that mmcifio throws. Stoked :(
#? Relevant fields:
# usually fall through
# _atom_site.auth_comp_id
# _atom_site.auth_asym_id
# _atom_site.auth_atom_id
# _atom_site.pdbx_PDB_model_num

# -----------------------------
# _entity_poly.pdbx_strand_id
#******* we are adding thesE:


PDBID   = '3j7z'.upper()
struct  = dops.get_dict(path=f'{PDBID}.cif')
cmap    = get_nom_cmap(PDBID)
strands = dops._get_strand_ids(struct)

for  strand in strands:
	if strand in cmap:
		...
	else:
		cmap.update({strand: None})
		print("Updated cmap. (Strand {} wasnt found)".format(strand))


del struct['_atom_site.auth_comp_id']
del struct['_atom_site.auth_asym_id']
del struct['_atom_site.auth_atom_id']
del struct['_atom_site.pdbx_PDB_model_num']

remark_ids  = list(cmap.keys())
remark_text = list(cmap.values())
for i, v in enumerate( remark_text ):
	if v == None:
		remark_text[i] = "Unclassified"



print(remark_text)

struct['_pdbx_database_remark.id']   = remark_ids
struct['_pdbx_database_remark.text'] = remark_text

print(remark_ids)
print(remark_text)
struct['_ubc_ribosomal_protein_name.entry_id']  = remark_ids
struct['_ubc_ribosomal_protein_name.entity_id'] = remark_text

# struct['_ubc_ribosomal_protein_name.name'] =  struct['_poly_entity.pdbx_strand_id']


print(struct.keys())
io.set_dict(struct)
io.save("here.cif")
# try:
# 	io.save("augmented.cif")
# except Exception as e:
# 	# logging.error(f'Struct {PDBID}', e)
# 	print(e)








