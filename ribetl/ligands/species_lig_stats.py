from ctypes import Union
from fnmatch import translate
import os
import re
import sys
from pprint import pprint
import json
from typing import List
from Bio import Align
from Bio.SeqRecord import SeqRecord
import dotenv
from pkg_resources import FileMetadata
from pymol import cmd
import argparse
from typing import TypedDict, Union
import glob
from ribetl.ciftools.bsite_mixed import BindingSite
from Bio.Align import MultipleSeqAlignment
from enum import Enum
import numpy as np
from Bio import pairwise2
import matplotlib.pyplot as plt

def lig_cat(description:str)->str:

	if "[(" in description.lower() :
		return "Mixed"

	if len(re.findall(r"(\w*(?<!(cha|pro|dom|str|pla|tox))in\b|(\b\w*zyme\b))", description.lower()))> 0:

		return "Antibiotics"

	if len(re.findall(r"(factor)", description.lower()))> 0:
        # elongation
        # release
        # recycling
        # translation
        # maturation 
		return "Factors"

	if "mrna" in description.lower() or "messenger" in description.lower():

		return "mRNA"

	if "trna" in description.lower() or "transfer" in description.lower():
        # P-
        # E-
        # A-
        # fmet
        # phe


		return "tRNA"

	return "Mixed"

def get_noms(a: dict):
    _ = []
    for i in a:
        for name in a[i]['nomenclature']:
            if name in _:
                print("DUPLICATE")
            else:
                _ = [*_, [name, len(a[i]['residues'])]]
    return _

sys.path.append('/home/rxz/dev/riboxyzbackend/')
dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')

STATIC_ROOT = os.environ.get("STATIC_ROOT")
bsites      = []
profiles    = []

SPECIES = int(sys.argv[1])

class LigandCategory(Enum):
    tRNA        = 'tRNA'
    mRNA        = 'mRNA'
    Factor      = 'Factor'
    Antibiotics = 'Antibiotics'
    Mixed       = 'Mixed'

class LigandLike(TypedDict): 
      description          : str
      parent               : str
      chain                : str
      src_organism_ids     : List[int]
      path                 : str
      category             : LigandCategory

for _ in os.listdir(STATIC_ROOT):
    if len(_) > 4:
        continue
    else:
        profiles = [*profiles,
                    os.path.join(STATIC_ROOT, _, "{}.json".format(_))]




def profile_ligandlikes_by_species(profile_path: str, sought_spec: int)-> List[LigandLike]:
    profile = {}  
    _       = []
    with open(profile_path, 'rb') as infile:

        profile     = json.load(infile)
        polymers    = []
        nonpolymers = []

        if sought_spec not in profile['src_organism_ids']:
            return []

        if profile['proteins'] != None:
            polymers = [*polymers, *profile['proteins']]

        if profile['rnas'] != None:
            polymers = [*polymers, *profile['rnas']]

        try:
            nonpolymers = [*profile['ligands']]

        except:
            ...

    for poly in polymers:
        if bool(poly['ligand_like']) == True:
            _.append({
                'description'     : poly['rcsb_pdbx_description'],
                'category'        : lig_cat(poly['rcsb_pdbx_description']),
                'parent'          : profile['rcsb_id'],
                'chain'           : poly['auth_asym_id'],
                'src_organism_ids': poly['src_organism_ids'],
                'path'            : os.path.join(STATIC_ROOT, profile['rcsb_id'].upper(), "POLYMER_{}.json".format(poly['auth_asym_id']))})

    for np in nonpolymers:
        if not "ion" in np['chemicalName'].lower():
            _.append({
                'description': np['chemicalName'],
                'category'   : lig_cat(np['chemicalName']),
                'parent'     : profile['rcsb_id'],
                'path'       : os.path.join(STATIC_ROOT, profile['rcsb_id'].upper(), "LIGAND_{}.json".format(np['chemicalId']))})
    return _



# p_counts= {
#     "antibiotic"      : 0,
#     "mrna_trna_factor": 0,
#     "both"            : 0
# }

category_counts= {
    'mRNA'       : 0,
    'tRNA'       : 0,
    'Factors'    : 0,
    'Antibiotics': 0,
}


# category_counts= {
#     'mRNA'       : 0,
#     'tRNA'       : 0,
#     'Factors'    : 0,
#     'Antibiotics': 0,
# }


for ppath in profiles: 
    ligandlike = profile_ligandlikes_by_species(ppath, SPECIES)
    counts =  [ 0,0 ]
    
    if len( ligandlike) > 0 :
        for l in ligandlike:
            path     = l['path']
            desc     = l['description']
            category = l['category']

            if category == 'Mixed':
                continue

            # if category == 'Factors':
            #     print("Factors: ", desc)
            if category == 'Antibiotics':
                print("Antibiotics: ", desc)



# labels  = 'Antibiotics', 'mRNA/tRNA/Factor', 'Both'
# sizes   = [p_counts['antibiotic'], p_counts['both'], p_counts['mrna_trna_factor']]
# explode = (0, 0.1, 0)
# # +
# # labels  = 'Antibiotics', 'mRNA', 'tRNA', 'Factors'
# # sizes   = [ category_counts['Antibiotics'], category_counts['mRNA'], category_counts['tRNA'],category_counts['Factors'] ]
# # explode = (0, 0.1, 0,0)  # only "explode" the 2nd slice (i.e. 'Hogs')


# fig1, ax1 = plt.subplots()
# ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
#         shadow=True, startangle=90)
# ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

# plt.show()
# !------------------------------------------------------------------------------#
# fig, ax = plt.subplots()

# size = 0.3

# vals = np.array([[60., 32.], [37., 300], [29., 10.]])
# cmap = plt.get_cmap("tab20c")

# outer_colors = cmap(np.arange(3)*4)
# inner_colors = cmap(np.array([1, 2, 5, 6, 9, 10]))

# ax.pie(vals.sum    (axis=1), radius=1-size     , colors=outer_colors,wedgeprops=dict(width=size, edgecolor='w'))
# ax.pie(vals.flatten(      ), radius=1, colors=inner_colors,wedgeprops=dict(width=size, edgecolor='w'))

# ax.set(aspect="equal", title='Pie plot with `ax.pie`')
# plt.show()
# !------------------------------------------------------------------------------#


#? taxids of interest
# Curate antibiotics manually for the 4 species; trna/factor categories are a little easier to curate; mrna is one solid category
# 
# human 9606, antibiotics: 

# {
#     'streptomycin': 0,
#     'bystin'      : 0,
#     'cadherin'    : 0,
#     'listerin'    : 0,
#     'hygromycin'  : 0
# }


# ecoli 562
# {
#     'apidaesin'    : 0,
#     'midasin'      : 0,
#     'viomycin'     : 0,
#     'paromomycin'  : 0,
#     'blasticidin'  : 0,
#     'kirromycin'   : 0,
#     'puromycin'    : 0,
#     'virginiamycin': 0,
#     'titin'        : 0,
#     'neomycin'     : 0,
#     'spectinomycin': 0,
#     'colicin'      : 0,
#     'hygromycin'   : 0,
#     'erithromycin' : 0,
#     'cathelicidin' : 0
# }
# thermus 300852

# {
#     'paromomycin'  : 0,
#     'puromycin'    : 0,
#     'pactamycin'   : 0,
#     'dityromycin'  : 0,
#     'metalnikowin' : 0,
#     'madumycin'    : 0,
#     'capreomycin'  : 0,
#     'capreomycin'  : 0,
#     'kirromycin'   : 0,
#     'neomycin'     : 0,
#     'erithromycin' : 0,
#     'quinupristin' : 0,
#     'viomycin'     : 0,
#     'negamycin'    : 0,
#     'linopristin'  : 0,
#     'azithromycin' : 0,
#     'telithromycin': 0,
#     'hygromycin'   : 0,
#     'oncocin'      : 0,
#     'colicin'      : 0,
# }
# s.cerevisea 4932


