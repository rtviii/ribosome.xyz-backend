import os
import re
import sys
from pprint import pprint
import json
from typing import List
from Bio import Align
from Bio.SeqRecord import SeqRecord
import dotenv
from pymol import cmd
import argparse
import glob
from ribetl.ciftools.bsite_mixed import BindingSite
from Bio.Align import MultipleSeqAlignment
import numpy as np
from Bio import pairwise2
import matplotlib.pyplot as plt





def lig_cat(description:str)->str:

	if len(re.findall(r"(\w*(?<!(cha|pro|dom|str|pla))in\b|(\b\w*zyme\b))", description.lower()))> 0:
		return "Antibiotics"

	if len(re.findall(r"(factor)", description.lower()))> 0:
		return "Factors"

	if "mrna" in description.lower() or "messenger" in description.lower():
		return "mRNA"

	if "trna" in description.lower() or "transfer" in description.lower():
		return "tRNA"

	return "Mixed Ligands"

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

for _ in os.listdir(STATIC_ROOT):
    if len(_) > 4:
        continue
    else:
        profiles = [*profiles,
                    os.path.join(STATIC_ROOT, _, "{}.json".format(_))]

def get_liglike_species(profile_path: str, sought_spec: int):
    profile = {}
    _ = []

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
                'parent'          : profile['rcsb_id'],
                'chain'           : poly['auth_asym_id'],
                'src_organism_ids': poly['src_organism_ids'],
                'path'            : os.path.join(STATIC_ROOT, profile['rcsb_id'].upper(), "POLYMER_{}.json".format(poly['auth_asym_id']))})

    for np in nonpolymers:
        if not "ion" in np['chemicalName'].lower():
            _.append({
                'description': np['chemicalName'],
                'parent'     : profile['rcsb_id'],
                'path'       : os.path.join(STATIC_ROOT, profile['rcsb_id'].upper(), "LIGAND_{}.json".format(np['chemicalId']))})
    return _

def collect_neighborhood(ligpath: str):
    """We want to extract all nomenclature classes that this """
    try:
        with open(ligpath, 'rb') as lig_infile:
            data = json.load(lig_infile)
            # pprint(data)
    except:
        return []

    # for chain in data:
    # 	if nomclass in data[chain]['nomenclature']:
    # 		# return data[chain]['sequence']
    # 		return  [description, data[chain]['residues'] ]

    # return []


# ※ p-charts


for ppath in profiles: 
    polymer_ligand_map= {
        # 'ligand_name' :{
        #     'poly_class': 0,
        #     'poly_class': 0,
        #     'poly_class': 0,
        #     'poly_class': 0,
        #     'poly_class': 0,
        # }
    }

    # a list of ligandlikes in a [human] struct
    ligandlike = get_liglike_species(ppath, 9606)

    if len( ligandlike )  > 0:
        path       = ligandlike[0]['path']
        desc       = ligandlike[0]['description']
        collect_neighborhood(path)

        print("<<<<<<<<<<>>>>>>>>>>>>hh", desc)

        exit()



