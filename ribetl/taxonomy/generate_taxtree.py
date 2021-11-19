from dataclasses import dataclass
from importlib.resources import path
import json
import os
from pathlib import Path
from tkinter import ARC
from typing import Dict, List, Tuple
import dotenv
import numpy as np
from ete3 import NCBITaxa, PhyloTree, TreeStyle, Tree, faces, AttrFace, NodeStyle, TextFace, RectFace, CircleFace

BACTERIA = 2
ARCHAEA  = 2157
EUKARYA  = 2759

ncbi = NCBITaxa()
# ncbi.update_taxonomy_database()  #<--------------------------- NCBI Database has to be updated from time to time.

# ? These tax ids are a unique scan over all structures in the current instance of neo4j.
"""This script is used to generate taxonomy files that contain the structure's species present in the database:

bacteria  .json
eukaryota .json
archaea   .json

Most filters on the site rely on these files for situating individual tax ids. 
Hence, these files should be generated anew when new structures are added to the database or extended with tax ids of the new structures.

match (r:RibosomeStructure) 
unwind r._organismId as orgs
return  collect(distinct orgs);
"""


dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get('STATIC_ROOT')


def __to_names(d:dict):
    nu = {}
    for k in d.keys():
        if type(d[k]) == dict:
            nu[ [*ncbi.get_taxid_translator([int(k)]).values() ] [0]]=__to_names(d[k])
        else:
            nu[[* ncbi.get_taxid_translator([int(k)]).values() ][0]] = d[k]
    return nu

def generate_tax_browser_templates(unique_taxa:List[int]):
    taxid2name = ncbi.get_taxid_translator(unique_taxa)
    b          = ncbi.get_name_translator (['Bacteria'] )['Bacteria'][0]
    a          = ncbi.get_name_translator (['Archaea']  )['Archaea'][0]
    e          = ncbi.get_name_translator (['Eukaryota'])['Eukaryota'][0]
    v          = ncbi.get_name_translator (['Viruses']  )['Viruses'][0]

    e_arr =[]
    b_arr =[]
    a_arr =[]

    for i in unique_taxa:
        lin = ncbi.get_lineage(i)
        sortedd = False
        if b in lin:
            b_arr.append(i)
            sortedd=True
        if a in lin:
            assert(sortedd==False)
            a_arr.append(i)
            sortedd=True
        if e in lin:
            assert(sortedd==False)
            e_arr.append(i)
            sortedd=True


    e_dict=[]
    b_dict=[]
    a_dict=[]


    for t in e_arr:
        id  = [* ncbi.get_taxid_translator( [str(t)] ).keys() ][0]
        tax = [*ncbi.get_taxid_translator( [str(t)] ).values()][0]
        e_dict.append({
            "label":tax,
            "value":[ id ],
            "checked":False
        })

    for f in a_arr:
        id  = [*ncbi.get_taxid_translator( [str(f)] ).keys() ][0]
        tax = [*ncbi.get_taxid_translator( [str(f)] ).values()][0]
        a_dict.append({
        "label":tax,
        "value":[ id ],
        "checked":False
        })

    for c in b_arr:
        id  = [* ncbi.get_taxid_translator( [str(c)] ).keys() ][0]
        tax = [*ncbi.get_taxid_translator( [str(c)] ).values()][0]
        b_dict.append({
            "label":tax,
            "value":[ id ],
            "checked":False
        })

    b_dict= {
        "label":"Bacteria",
        "value":b_arr,
        "checked":False,
        "children":b_dict
    }
    a_dict= {
        "label":"Archaea",
        "value":a_arr,
        "checked":False,
        "children":a_dict
    }
    e_dict= {
        "label":"Eukaryota",
        "value":e_arr,
        "checked":False,
        "children":e_dict
    }


    with open('bacteria_browser.json','w') as infile:
        json.dump(b_dict,infile)
    with open('eukaryota_browser.json','w') as infile:
        json.dump(e_dict,infile)
    with open('archaea_browser.json','w') as infile:
        json.dump(a_dict,infile)

# Want to know each structure's tax ids.
# Want to know each tax id's structures.

structs       = [*filter(lambda x: os.path.isdir(os.path.join(STATIC_ROOT,x)) ,os.listdir(STATIC_ROOT))]
profiles      = list(map(lambda _: os.path.join(STATIC_ROOT,_,f"{_}.json"),structs))

# def gather_taxa(profiles:List[Path])->List[dict]:
#     """Given a list of profile paths, returns an array of organism ids."""

#     org_id_arrays = []

#     for profile in profiles:
#         org_id_arrays.append(profile_taxa(str(profile)))

#     return org_id_arrays


@dataclass
class TaxaProfile(): 
      rcsb_id      : str
      classified_as: int
      src_orgs     : Dict[int,int]
      host_orgs    : Dict[int,int]


# Path --> Dict[rcsb_id,
# src_orgs : Array[int],
# host_orgs: Array[int]]
def profile_taxa(path:str)->TaxaProfile:
    """"Provided a profile path, extract source and host organisms in a given profile."""

    # Dict[...] -->  TaxaProfile
    def classify_profile(d:dict)->TaxaProfile:
        """Given a an organism id dictionary of the form(above), classifies a given structure.
        {  
                "rcsb_id"         : str,
                "source_organisms": List[str],
                "host_organisms"  : List[str] 
                }
        """
        s={}
        for _ in d['source_organisms']:
            if _ in s:
                s[_]+=1
            else:
                s[_]= 1

        h={}
        for _ in d['host_organisms']:
            if _ in h:
                h[_]+=1
            else:
                h[_]= 1

        top_org = [* sorted(s.items(), key=lambda x: x[1], reverse=True) ][0][0]

        return TaxaProfile(
            rcsb_id       = d["rcsb_id"],
            classified_as = top_org,
            src_orgs      = s,
            host_orgs     = h,
        )

    with open(path,'r') as infile:
        profile   = json.load(infile)
        rcsb_id   = profile['rcsb_id']
        src_orgs  = []
        host_orgs = []

        for prot in profile['proteins']:
            src_orgs  = [*src_orgs,*prot['src_organism_ids']]
            host_orgs = [*host_orgs,*prot['host_organism_ids']]

        for rna in profile['rnas']:
            src_orgs  = [*src_orgs,*rna['src_organism_ids']]
            host_orgs = [*host_orgs,*rna['host_organism_ids']]

        return classify_profile({  
            "rcsb_id"         : rcsb_id,
            "source_organisms": src_orgs,
            "host_organisms"  : host_orgs })



#! You probably want something like   map (classify_profile(map(profile_taxa, profiles)))