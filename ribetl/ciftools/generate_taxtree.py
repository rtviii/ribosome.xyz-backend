from asyncio import gather
from calendar import c
from dataclasses import dataclass
from distutils.util import copydir_run_2to3
import json
import os
from pprint import pprint
from re import L
import sys
from tkinter import ARC
from typing import Dict, List, Tuple
from unicodedata import bidirectional
from xml.dom.minicompat import NodeList
import dotenv
import numpy as np
# from ete2 import Tree, TreeStyle, TextFace
from ete3 import NCBITaxa, PhyloTree, TreeStyle, Tree, faces, AttrFace, NodeStyle,TextFace,RectFace, CircleFace


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

# with open('unique_taxa.txt','w') as infile:
#     for i in unique_taxa:
#         infile.write(str(i)+"\n")
#     infile.close()

dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get('STATIC_ROOT')


def __to_names(d:dict):
    nu = {
    }
    for k in d.keys():
        if type(d[k]) == dict:
            nu[ [*ncbi.get_taxid_translator([int(k)]).values() ] [0]]=__to_names(d[k])
        else:
            nu[[* ncbi.get_taxid_translator([int(k)]).values() ][0]] = d[k]
    return nu

def generate_tax_trees(unique_taxa:List[int]):
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

# Void --> List[Path]
def gather_taxa():
    """"""
    structs       = [*filter(lambda x: os.path.isdir(os.path.join(STATIC_ROOT,x)) ,os.listdir(STATIC_ROOT))]
    profiles      = list(map(lambda _: os.path.join(STATIC_ROOT,_,f"{_}.json"),structs))

    org_id_arrays = []

    for profile in profiles:
        org_id_arrays.append(profile_taxa(profile))
    return org_id_arrays

# Path --> Dict[rcsb_id,
# src_orgs : Array[int],
# host_orgs: Array[int]]
def profile_taxa(path:str):
    """"Provided a profile path, extract source and host organisms in a given profile."""

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

        return {  
            "rcsb_id"         : rcsb_id,
            "source_organisms": src_orgs,
            "host_organisms"  : host_orgs }


@dataclass
class TaxaProfile(): 
      rcsb_id      : str
      classified_as: int
      src_orgs     : Dict[int,int]
      host_orgs    : Dict[int,int]


# Dict[...] -->  TaxaProfile
def classify_profile(d:dict)->TaxaProfile:

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

BACTERIA = 2
ARCHAEA  = 2157
EUKARYA  = 2759

# ncbi.get_taxid_translator([BACTERIA,EUKARYA,ARCHAEA])
# tax_profiles = [classify_profile(_).classified_as for _ in gather_taxa()]
# t            = ncbi.get_topology(tax_profiles)

# def node_lineage(node):
#     return  ncbi.get_lineage( node.taxid )

# # def get_count(node):




# lins = [*sorted([*map(lambda x: ncbi.get_lineage(x), tax_profiles)], key=lambda x: len(x))] 


# # [print(l) for l in lins]

# # [561,562]
# flat_lincount = {
# }

# for i in range(100):
#     for l in lins:
#         if i >len(l)-1:
#             continue
#         else:
#             if l[i] not in flat_lincount:
#                 flat_lincount[l[i]] = 1
#             else:
#                 flat_lincount[l[i]] += 1



# pprint(flat_lincount)


# root= t.get_tree_root()
# def layout(n):

#     # rd = int(t.get_distance(n,root))
#     # if rd ==2:

#     faces.add_face_to_node(TextFace (fsize =20 , text=f"{flat_lincount[n.taxid]}", fgcolor= f"blue"), n, column=2)

#     _ = [* (ncbi.get_taxid_translator([n.taxid])).values() ][0].lower()
#     if "homo" in _ or "coli" in _ or 'thermus' in _:
#         faces.add_face_to_node(AttrFace ( "sci_name" , fsize =20,                  fgcolor= 'black' ), n, column=0)

#     # faces.add_face_to_node(TextFace (              fsize =14 , text=f"distance to root : {t.get_distance(n,root)}", fgcolor= f"purple"), n, column=2)



# ts           = TreeStyle()
# tshz_line_width = 4
# ts.mode      = "r"
# ts.layout_fn = layout

# ts.show_leaf_name = False
# ts.legend.add_face(TextFace("Structure Count", fsize=24), column=1)

# for n in t.traverse():
#     # n.img_style = custom_style
#     nstyle = NodeStyle()

#     nstyle["fgcolor"]              = "#0f0f0f"
#     nstyle["vt_line_color"]        = "black"
#     nstyle["hz_line_color"]        = "black"
#     nstyle["vt_line_width"]        = 4
#     nstyle["hz_line_width"]        = 4
#     nstyle["vt_line_type"]         = 0 # 0 solid, 1 dashed, 2 dotted
#     nstyle["hz_line_type"]         = 1
#     n.set_style(nstyle)
#     if BACTERIA in node_lineage(n):

#         nstyle["fgcolor"]              = "#0f0f0f"
#         nstyle["vt_line_color"]        = "black"
#         nstyle["hz_line_color"]        = "black"
#         nstyle["vt_line_width"]        = 4
#         nstyle["hz_line_width"]        = 4
#         nstyle["vt_line_type"]         = 0 # 0 solid, 1 dashed, 2 dotted
#         nstyle["hz_line_type"]         = 1
#         nstyle["bgcolor"]              = "PaleGreen"
#         n.set_style(nstyle)

#     if ARCHAEA in node_lineage(n):

#         # nstyle["vt_line_width"] = 4
#         # nstyle["hz_line_width"] = 4
#         nstyle["fgcolor"]              = "#0f0f0f"
#         nstyle["vt_line_color"]        = "black"
#         nstyle["hz_line_color"]        = "black"
#         nstyle["vt_line_width"]        = 4
#         nstyle["hz_line_width"]        = 4
#         nstyle["vt_line_type"]         = 0 # 0 solid, 1 dashed, 2 dotted
#         nstyle["hz_line_type"]         = 1
#         nstyle["bgcolor"] = "PowderBlue"
#         n.set_style(nstyle)

#     if EUKARYA in node_lineage(n):
#         nstyle["fgcolor"]              = "#0f0f0f"
#         nstyle["vt_line_color"]        = "black"
#         nstyle["hz_line_color"]        = "black"
#         nstyle["vt_line_width"]        = 4
#         nstyle["hz_line_width"]        = 4
#         nstyle["vt_line_type"]         = 0 # 0 solid, 1 dashed, 2 dotted
#         nstyle["hz_line_type"]         = 1
#         nstyle["bgcolor"] = "OldLace"
#         n.set_style(nstyle)


# t.show(tree_style=ts)




# # ========================================================================================
# bydomain={}
# #depth 1: domains
# for taxon in tax_profiles: 

#     lineage            = ncbi.get_lineage(taxon)
#     lineage_name_pairs = ncbi.get_taxid_translator(lineage)
#     kvps               = [*lineage_name_pairs.keys()]

#     if kvps[1] in bydomain:
#         bydomain[kvps[1]] +=1
#     else:
#         bydomain[kvps[1]] =1

# byclass = {

#     2:{

#     },
#     2759:{
    
#     },
#     2157:{

#     },
#     10239:{}
# }
# #depth 1: domains
# for taxon in tax_profiles: 

#     lineage            = ncbi.get_lineage(taxon)
#     lineage_name_pairs = ncbi.get_taxid_translator(lineage)
#     kvps               = [*lineage_name_pairs.keys()]
#     domain             = kvps[1]
#     cls                = kvps[2]

#     if cls in byclass[domain]:
#         byclass[domain][cls] +=1
#     else:
#         byclass[domain][cls] =1
# pprint(__to_names(bydomain))



# pprint(to_names(byclass))
# generate_tax_trees(list(set(tax_profiles)))


# print(len(tax_profiles))