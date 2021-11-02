from asyncio import gather
from cgi import print_environ
from dataclasses import dataclass
import json
import os
from pprint import pprint
import sys
from typing import Dict, List, Tuple
import dotenv
import numpy as np
import requests
import pandas as pd
from ete3 import NCBITaxa,  TreeStyle,  faces, AttrFace, TextFace, NodeStyle
from .generate_taxtree import classify_profile, profile_taxa
BACTERIA        = 2
ARCHAEA         = 2157
EUKARYA         = 2759

dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get('STATIC_ROOT')

structs       = [*filter(lambda x: os.path.isdir(os.path.join(STATIC_ROOT,x)) ,os.listdir(STATIC_ROOT))]
profiles      = list(map(lambda _: os.path.join(STATIC_ROOT,_,f"{_}.json"),structs))

org_id_arrays = []
ncbi = NCBITaxa()

def node_lineage(node):
    return  ncbi.get_lineage( node.taxid )


def layout(n):
    taxid = n.taxid
    count = 0
    for x in taxons:
        if int(x) == int(taxid):
            count+=1

        faces.add_face_to_node(AttrFace ( "sci_name" , fsize =15,                  fgcolor= 'black' ), n, column=0)
        faces.add_face_to_node(AttrFace ( "rank" , fsize =15,                  fgcolor= 'black' ), n, column=0)
        # faces.add_face_to_node(TextFace ( f"{count}" , fsize =30,                  fgcolor= 'red' ), n, column=0)



def lift_rank(taxid:int)->int:
	"""Given a taxid, make sure that it's a SPECIES (as opposed to strain, subspecies, isolate, norank etc.)"""
	if ncbi.get_rank([taxid])[taxid] == 'species':
		return taxid
	else:
		lin = iter(ncbi.get_lineage(taxid))
		node = 1
		while ncbi.get_rank([node])[node] != 'species':
			node = next(lin)
		return node
		
	
# for i in taxons:
# 	print("\n")
# 	print(i)
# 	i = lift_rank(i)
# 	print(ncbi.get_rank([ i ]))
# 	rank = ncbi.get_rank([i])[i]
# 	print(ncbi.get_lineage( i ))
# 	print(ncbi.get_taxid_translator([ i ]))
# 	print(lift_rank(i))

taxprofiles = list(map(classify_profile,map(profile_taxa, profiles)))
taxons = list(map( lift_rank, map(lambda _: _.classified_as, taxprofiles)))

ncbi = NCBITaxa()
t    = ncbi.get_topology(taxons)



# print(ncbi.get_taxid_translator([ lift_rank(9606) ]))


# for  t in taxons:
# 	print(t)
# 	print(ncbi.get_rank([t]))
# 	print(ncbi.get_lineage(t))



# ts                = TreeStyle()
# ts.mode           = "c"
# ts.layout_fn      = layout
# ts.show_leaf_name = False


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

#         #      nstyle["vt_line_width"] = 4
#         #      nstyle["hz_line_width"] = 4
#         nstyle["fgcolor"]              = "#0f0f0f"
#         nstyle["vt_line_color"]        = "black"
#         nstyle["hz_line_color"]        = "black"
#         nstyle["vt_line_width"]        = 4
#         nstyle["hz_line_width"]        = 4
#         nstyle["vt_line_type"]         = 0 # 0 solid, 1 dashed, 2 dotted
#         nstyle["hz_line_type"]         = 1
#         nstyle["bgcolor"]              = "PowderBlue"
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