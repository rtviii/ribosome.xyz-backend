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
import re

dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get('STATIC_ROOT')

structs       = [*filter(lambda x: os.path.isdir(os.path.join(STATIC_ROOT,x)) ,os.listdir(STATIC_ROOT))]
profiles      = list(map(lambda _: os.path.join(STATIC_ROOT,_,f"{_}.json"),structs))

org_id_arrays = []


ligands    = pd.read_csv('ribetl/ciftools/ligands.csv')
ligandlike = pd.read_csv('ribetl/ciftools/ligandlike.csv')


factor_descs  = ligandlike['description'].to_list()
factor_taxids = ligandlike['ids'].to_list()

ligand_descs  = ligands['name'].to_list()
ligand_taxids = ligandlike['ids'].to_list()

antibios_ids = []
factors_ids  = []

antibio_reg = r"(\w*(?<!(cha|pro|dom|str|pla))in\b|(\b\w*zyme\b))"

for i, desc in enumerate(factor_descs):
    factor_matches = re.search('factor', desc, re.IGNORECASE)
    factors_ids = [*factors_ids, *factor_taxids[i].split(',')]



for i, desc in enumerate(ligand_descs):
    anti_matches = re.match(antibio_reg, desc, re.IGNORECASE)
    antibios_ids = [*antibios_ids, *ligand_taxids[i].split(',')]


# Paromomycin




print(antibios_ids)
print(factors_ids)

ncbi = NCBITaxa()
t    = ncbi.get_topology(factors_ids)
def node_lineage(node):
    return  ncbi.get_lineage( node.taxid )


BACTERIA        = 2
ARCHAEA         = 2157
EUKARYA         = 2759

def layout(n):
    taxid = n.taxid
    count = 0
    for x in factors_ids:
        if int(x) == int(taxid):
            count+=1

    if count !=0:

        if count > 100:
            fc = 'red'
        elif 100>  count >50 :
            fc = 'orange'
        elif 50>  count >25:
            fc = 'cyan'
        elif 25>  count:
            fc = 'black'
        faces.add_face_to_node(AttrFace ( "sci_name" , fsize =15,                  fgcolor= 'black' ), n, column=0)
        faces.add_face_to_node(TextFace ( f"{count}" , fsize =30,                  fgcolor= fc ), n, column=0)

ts                = TreeStyle()
ts.mode           = "r"
ts.layout_fn      = layout
ts.show_leaf_name = False

ts.legend.add_face(TextFace("Initiation, Elongation and Translation Factors", fsize=100), column=1)
# ts.legend.add_face(TextFace("Distribution of Antibiotics by species", fsize=100), column=1)

for n in t.traverse():
    # n.img_style = custom_style
    nstyle = NodeStyle()

    nstyle["fgcolor"]              = "#0f0f0f"
    nstyle["vt_line_color"]        = "black"
    nstyle["hz_line_color"]        = "black"
    nstyle["vt_line_width"]        = 4
    nstyle["hz_line_width"]        = 4
    nstyle["vt_line_type"]         = 0 # 0 solid, 1 dashed, 2 dotted
    nstyle["hz_line_type"]         = 1
    n.set_style(nstyle)
    if BACTERIA in node_lineage(n):

        nstyle["fgcolor"]              = "#0f0f0f"
        nstyle["vt_line_color"]        = "black"
        nstyle["hz_line_color"]        = "black"
        nstyle["vt_line_width"]        = 4
        nstyle["hz_line_width"]        = 4
        nstyle["vt_line_type"]         = 0 # 0 solid, 1 dashed, 2 dotted
        nstyle["hz_line_type"]         = 1
        nstyle["bgcolor"]              = "PaleGreen"
        n.set_style(nstyle)

    if ARCHAEA in node_lineage(n):

        #      nstyle["vt_line_width"] = 4
        #      nstyle["hz_line_width"] = 4
        nstyle["fgcolor"]              = "#0f0f0f"
        nstyle["vt_line_color"]        = "black"
        nstyle["hz_line_color"]        = "black"
        nstyle["vt_line_width"]        = 4
        nstyle["hz_line_width"]        = 4
        nstyle["vt_line_type"]         = 0 # 0 solid, 1 dashed, 2 dotted
        nstyle["hz_line_type"]         = 1
        nstyle["bgcolor"]              = "PowderBlue"
        n.set_style(nstyle)

    if EUKARYA in node_lineage(n):
        nstyle["fgcolor"]              = "#0f0f0f"
        nstyle["vt_line_color"]        = "black"
        nstyle["hz_line_color"]        = "black"
        nstyle["vt_line_width"]        = 4
        nstyle["hz_line_width"]        = 4
        nstyle["vt_line_type"]         = 0 # 0 solid, 1 dashed, 2 dotted
        nstyle["hz_line_type"]         = 1
        nstyle["bgcolor"] = "OldLace"
        n.set_style(nstyle)


t.show(tree_style=ts)




