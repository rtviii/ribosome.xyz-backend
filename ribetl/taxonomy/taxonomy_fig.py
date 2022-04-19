from dataclasses import dataclass
import json
from lib2to3.pytree import Node
import os
from typing import List
import dotenv
import numpy as np
import pandas as pd
from ete3 import NCBITaxa,  TreeStyle,  faces, AttrFace, TextFace, NodeStyle, CircleFace, ImgFace, RectFace, SeqMotifFace
from generate_taxtree import TaxaProfile, profile_taxa

BACTERIA = 2
ARCHAEA = 2157
EUKARYA = 2759


ncbi = NCBITaxa()

INCLUDE_ONLY = [209285, 3562, 5691, 9823, 9606,
                9986, 209285, 4932, 562, 470, 1772, 274, 274]


def exclude_nonmain(taxid):
    Lineage = ncbi.get_lineage(taxid)
    print(Lineage)
    for i in INCLUDE_ONLY:
        if i in Lineage:
            return True
    return False


dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get('STATIC_ROOT')

structs = [
    *filter(lambda x: os.path.isdir(os.path.join(STATIC_ROOT, x)), os.listdir(STATIC_ROOT))]
profiles = list(map(lambda _: os.path.join(
    STATIC_ROOT, _, f"{_}.json"), structs))
org_id_arrays = []


def node_lineage(node):
    return ncbi.get_lineage(node.taxid)


def lift_rank(taxid: int) -> int:
    """Given a taxid, make sure that it's a SPECIES (as opposed to strain, subspecies, isolate, norank etc.)"""
    if ncbi.get_rank([taxid])[taxid] == 'species':

        return taxid
    else:
        lin = iter(ncbi.get_lineage(taxid))
        node = 1
        while ncbi.get_rank([node])[node] != 'species':
            node = next(lin)
        return node


taxprofiles: List[TaxaProfile] = list(map(profile_taxa, profiles))
taxons = list(map(lift_rank, map(lambda _: _.classified_as, taxprofiles)))
ncbi = NCBITaxa()
print("len taxons", len(taxons))
taxons = list(filter(exclude_nonmain, taxons))
# print("len mainonly:", len(main_taxa_only))



t = ncbi.get_topology(taxons)

def taxid_to_linnaean(spec_taxid: int): _ = str(ncbi.get_taxid_translator(
    [spec_taxid])[spec_taxid]).split(" "); return _[0][0] + ". " + _[1]


def layout(n):
    if n.name == 'Archaea': n.taxid = 2157

    tid   = n.taxid
    count = 0
    for s in taxons:
        if s == tid:
            count += 1

    if tid in taxons:
        taxnamenode   = TextFace(f"{taxid_to_linnaean(tid)}",fsize=17,penwidth=2,fstyle='italic',fgcolor='black')
        textface_node = TextFace(f"{count}", fsize=16, fgcolor='black')

        # taxnamenode.margin_top    = 20
        taxnamenode.margin_bottom = 60
        taxnamenode.margin_left   = 10
        # taxnamenode.margin_right  = 40
        taxnamenode.rotation      = -30

        textface_node.rotation      = 270
        textface_node.margin_bottom = 60

        faces.add_face_to_node(taxnamenode, n, column=2, aligned=True)
        faces.add_face_to_node(textface_node, n,    column=1, aligned=True)

    elif tid == 2157:

        taxnamenode   = TextFace(f"Archaea\n(4 species)",fsize    = 17,penwidth = 2,fstyle   = 'italic',fgcolor  = 'black')
        textface_node = TextFace("8", fsize=16, fgcolor='black')

        # taxnamenode.margin_top    = 20
        taxnamenode.margin_bottom = 60
        taxnamenode.margin_left   = 10
        # taxnamenode.margin_right  = 40
        taxnamenode.rotation      = -30

        textface_node.rotation      = 270
        textface_node.margin_bottom = 60

        faces.add_face_to_node(taxnamenode, n, column=2, aligned=True)
        faces.add_face_to_node(textface_node, n, column=1, aligned=True)

ts                     = TreeStyle()
ts.mode                = "r"
ts.layout_fn           = layout
ts.show_branch_length  = False
ts.show_branch_support = False
ts.show_scale          = False
ts.show_leaf_name      = False
ts.draw_guiding_lines  = True
ts.rotation            = -270


# -/-/-/-/--/-/-/-/-/- Archaeal node
treeroot = t.get_tree_root()
archnode = treeroot.add_child(name="Archaea")
nstyle   = NodeStyle()

nstyle["bgcolor"]       = "#DDF2FA"
nstyle["shape"]         = "square"
nstyle["size"]          = 1
nstyle["fgcolor"]       = "blue"
nstyle["hz_line_width"] = 2
nstyle["vt_line_width"] = 2
nstyle["fgcolor"]       = "blue"
nstyle["fgcolor"]       = "blue"
archnode.set_style(nstyle)
# -/-/-/-/--/-/-/-/-/- Archaeal node


for n in t.traverse():
    nstyle = NodeStyle()
    if n.name == 'Archaea':
        continue
    

    if BACTERIA in node_lineage(n):
        nstyle["bgcolor"] = "#DCFCDE"
    if ARCHAEA in node_lineage(n):
        nstyle["bgcolor"] = "#DDF2FA"
    if EUKARYA in node_lineage(n):
        nstyle["bgcolor"] = "#FAF8DD"

    nstyle["shape"]         = "square"
    nstyle["size"]          = 1
    nstyle["fgcolor"]       = "blue"
    nstyle["hz_line_width"] = 2
    nstyle["vt_line_width"] = 2
    nstyle["fgcolor"]       = "blue"
    nstyle["fgcolor"]       = "blue"
    n.set_style(nstyle)





t.show(tree_style=ts)
# t.render("tree_flipped.svg", tree_style=ts)
