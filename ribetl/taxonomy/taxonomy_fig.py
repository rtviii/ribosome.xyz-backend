from asyncio import gather
from cgi import print_environ
from dataclasses import dataclass
import json
from lib2to3.pytree import Node
import os
from pprint import pprint
from re import I
import sys
from turtle import bgcolor, position
from typing import Dict, List, Text, Tuple
import dotenv
import numpy as np
import requests
import pandas as pd
from ete3 import NCBITaxa,  TreeStyle,  faces, AttrFace, TextFace, NodeStyle, CircleFace, ImgFace, RectFace,SeqMotifFace
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
taxons      = list(map( lift_rank, map(lambda _: _.classified_as, taxprofiles)))

ncbi        = NCBITaxa()
t           = ncbi.get_topology(taxons)

def taxid_to_linnaean(spec_taxid:int):_ = str(ncbi.get_taxid_translator([spec_taxid])[spec_taxid]).split(" "); return _[0][0] + ". "+ _[1]
def layout(n):
	tid = n.taxid
	count = 0
	for s in taxons:
		if s == tid: count+=1

	if tid in taxons:
		taxname = TextFace( 
			f"{taxid_to_linnaean(tid)}" ,
			fsize =10,
			fstyle='italic',
			fgcolor= 'black' )

		taxname.margin_left  = 15
		taxname.margin_right = 15



		def getcolor(c):
			if  0 > c > 10:
				return "BlanchedAlmond"
			if  10 > c > 25:
				return "NavajoWhite"
			if  25 > c > 50:
				return "BurlyWood"
			if  50 > c > 100:
				return "SandyBrown"
			if  100 > c > 150:
				return "Sienna"
			if   c > 150:
				return "Maroon"

		C = CircleFace(radius= 4 +  10 * count/250, color=getcolor(count), style="circle")
		C.opacity = 0.6
		C.margin_top         = 5
		C.margin_right       = 5
		C.margin_left        = 5
		C.margin_bottom      = 5
		
		ligands_svg = '/home/rxz/dev/riboxyzbackend/ribetl/taxonomy/ligand_icon.svg'
		imgf= ImgFace(ligands_svg, 30,30)



		rna_rect                    = TextFace("87", fsize=8, fgcolor='PaleTurquoise')
		factors_rect                = TextFace("25", fsize=8, fgcolor='PaleTurquoise')

		rna_rect.margin_top           = 5
		rna_rect.margin_right         = 5
		rna_rect.margin_left          = 5
		rna_rect.margin_bottom        = 5
		rna_rect.opacity              = 0.9
		rna_rect.border.width         = 1
		rna_rect.background.color     = "SteelBlue"

		factors_rect.margin_top         = 5
		factors_rect.margin_right       = 5
		factors_rect.margin_left        = 5
		factors_rect.margin_bottom      = 5

		factors_rect.opacity            = 0.9
		factors_rect.border.width         = 1
		# factors_rect.inner_border.width = 1
		# factors_rect.inner_border.type  = 1
		# factors_rect.border.width       = 1
		factors_rect.background.color   = "DimGray"
		

		rna_rect.hz_align               = 1
		factors_rect.hz_align               = 3

		
		
		seq = ("-----------------------------------------------AQAK---IKGSKKAIKVFSSA---"
		"APERLQEYGSIFTDA---GLQRRPRHRIQSK-------ALQEKLKDFPVCVSTKPEPEDDAEEGLGGLPSN"
		"ISSVSSLLLFNTTENLYKKYVFLDPLAG----THVMLGAETEEKLFDAPLSISKREQLEQQVPENYFYVPD"
		"LGQVPEIDVPSYLPDLPGIANDLMYIADLGPGIAPSAPGTIPELPTFHTEVAEPLKVGELGSGMGAGPGTP"
		"AHTPSSLDTPHFVFQTYKMGAPPLPPSTAAPVGQGARQDDSSSSASPSVQGAPREVVDPSGGWATLLESIR"
		"QAGGIGKAKLRSMKERKLEKQQQKEQEQVRATSQGGHL--MSDLFNKLVMRRKGISGKGPGAGDGPGGAFA"
		"RVSDSIPPLPPPQQPQAEDED----")


		seqFace = SeqMotifFace(seq, gapcolor="blue")

		seqFace = SeqMotifFace(seq[40:100], seq_format="seq")

		faces.add_face_to_node(taxname                   ,                                 n,    column   =0      ,             )
		faces.add_face_to_node(TextFace     ( f"{count}" , fsize = 14, fgcolor= 'black' ), n,    column   =1      , aligned=True)
		# faces.add_face_to_node(C                         ,                                 n, 2, position ="float", aligned=True)
		# faces.add_face_to_node(rna_rect                  ,                                 n, 3, position ="float", aligned=True)
		# faces.add_face_to_node(factors_rect              ,                                 n, 4, position ="float", aligned=True)
		faces.add_face_to_node(seqFace              ,                                 n, 2, position ="float", aligned=True)




ts                = TreeStyle()
ts.mode           = "r"
ts.layout_fn      = layout
ts.show_leaf_name = False



ts.show_branch_length  = False
ts.show_branch_support = False

for n in t.traverse():
	nstyle = NodeStyle()

	if BACTERIA in node_lineage(n):
		nstyle["bgcolor"]              = "PaleGreen"
	if ARCHAEA in node_lineage(n):
		nstyle["bgcolor"]              = "PowderBlue"
	if EUKARYA in node_lineage(n):
		nstyle["bgcolor"] = "OldLace"

	#      nstyle['shape'] = 'circle'
	nstyle["shape"]        = "square"
	nstyle["size"]         = 0
	nstyle["fgcolor"]      = "blue"


	# nstyle["vt_line_color"] = 10
	# nstyle["hz_line_color"] = 10
	nstyle["hz_line_width"]  = 2
	nstyle["vt_line_width"]  = 2

	nstyle["fgcolor"]       = "blue"
	nstyle["fgcolor"]       = "blue"
	n.set_style(nstyle)

t.show(tree_style=ts)