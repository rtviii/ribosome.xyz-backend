from dataclasses import dataclass
import json
from lib2to3.pytree import Node
import os
from pprint import pprint
from typing import List
import dotenv
import numpy as np
import pandas as pd
from ete3 import NCBITaxa,  TreeStyle,  faces, AttrFace, TextFace, NodeStyle, CircleFace, ImgFace, RectFace,SeqMotifFace
from generate_taxtree import TaxaProfile, profile_taxa

BACTERIA = 2
ARCH  = 2157
EUKARYA  = 2759

# BACTERIA_other        = 9992
# ARCHAEA_other        = 9992157
# EUKARYA_other        = 9992759

dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get('STATIC_ROOT')

structs       = [*filter(lambda x: os.path.isdir(os.path.join(STATIC_ROOT,x)) ,os.listdir(STATIC_ROOT))]
profiles      = list(map(lambda _: os.path.join(STATIC_ROOT,_,f"{_}.json"),structs))
org_id_arrays = []
ncbi          = NCBITaxa()



MAIN_IDS=['S. oleracea', 'T. brucei', 'S. scrofa', 'H. sapiens', 'O. cuniculus', 'K. lactis', 'S. cerevisiae', 'C. thermophilum', 'E. coli', 
        'A. baumanii', 'P. aeruginosa', 'B. subtilis', 'F. faecalis', 'T. thermophilus', 'D. radiodurans','M. smegmatis', 'S. aureus']


def node_lineage(node):
    # print(node.taxid," is lineage ", ncbi.get_lineage(node.taxid ));
    return ncbi.get_lineage(node.taxid )

def taxid_to_linnaean(spec_taxid:int):_ = str(ncbi.get_taxid_translator([spec_taxid])[spec_taxid]).split(" "); return _[0][0] + ". "+ _[1]

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
		
def tail_species_to_other(taxid:int ):
	...
	
taxprofiles:List[TaxaProfile] = list(map(profile_taxa, profiles))

taxons      = list(map( lift_rank, map(lambda _: _.classified_as, taxprofiles)))

ncbi          = NCBITaxa()
ncbi_topology = ncbi.get_topology(taxons)
# print("gto topoloyg", ncbi_topology)



print("----------------------++")

# def node_get_upto_kingdom(node):
	# ncbi    = NCBITaxa()
	# nodeid  = node.taxid;
	# lineage = list( ncbi.get_lineage(nodeid) );
	# print("Got lineage", lineage)
	# if len(lineage) < 3:
	# 	return node
	# kingdom = lineage[2]

	# if kingdom == BACTERIA:
	# 	print("adding bacteria other", BACTERIA_other)
	# 	lineage_extended = [*lineage, BACTERIA_other]
	# 	node.taxid = BACTERIA_other
	# if kingdom == EUKARYA:
	# 	print("adding bacteria other", EUKARYA_other)
	# 	lineage_extended = [*lineage, EUKARYA_other]
	# 	node.taxid = EUKARYA_other
	# if kingdom == ARCHAEA:
	# 	print("adding bacteria other", ARCHAEA_other)
	# 	lineage_extended = [*lineage, ARCHAEA_other]
	# 	node.taxid = ARCHAEA_other


		
	# return node 
     

	


 

def find_main_nodes(node):
	nodeid = node.taxid
	lineage = ncbi.get_lineage(nodeid)
	names   = ncbi.translate_to_names(lineage)

	for ( tid,n ) in zip(lineage,names):
		if list( ncbi.get_rank([ tid ]).values() )[0] == 'species':
			linaen = taxid_to_linnaean(tid)
			print(linaen, tid)
			if linaen in MAIN_IDS:
				print("Detected main node ", node.taxid)
				...
			else:
				...

       
		

# for i in [1, 131567, 2759, 33154, 33208, 6072, 33213, 33511, 7711, 89593, 7742, 7776, 117570, 117571, 8287, 1338369, 32523, 32524, 40674, 32525, 9347, 1437010, 314146, 314147, 9989, 1963758, 337687, 10066, 39107, 10088, 862507, 10090]:
#     print(ncbi.get_rank([ i ]))
#     tid_ = list(ncbi.get_rank([ i ]).keys())
#     print(ncbi.get_common_names(tid_))
#     print()
print("----------------------++")



def layout(n):

	if 'manual' in n.name:

		taxnamenode = TextFace( 
			f"Manual X" ,
			fsize    = 16,
			penwidth = 2,
			fstyle   = 'italic',
			fgcolor  = 'black' )

		taxnamenode.margin_bottom = 80
		taxnamenode.rotation      = -50
		faces.add_face_to_node(taxnamenode, n,column=2 , aligned=True)

		textface_node               = TextFace(f"COUHNt" , fsize = 16, fgcolor= 'black' )
		textface_node.rotation      = -90
		textface_node.margin_top    = 10
		textface_node.margin_bottom = 14
		faces.add_face_to_node(textface_node, n,    column   =1, aligned=True)
		return
	tid   = n.taxid
	count = 0
	for s in taxons:
		if s == tid: count+=1

	if tid in taxons:
		taxnamenode = TextFace( 
			f"{taxid_to_linnaean(tid)}" ,
			fsize    = 16,
			penwidth = 2,
			fstyle   = 'italic',
			fgcolor  = 'black' )

		taxnamenode.margin_bottom = 80
		taxnamenode.rotation      = -50
		faces.add_face_to_node(taxnamenode, n,column=2 , aligned=True)

		textface_node               = TextFace(f"{count}" , fsize = 16, fgcolor= 'black' )
		textface_node.rotation      = -90
		textface_node.margin_top    = 10
		textface_node.margin_bottom = 14
		faces.add_face_to_node(textface_node, n,    column   =1, aligned=True)

ts                     = TreeStyle()
ts.mode                = "r"
ts.layout_fn           = layout
ts.show_branch_length  = False
ts.show_branch_support = False
ts.show_scale          = False
ts.show_leaf_name      = False
ts.draw_guiding_lines  = True
ts.rotation            = -270
# ts.rotation            = -90

ncbi_topology.search_nodes(taxid=BACTERIA)[0].add_child(name='Bac_manual')
ncbi_topology.search_nodes(taxid=EUKARYA)[0].add_child(name='Euk_manual')
ncbi_topology.search_nodes(taxid=131567)[0].add_child(name='Arch_manual')

for n in ncbi_topology.traverse():
	# find_main_nodes(n)
	nstyle = NodeStyle()
	
	if 'manual' in n.name:
		continue
	if BACTERIA in node_lineage(n):
		nstyle["bgcolor"] ="#DCFCDE"
	if EUKARYA in node_lineage(n):
		nstyle["bgcolor"] = "#FAF8DD"
	if ARCH in node_lineage(n):
		nstyle["bgcolor"] = "#DDF2FA"
     
	# n = node_get_upto_kingdom(n)

	nstyle["shape"]   = "square"
	nstyle["size"]    = 0
	nstyle["fgcolor"] = "blue"


	# nstyle["vt_line_color"] = 10
	# nstyle["hz_line_color"] = 10
	nstyle["hz_line_width"]  = 1
	nstyle["vt_line_width"]  = 1

	nstyle["fgcolor"]       = "blue"
	nstyle["fgcolor"]       = "blue"


	print(n.taxid)
	print(ncbi.get_lineage(n.taxid))
	n.set_style(nstyle)



# print(ncbi_topology.search_nodes(taxid=ARCH))

# ncbi_topology.search_nodes(taxid=BACTERIA)[0].add_child(name='OtherB')
# ncbi_topology.search_nodes(taxid=EUKARYA)[0].add_child(name='OtherE')
# ncbi_topology.search_nodes(taxid=ARCHAEA)[0].add_child(name='OtherA')



ncbi_topology.show(tree_style=ts)
# ncbi_topology.render("tree_flipped.svg", tree_style=ts)