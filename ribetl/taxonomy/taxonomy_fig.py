from dataclasses import dataclass
import json
from lib2to3.pytree import Node
import os
from platform import architecture
from pprint import pprint
from tokenize import Special
from typing import List
import dotenv
import numpy as np
import pandas as pd
from ete3 import NCBITaxa,  TreeStyle,  faces, AttrFace, TextFace, NodeStyle, CircleFace, ImgFace, RectFace,SeqMotifFace
from sqlalchemy import false
from generate_taxtree import TaxaProfile, profile_taxa

BACTERIA = 2
ARCH     = 2157
EUKARYA  = 2759


global bact_other_counter;
bact_other_counter = 0
global euk_other_counter;
euk_other_counter =0
global arch_other_counter;
arch_other_counter =0

dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get('STATIC_ROOT')

structs       = [*filter(lambda x: os.path.isdir(os.path.join(STATIC_ROOT,x)) ,os.listdir(STATIC_ROOT))]
profiles      = list(map(lambda _: os.path.join(STATIC_ROOT,_,f"{_}.json"),structs))
org_id_arrays = []
ncbi          = NCBITaxa()



MAIN_IDS=[ 'H. sapiens', 'O. cuniculus',  'S. cerevisiae',  'E. coli',  'T. thermophilus',  'S. aureus']


def node_lineage(node):
    # print(node.taxid," is lineage ", ncbi.get_lineage(node.taxid ));
    return ncbi.get_lineage(node.taxid )

def taxid_to_linnaean(spec_taxid:int):
	print("got taxid", spec_taxid)
	_ = str(ncbi.get_taxid_translator([spec_taxid])[spec_taxid]).split(" ")
	if len(_) <2:
		print("[{}] is not a species likely.".format(_))
		return _
	return _[0][0] + ". "+ _[1]




def taxon_remove_if_not_main(t):
	"""get taxon from list, see if it's in the main subset. if not increment the [kingdom]_other counter and pop from list"""
	ncbi    = NCBITaxa()
	lineage = list( ncbi.get_lineage(t) );
	# species_tuple = list(filter(lambda kv: kv[1]=='species' ,ncbi.get_rank(lineage).items()))[0][0]

	if taxid_to_linnaean(t) in MAIN_IDS:
		return False

	else:
		if len(lineage) < 3:
			return False
		kingdom = lineage[2]

		if kingdom == BACTERIA:
			global bact_other_counter
			bact_other_counter = 1
			return True

		if kingdom == EUKARYA:
			global euk_other_counter
			euk_other_counter+=1
			return True
			
		if kingdom == ARCH:
			global arch_other_counter
			arch_other_counter +=1
			return True
		print("error, should reach here")
		exit(2)





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
		
	
taxprofiles:List[TaxaProfile] = list(map(profile_taxa, profiles))
taxons        = list(map( lift_rank, map(lambda _: _.classified_as, taxprofiles)))


print("TAXONS:",taxons)
     
print("Got taxons: ", len(taxons))
for (i,t) in enumerate(taxons):
	print(i,t)
	if taxon_remove_if_not_main(t):

		print("Removing taxon {}: {}".format(t, taxid_to_linnaean(t)))
			
		del taxons[i]
	# print(t)

print("Got taxons after clean ", len(taxons))
print(arch_other_counter)
print(euk_other_counter)
print(bact_other_counter)

ncbi          = NCBITaxa()
ncbi_topology = ncbi.get_topology(taxons)
# print("gto topoloyg", ncbi_topology)




	
print("----------------------++")





def layout(n):

	if 'manual' in n.name:

		taxnamenode = TextFace( 
			f"Species 1\nSpecies2\n..." ,
			fsize    = 16,
			penwidth = 2,
			fstyle   = 'italic',
			fgcolor  = 'black' )

		taxnamenode.margin_bottom = 80
		taxnamenode.rotation      = -50
		faces.add_face_to_node(taxnamenode, n,column=2 , aligned=True)

		textface_node               = TextFace(f"24" , fsize = 16, fgcolor= 'black' )
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

	if 'taxid' not in n.__dict__:
		continue
	# find_main_nodes(n)
	nstyle = NodeStyle()

	# if 'taxid' in n.__dict__:
	# 	print(taxid_to_linnaean(n.taxid))
	
	# if 'manual' in n.name:
	# 	continue
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


	# print(n.taxid)
	# print(ncbi.get_lineage(n.taxid))
	n.set_style(nstyle)



# print(ncbi_topology.search_nodes(taxid=ARCH))

# ncbi_topology.search_nodes(taxid=BACTERIA)[0].add_child(name='OtherB')
# ncbi_topology.search_nodes(taxid=EUKARYA)[0].add_child(name='OtherE')
# ncbi_topology.search_nodes(taxid=ARCHAEA)[0].add_child(name='OtherA')



ncbi_topology.show(tree_style=ts)
# ncbi_topology.render("tree_flipped.svg", tree_style=ts)