import json
from ete3 import NCBITaxa
from pprint import pprint
unique_taxa =[ 5661, 9606, 562, 2238, 1299, 68212, 300852, 1961, 1885, 146537, 4932, 36329, 9823, 83333, 1280, 28985, 9986, 224308, 1351, 559292, 312017, 5911, 262724, 32630, 243230, 274, 67351, 1110693, 316385, 272569, 1883, 1931, 585, 93061, 367830, 544404, 7536, 9913, 55431, 37000, 3562, 1177187, 585035, 331111, 83334, 93062, 285006, 10665, 663, 7460, 246196, 759272, 5693, 1773, 574, 1247190, 5811, 5722, 287, 1772, 9739, 284590, 194966, 679895, 5702, 469008, 9615, 6039, 209285, 311400, 10090, 272844, 273057, 69014, 1293037, 2287, 1223565, 1144670, 480119, 1977881, 470, 1960940, 52133, 1217710, 1217649, 1310637, 421052, 1144663, 1310678, 474186, 3702, 575584, 3039, 235221, 5664, 5689, 5691, 353153, 511145, 986, 2697049, 1283, 169963]
ncbi = NCBITaxa()
unique_taxa = list(set(unique_taxa))

# ? These tax ids are a unique scan over all structures in the current instance of neo4j.
"""This script is used to generate taxonomy files that contain the structure's species present in the database:

viruses   .json
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





taxid2name = ncbi.get_taxid_translator(unique_taxa  )
b          = ncbi.get_name_translator (['Bacteria'] )['Bacteria'][0]
a          = ncbi.get_name_translator (['Archaea']  )['Archaea'][0]
e          = ncbi.get_name_translator (['Eukaryota'])['Eukaryota'][0]
v          = ncbi.get_name_translator (['Viruses']  )['Viruses'][0]


e_arr =[]
b_arr =[]
a_arr =[]
v_arr =[]
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
    if v in lin:
        assert(sortedd==False)
        v_arr.append(i)
        sortedd=True
# pprint(taxid2name)


e_dict=[]
b_dict=[]
a_dict=[]
v_dict=[]


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
for g in v_arr:
    id  = [* ncbi.get_taxid_translator( [str(g)] ).keys  () ][0]
    tax = [* ncbi.get_taxid_translator( [str(g)] ).values() ][0]
    v_dict.append({
        "label":tax,
        "value":[ id ],
        "checked":False
    })

v_dict= {
    "label":"Viruses",
    "value":v_arr,
    "checked":False,
    "children":v_dict
}
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


with open('viruses_browser.json','w') as infile:
    json.dump(v_dict,infile)
with open('bacteria_browser.json','w') as infile:
    json.dump(b_dict,infile)
with open('eukaryota_browser.json','w') as infile:
    json.dump(e_dict,infile)
with open('archaea_browser.json','w') as infile:
    json.dump(a_dict,infile)

# pprint(ncbi.get_name_translator(['Eukaryota']))
# pprint(ncbi.get_name_translator(['Archaea']))
# pprint(ncbi.get_name_translator(['Viruses']))