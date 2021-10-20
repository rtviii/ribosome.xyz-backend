from asyncio import gather
import json
import os
import pprint
import dotenv
import numpy as np
from ete3 import NCBITaxa
ncbi = NCBITaxa()
# unique_taxa = list(set(unique_taxa))

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


# Void --> List[Path]
def gather_taxa():
    """"""
    structs       = [*filter(lambda x: os.path.isdir(os.path.join(STATIC_ROOT,x)) ,os.listdir(STATIC_ROOT))]
    profiles      = list(map(lambda _: os.path.join(STATIC_ROOT,_,f"{_}.json"),structs))

    org_id_arrays = []

    for profile in profiles:
        org_id_arrays.append(profile_taxa(profile))
    return org_id_arrays

# Path --> Tuple[src_organisms, host_organisms]
def profile_taxa(path:str):
    """"Provided a profile path, extract source and host organisms in a given profile."""

    with open(path,'r') as infile:
        profile   = json.load(infile)
        src_orgs  = []
        host_orgs = []

        for prot in profile['proteins']:
            src_orgs  = [*src_orgs,*prot['src_organism_ids']]
            host_orgs = [*host_orgs,*prot['host_organism_ids']]

        for rna in profile['rnas']:
            src_orgs  = [*src_orgs,*rna['src_organism_ids']]
            host_orgs = [*host_orgs,*rna['host_organism_ids']]

        return [src_orgs, host_orgs]

        
# Void --> List[Int]
def unique_taxa_static():
    src_tax_ids = np.array(gather_taxa())[:,0]
    flat        = []
    for i in src_tax_ids:
        flat = [*flat, *i]
    return list(set(flat))


def generate_tax_trees():
    taxid2name = ncbi.get_taxid_translator(unique_taxa_static() )
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



pprint.pprint(gather_taxa())
# pprint.pprint(unique_taxa_static())
# pprint.pprint(unique_taxa_static())