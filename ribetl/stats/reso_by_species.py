import json
import os
from pprint import pprint
import re
import dotenv
import dotenv
import numpy as np
from ete3 import NCBITaxa


# ----------------------------- SETUP -----------------------------------------
# Paper additions: 
# - 3 years of data -- percent between xrays/em
# - 5 years of data -- number of new species

LOOK_AT_LAST_N_YEARS = 5

BACTERIA = 2
ARCHAEA  = 2157
EUKARYA  = 2759

ncbi = NCBITaxa()

dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get('STATIC_ROOT')
structs       = [*filter(lambda x: os.path.isdir(os.path.join(STATIC_ROOT,x)) ,os.listdir(STATIC_ROOT))]
profiles      = list(map(lambda _: os.path.join(STATIC_ROOT,_,f"{_}.json"),structs))


# taxid2name = ncbi.get_taxid_translator(profile['src_organism_ids'])



# ----------------------------- CONTAINERS -----------------------------------------

species_dict = {}
years_dict   = {}

for path in profiles:
    with open(path,'r') as infile:
        profile = json.load(infile)

        

        src_ids      = profile['src_organism_ids']
        organism     = profile['src_organism_ids'][0]
        reso         = profile['resolution']
        method       = profile['expMethod']
        release_date = profile['citation_year']

        
        if len(src_ids)<1 or release_date == None: # ----------------- If no organisms or release date
            continue

        if release_date not in years_dict:
            years_dict[release_date] =[*src_ids]

        else:
            years_dict[release_date] =[*src_ids, *years_dict[release_date]]

# -------------------------------------------------------------------------------- YEARS DICT
sorted_by_year = sorted(years_dict.items(), key=lambda x: x[0])
last10years    = sorted_by_year[-10:]
all            = sorted_by_year

allset = set();
for x in all:
    allset.update(x[1])

last10yearsset = set();
for i in last10years:
    last10yearsset.update(i[1])

print(len(allset))
print(len(last10yearsset))
print(allset.difference(last10yearsset))
print(ncbi.get_taxid_translator(list(allset.difference(last10yearsset))))



# -------------------------------------------------------------------------------- YEARS DICT
########
#######
######
####
##
#

# species_dict[562] ={
#     'count'               :   species_dict[562]['count'              ] + species_dict[83333]['count'               ],
#     'X-RAY DIFFRACTION'   : [*species_dict[562]['X-RAY DIFFRACTION'  ], *species_dict[83333]['X-RAY DIFFRACTION'  ]],
#     'ELECTRON MICROSCOPY' :[* species_dict[562]['ELECTRON MICROSCOPY'], *species_dict[83333]['ELECTRON MICROSCOPY']], 
# }

# del(species_dict[83333])

# sorted_by_count = sorted(species_dict.items(), key=lambda x: x[1]['count'])
# for item in sorted_by_count[-5:]:
#     assert(len(item[1]['X-RAY DIFFRACTION']) +len(item[1]['ELECTRON MICROSCOPY']) == item[1]['count']) 
#     print("Count:",item[1]['count'])
#     print("Species: ", ncbi.get_taxid_translator([ item[0] ]))
#     print(f"Avg xray({len(item[1]['X-RAY DIFFRACTION'])})",  np.mean(item[1]['X-RAY DIFFRACTION']))
#     print(f"Avg em({len(item[1]['ELECTRON MICROSCOPY'])})",  np.mean(item[1]['ELECTRON MICROSCOPY']))
#     print("\n")


# rest_xrays = []
# rest_ems   = []
# rest_count = 0

# for rest_item in sorted_by_count[:-5]:
#     rest_count +=              rest_item[1]['count'               ]
#     rest_xrays =[*rest_xrays, *rest_item[1]['X-RAY DIFFRACTION'  ]]
#     rest_ems   =[*rest_ems, *rest_item[1]['ELECTRON MICROSCOPY']]


# print("\n")
# print("rest Count:", rest_count)
# print(f"rest Avg xray({len(rest_xrays)})",  np.mean(rest_xrays))
# print(f"rest Avg em({len(rest_ems)})",  np.mean(rest_ems))
