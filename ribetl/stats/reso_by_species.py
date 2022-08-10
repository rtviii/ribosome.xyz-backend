import json
import os
from pprint import pprint
import re
import dotenv
import dotenv
import numpy as np
from ete3 import NCBITaxa

BACTERIA = 2
ARCHAEA  = 2157
EUKARYA  = 2759

ncbi = NCBITaxa()

dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get('STATIC_ROOT')
structs       = [*filter(lambda x: os.path.isdir(os.path.join(STATIC_ROOT,x)) ,os.listdir(STATIC_ROOT))]
profiles      = list(map(lambda _: os.path.join(STATIC_ROOT,_,f"{_}.json"),structs))


# taxid2name = ncbi.get_taxid_translator(profile['src_organism_ids'])


species_dict = {}



for path in profiles:
    with open(path,'r') as infile:
        profile = json.load(infile)
        src_ids = profile['src_organism_ids']
        if len(src_ids)<1:continue
        
        organism = src_ids[0]
        reso     = profile['resolution']
        method   = profile['expMethod']

        if organism not in species_dict:
            species_dict[organism] = {
                'X-RAY DIFFRACTION'  : [],
                'ELECTRON MICROSCOPY': [],
                'count':1
            }
        else:
            species_dict[organism]['count']+=1

        species_dict[organism][method].append(reso)




sorted_by_count = sorted(species_dict.items(), key=lambda x: x[1]['count'])
pprint(sorted_by_count[-5:])