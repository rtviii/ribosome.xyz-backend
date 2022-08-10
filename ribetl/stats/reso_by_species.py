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


species_dict = {

-1: {
    'count':0,
    'xray_resos':[],
    'cryoem_resos':[] # assert len(xray_resos) + len(cryoem_resos) == count
}
    
}



for path in profiles:
    with open(path,'r') as infile:
        profile = json.load(infile)
        src_ids = profile['src_organism_ids']
        if len(src_ids)<1:continue
        
        # print(profile['host_organism_ids'])