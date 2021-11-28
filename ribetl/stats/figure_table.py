from cgi import print_form
from dataclasses import dataclass
import json
import os
from pprint import pprint
from typing import List
import re
import dotenv
import numpy as np


BACTERIA = 2
ARCHAEA  = 2157
EUKARYA  = 2759

COMPONENTS = [
    # ... Ribosomal Proteins 
    # ... RNA Categories
    "P-Site tRNA" ,
    "E-Site tRNA" ,
    "A-Site tRNA" ,
    "fMet tRNA"   ,
    "Phe tRNA"    ,
    "tRNA"        ,
    "mRNA"        ,
    "Rescue Factor"       ,
    "Elongation Factor"   ,
    "Initiation Factor"   ,
    "Recycling Factor"    ,
    "Release Factor"      ,
    "Transcription Factor",
    "Antibiotic"
]

dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get('STATIC_ROOT')
structs       = [*filter(lambda x: os.path.isdir(os.path.join(STATIC_ROOT,x)) ,os.listdir(STATIC_ROOT))]
profiles      = list(map(lambda _: os.path.join(STATIC_ROOT,_,f"{_}.json"),structs))



def lig_cat(description: str) -> List[str]:

    if "[(" in description.lower():
        return ["Mixed", 'other']

    if len(re.findall(r"(\w*(?<!(cha|pro|dom|str|pla|tox))in\b|(\b\w*zyme\b))", description.lower())) > 0:
        for t in COMPONENTS:
            if t.lower() in description.lower():
                return ['Antibiotics', t]
        return ['Antibiotics', 'Other Antibiotics']

    if len(re.findall(r"(factor)", description.lower())) > 0:
        for f in COMPONENTS:
            if f.lower() in description.lower():
                return ['Factors', f + " Factor"]
        return ['Factors', 'Other Factors']

    if "mrna" in description.lower() or "messenger" in description.lower():
        return ["mRNA", 'Other mRNA']

    if "trna" in description.lower() or "transfer" in description.lower():

        if 'p-' in description.lower():
            ab_class = 'P-Site tRNA'
        elif 'e-' in description.lower():
            ab_class = 'E-Site tRNA'
        elif 'a-' in description.lower():
            ab_class = 'A-Site tRNA'
        elif 'fmet' in description.lower():
            ab_class = "fMet tRNA"
        elif 'phe' in description.lower():
            ab_class = "Phe tRNA"
        else:
            ab_class = 'Other tRNA'

        return ["tRNA", ab_class]

    else:
        return ["Mixed", 'other']

def prof_elemets(path):
    # elements = {
    #     *COMPONENTS
    # }
    cntr  = { 
        "large" : 0,
        "small" : 0,
        "other" : 0
    }


    with open(path,'r') as infile:
        profile   = json.load(infile)
        for prot in profile['proteins']:
            # print(prot['nomenclature'])
            for name in prot['nomenclature']:
                if 'L' in name:
                   cntr [ 'large' ]+=1
                elif 'S' in name:
                   cntr [ 'small' ] += 1
                else:
                   cntr [ 'other' ] +=1
        return [ profile, cntr["large"], cntr["small"], cntr['other'] ]

i = 0
s = 0

LSU  = 0
SSU  = 0
BOTH = 0

for p in profiles:
    i += 1
    [profile, large, small, other]  = prof_elemets(p)

    if ( large ==0 and small ==0 ) or ( large + small < 25 ): 
        print("skipped ", p)
        print(large, small, other)
        s+=1
        continue 

    [ larger, smaller ] = [ large, small ] if large > small else [ small, large ]

    if large == 0:
        SSU +=1 
        continue

    elif small == 0:
        LSU +=1
        continue

    else:
        if larger == large  and larger/smaller > 8:
            LSU +=1
        elif larger == small  and smaller/larger > 8:
            SSU +=1
        else:
            BOTH +=1


print(BOTH, SSU, LSU)
print(BOTH + SSU + LSU)
print(i)
print("skiupped", s)




    

    
    
    


# if abs(small-large):
