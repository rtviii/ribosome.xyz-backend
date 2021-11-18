from nis import cat
import os
from pprint import pprint
import re
import sys
import json
from typing import List, Tuple, TypedDict
import dotenv
from enum import Enum
import numpy as np
import matplotlib.pyplot as plt


def lig_cat(description: str) -> List[str]:

    if "[(" in description.lower():
        return ["Mixed", 'other']

    if len(re.findall(r"(\w*(?<!(cha|pro|dom|str|pla|tox))in\b|(\b\w*zyme\b))", description.lower())) > 0:

        #! E coli implementation
        if 'apidaesin' in description.lower():
            ab_class = 'apidaesin'
        elif 'midasin' in description.lower():
            ab_class = 'midasin'
        elif 'viomycin' in description.lower():
            ab_class = 'viomycin'
        elif 'paromomycin' in description.lower():
            ab_class = 'paromomycin'
        elif 'blasticidin' in description.lower():
            ab_class = 'blasticidin'
        elif 'kirromycin' in description.lower():
            ab_class = 'kirromycin'
        elif 'puromycin' in description.lower():
            ab_class = 'puromycin'
        elif 'virginiamycin' in description.lower():
            ab_class = 'virginiamycin'
        elif 'titin' in description.lower():
            ab_class = 'titin'
        elif 'neomycin' in description.lower():
            ab_class = 'neomycin'
        elif 'spectinomycin' in description.lower():
            ab_class = 'spectinomycin'
        elif 'colicin' in description.lower():
            ab_class = 'colicin'
        elif 'hygromycin' in description.lower():
            ab_class = 'hygromycin'
        elif 'erithromycin' in description.lower():
            ab_class = 'erithromycin'
        elif 'cathelicidin' in description.lower():
            ab_class = 'cathelicidin'
        else:
            ab_class = 'other'

        return ["Antibiotics", ab_class]

    if len(re.findall(r"(factor)", description.lower())) > 0:

        print("Got factor:", description)
        if 'rescue' in description.lower():
            ab_class = 'Ribosome-Rescue Factor'
        elif 'elongation' in description.lower():
            ab_class = 'Elongation Factor'
        elif 'initiation' in description.lower():
            ab_class = 'Translation Inititation'
        elif 'recycling' in description.lower():
            ab_class = 'Ribosome Recycling Factor'
        elif 'releaase' in description.lower():
            ab_class = 'Peptide Chain Release Factor'
        elif 'transcription' in description.lower():
            ab_class = 'Transcription Factor'
        else:
            ab_class = 'other'

        return ["Factors", ab_class]

    if "mrna" in description.lower() or "messenger" in description.lower():

        return ["mRNA", 'other']

    if "trna" in description.lower() or "transfer" in description.lower():

        if 'p-' in description.lower():
            ab_class = 'P-Site'
        elif 'e-' in description.lower():
            ab_class = 'E-Site'
        elif 'a-' in description.lower():
            ab_class = 'A-Site'
        elif 'fmet' in description.lower():
            ab_class = 'Fmet'
        elif 'phe' in description.lower():
            ab_class = 'Phe'
        else:
            ab_class = 'other'

        return ["tRNA", ab_class]

    else:
        return ["Mixed", 'other']


def get_noms(a: dict):
    _ = []
    for i in a:
        for name in a[i]['nomenclature']:
            if name in _:
                print("DUPLICATE")
            else:
                _ = [*_, [name, len(a[i]['residues'])]]
    return _


sys.path.append('/home/rxz/dev/riboxyzbackend/')
dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')

STATIC_ROOT = os.environ.get("STATIC_ROOT")
bsites = []
profiles = []

SPECIES = int(sys.argv[1])


class LigandCategory(Enum):
    tRNA = 'tRNA'
    mRNA = 'mRNA'
    Factor = 'Factor'
    Antibiotics = 'Antibiotics'
    Mixed = 'Mixed'


class LigandLike(TypedDict):
    description: str
    parent: str
    chain: str
    src_organism_ids: List[int]
    path: str
    category: Tuple[LigandCategory, str]


for _ in os.listdir(STATIC_ROOT):
    if len(_) > 4:
        continue
    else:
        profiles = [*profiles,
                    os.path.join(STATIC_ROOT, _, "{}.json".format(_))]


def profile_ligandlikes_by_species(profile_path: str, sought_spec: int) -> List[LigandLike]:
    profile = {}
    _ = []
    with open(profile_path, 'rb') as infile:

        profile = json.load(infile)
        polymers = []
        nonpolymers = []

        if sought_spec not in profile['src_organism_ids']:
            return []

        if profile['proteins'] != None:
            polymers = [*polymers, *profile['proteins']]

        if profile['rnas'] != None:
            polymers = [*polymers, *profile['rnas']]

        try:
            nonpolymers = [*profile['ligands']]

        except:
            ...

    for poly in polymers:
        if bool(poly['ligand_like']) == True:
            _.append({
                'description': poly['rcsb_pdbx_description'],
                'category': lig_cat(poly['rcsb_pdbx_description']),
                'parent': profile['rcsb_id'],
                'chain': poly['auth_asym_id'],
                'src_organism_ids': poly['src_organism_ids'],
                'path': os.path.join(STATIC_ROOT, profile['rcsb_id'].upper(), "POLYMER_{}.json".format(poly['auth_asym_id']))})

    for np in nonpolymers:
        if not "ion" in np['chemicalName'].lower():
            _.append({
                'description': np['chemicalName'],
                'category': lig_cat(np['chemicalName']),
                'parent': profile['rcsb_id'],
                'path': os.path.join(STATIC_ROOT, profile['rcsb_id'].upper(), "LIGAND_{}.json".format(np['chemicalId']))})
    return _


category_counts = {
    'mRNA': {"all": 0},
    'tRNA': {},
    'Factors': {},
    'Antibiotics': {},
}


for ppath in profiles:
    ligandlike = profile_ligandlikes_by_species(ppath, SPECIES)
    counts = [0, 0]

    if len(ligandlike) > 0:
        for l in ligandlike:
            path = l['path']
            desc = l['description']
            category = l['category']
            if category[0] == 'Antibiotics':
                if category[1] in category_counts['Antibiotics']:
                    category_counts['Antibiotics'][category[1]] += 1
                else:
                    category_counts['Antibiotics'][category[1]] = 1
            if category[0] == 'tRNA':
                if category[1] in category_counts['tRNA']:
                    category_counts['tRNA'][category[1]] += 1
                else:
                    category_counts['tRNA'][category[1]] = 1

            if category[0] == 'mRNA':
                category_counts['mRNA']['all'] += 1

            if category[0] == 'Factors':
                if category[1] in category_counts['Factors']:
                    category_counts['Factors'][category[1]] += 1
                else:
                    category_counts['Factors'][category[1]] = 1


pprint(category_counts)
# labels  = 'Antibiotics', 'mRNA/tRNA/Factor', 'Both'
# sizes   = [p_counts['antibiotic'], p_counts['both'], p_counts['mrna_trna_factor']]
# explode = (0, 0.1, 0)
# # +
# # labels  = 'Antibiotics', 'mRNA', 'tRNA', 'Factors'
# # sizes   = [ category_counts['Antibiotics'], category_counts['mRNA'], category_counts['tRNA'],category_counts['Factors'] ]
# # explode = (0, 0.1, 0,0)  # only "explode" the 2nd slice (i.e. 'Hogs')


# fig1, ax1 = plt.subplots()
# ax1.pie(sizes, explode=explode, labels=labels, autopct='%1.1f%%',
#         shadow=True, startangle=90)
# ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

# plt.show()
# !------------------------------------------------------------------------------#


# Make data: I have 3 groups and 7 subgroups
group_names    = ['Antibiotics', 'mRNA', 'tRNA', 'Factors']
group_size     = []


subgroup_names = []
subgroup_size  = []

for name, count in category_counts['Antibiotics'].items():
    subgroup_names.append(name)
    subgroup_size.append(count)
group_size.append(sum(category_counts['Antibiotics'].values()))

for name, count in category_counts['mRNA'].items():
    subgroup_names.append(name)
    subgroup_size.append(count)
group_size.append(sum(category_counts['mRNA'].values()))

for name, count in category_counts['tRNA'].items():
    subgroup_names.append(name)
    subgroup_size.append(count)
group_size.append(sum(category_counts['tRNA'].values()))
    
for name, count in category_counts['Factors'].items():
    subgroup_names.append(name)
    subgroup_size.append(count)
group_size.append(sum(category_counts['Factors'].values()))
    
print(group_names)

# subgroup_names = ['A.1', 'A.2', 'A.3', 'B.1', 'B.2', 'C.1', 'C.2', 'C.3', 'C.4', 'C.5']
# subgroup_size  = [4,3,5,6,5,10,5,5,4,6]

# Create colors
# subgroup_names_legs=['A.1:a1desc', 'A.2:a2desc', 'A.3:a3desc', 
# 'B.1:b1desc', 'B.2:b2desc', 'C.1:c1desc', 'C.2:c2desc', 'C.3:c3desc', 
# 'C.4:c4desc', 'C.5:c5desc']
# plt.legend(subgroup_names_legs,loc='best')

# First Ring (outside)
fig, ax = plt.subplots()
ax.axis('equal')
a, b, c, d=[plt.cm.Blues, plt.cm.Reds, plt.cm.Greens, plt.cm.Reds]

mypie, _ = ax.pie(
    group_size,
    radius = 0.6,
    labels = group_names,
    labeldistance = 0.6,
    colors = [a(0.9), a(0.6), a(0.3), a(0.1) ]
    )

mypie2, _ = ax.pie(
    subgroup_size,
    radius        = 1.3,
    # labels        = subgroup_names,
    labeldistance = 0.5,
    # colors        = [a(0.5), a(0.4), a(0.3), b(0.5), b(0.4), c(0.6), c(0.5), c(0.4), c(0.3), c(0.2)]
    )

plt.setp( mypie, width=0.3, edgecolor='black')
plt.setp( mypie2, width=0.3, edgecolor='blue')
plt.margins(2,4)

plt.legend(loc=(0.9, 0.1))
handles, labels = ax.get_legend_handles_labels()

# ax.legend(handles[3:], subgroup_names_legs, loc=(0.9, 0.1))
plt.show()
# !------------------------------------------------------------------------------#


# ? taxids of interest
# Curate antibiotics manually for the 4 species; trna/factor categories are a little easier to curate; mrna is one solid category
#
# human 9606, antibiotics:

# {
#     'streptomycin': 0,
#     'bystin'      : 0,
#     'cadherin'    : 0,
#     'listerin'    : 0,
#     'hygromycin'  : 0
# }


# ecoli 562
# {
#     'apidaesin'    : 0,
#     'midasin'      : 0,
#     'viomycin'     : 0,
#     'paromomycin'  : 0,
#     'blasticidin'  : 0,
#     'kirromycin'   : 0,
#     'puromycin'    : 0,
#     'virginiamycin': 0,
#     'titin'        : 0,
#     'neomycin'     : 0,
#     'spectinomycin': 0,
#     'colicin'      : 0,
#     'hygromycin'   : 0,
#     'erithromycin' : 0,
#     'cathelicidin' : 0
# }
# thermus 300852

# {
#     'paromomycin'  : 0,
#     'puromycin'    : 0,
#     'pactamycin'   : 0,
#     'dityromycin'  : 0,
#     'metalnikowin' : 0,
#     'madumycin'    : 0,
#     'capreomycin'  : 0,
#     'capreomycin'  : 0,
#     'kirromycin'   : 0,
#     'neomycin'     : 0,
#     'erithromycin' : 0,
#     'quinupristin' : 0,
#     'viomycin'     : 0,
#     'negamycin'    : 0,
#     'linopristin'  : 0,
#     'azithromycin' : 0,
#     'telithromycin': 0,
#     'hygromycin'   : 0,
#     'oncocin'      : 0,
#     'colicin'      : 0,
# }
# s.cerevisea 4932
