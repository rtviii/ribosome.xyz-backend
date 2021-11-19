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
from ete3 import NCBITaxa,  TreeStyle,  faces, AttrFace, TextFace, NodeStyle, CircleFace, ImgFace, RectFace,SeqMotifFace
from ..taxonomy.generate_taxtree import TaxaProfile, profile_taxa


#! E coli implementation
ECOLI_TARGETS = [
    "Apidaesin"    ,
    "Midasin"      ,
    "Viomycin"     ,
    "Paromomycin"  ,
    "Blasticidin"  ,
    "Kirromycin"   ,
    "Puromycin"    ,
    "Virginiamycin",
    "Titin"        ,
    "Neomycin"     ,
    "Spectinomycin",
    "Colicin"      ,
    "Hygromycin"   ,
    "Erithromycin" ,
    "Cathelicidin" ,
]

HUMAN_TARGETS = [
]

YEAST_TARGETS = [
]

THERMUS_TARGETS = [
]


FACTORS = [
    'Rescue',
    'Elongation',
    'Initiation',

    'Recycling',
    'Release',
    'Transcription'
]
def lig_cat(description: str) -> List[str]:

    if "[(" in description.lower():
        return ["Mixed", 'other']

    if len(re.findall(r"(\w*(?<!(cha|pro|dom|str|pla|tox))in\b|(\b\w*zyme\b))", description.lower())) > 0:

        for t in ECOLI_TARGETS:
            if t.lower() in description.lower():
                return ['Antibiotics', t]
        return ['Antibiotics', 'other']


    if len(re.findall(r"(factor)", description.lower())) > 0:

        for f in FACTORS:
            if f.lower() in description.lower():
                return ['Factors', f + " Factor"]
        return ['Factors', 'other']

    if "mrna" in description.lower() or "messenger" in description.lower():

        return ["mRNA", 'other']

    if "trna" in description.lower() or "transfer" in description.lower():

        if 'p-' in description.lower():
            ab_class = 'P-Site tRNA'
        elif 'e-' in description.lower():
            ab_class = 'E-Site tRNA'
        elif 'a-' in description.lower():
            ab_class = 'A-Site tRNA'
        elif 'fmet' in description.lower():
            ab_class = 'Fmet tRNA'
        elif 'phe' in description.lower():
            ab_class = 'Phe tRNA'
        else:
            ab_class = 'Other tRNA'

        return ["tRNA", ab_class]

    else:
        return ["Mixed", 'other']


def get_noms(a: dict):
    _ = []
    for i in a:
        for name in a[i]['nomenclature']:
            if name in _:
                # print("DUPLICATE")
                ...
            else:
                _ = [*_, [name, len(a[i]['residues'])]]
    return _


sys.path.append('/home/rxz/dev/riboxyzbackend/')
dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')

STATIC_ROOT = os.environ.get("STATIC_ROOT")
bsites      = []
profiles    = []

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
    """Return the a list of the ligand-like elements associated with a given ribosome structure profile.
    Sought species is looked up in the tax tree to correct for structures which contain substrain ids.
    """

    tax_profile = profile_taxa(profile_path)
    ncbi        = NCBITaxa()
    # unroll full lineage of the given structure's source 
    lng         = ncbi.get_lineage( tax_profile.classified_as )

    profile = {}
    _       = []

    with open(profile_path, 'rb') as infile:

        
        profile     = json.load(infile)
        polymers    = []
        nonpolymers = []

        # verify that the sought species taxid is present in the lineage
        if sought_spec not in lng:return []

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
                'description'     : poly['rcsb_pdbx_description'],
                'category'        : lig_cat(poly['rcsb_pdbx_description']),
                'parent'          : profile['rcsb_id'],
                'chain'           : poly['auth_asym_id'],
                'src_organism_ids': poly['src_organism_ids'],
                'path'            : os.path.join(STATIC_ROOT, profile['rcsb_id'].upper(), "POLYMER_{}.json".format(poly['auth_asym_id']))})

    for np in nonpolymers:
        if not "ion" in np['chemicalName'].lower():
            _.append({
                'description': np['chemicalName'],
                'category'   : lig_cat(np['chemicalName']),
                'parent'     : profile['rcsb_id'],
                'path'       : os.path.join(STATIC_ROOT, profile['rcsb_id'].upper(), "LIGAND_{}.json".format(np['chemicalId']))})
    return _


category_counts = {
    'mRNA'       : {"all": 0},
    'tRNA'       : {},
    'Factors'    : {},
    'Antibiotics': {},
}


for ppath in profiles:
    ligandlike = profile_ligandlikes_by_species(ppath, SPECIES)
    counts     = [0, 0]

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


# pprint(category_counts)
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
# !----------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# Make data: I have 3 groups and 7 subgroups

group_names    = ['mRNA',  'tRNA', 'Factors', 'Antibiotics']
group_size     = []

subgroup_names  = []
subgroup_size   = []
subgroup_colors = []

a, b, c, d = [ plt.cm.Blues,  plt.cm.Greens,plt.cm.Reds, plt.cm.Oranges ]
ab_cs      = 'blue'
mrna_cs    = 'green'
trna_cs    = 'red'
factors_cs = 'orange'


for name, count in category_counts['mRNA'].items():
    subgroup_names.append(name)
    subgroup_size.append(count)
group_size.append(sum(category_counts['mRNA'].values()))

aval =0.1
for name, count in category_counts['tRNA'].items():
    subgroup_names.append(name)
    subgroup_size.append(count)
    subgroup_colors.append(b(aval))
    aval+=0.1
group_size.append(sum(category_counts['tRNA'].values()))
    
bval = 0.1
for name, count in category_counts['Factors'].items():
    subgroup_names.append(name)
    subgroup_size.append(count)
    subgroup_colors.append(c(bval))
    bval+=0.1
group_size.append(sum(category_counts['Factors'].values()))

cval =0.1
for name, count in category_counts['Antibiotics'].items():
    subgroup_names.append(name)
    subgroup_size.append(count)
    subgroup_colors.append(d(cval))
    cval+=0.1
group_size.append(sum(category_counts['Antibiotics'].values()))

# ※----------------------------------------------------------------------------------------------------------------------------------------------------------------------#
    
fig, ax = plt.subplots()
ax.axis('equal')

mypie, text = ax.pie(
    group_size,
    radius        = 0.7,
    labels        = group_names,
    labeldistance = 0.7,
    colors        = [ab_cs, mrna_cs, trna_cs, factors_cs ])

mypie2, text2 = ax.pie(
    subgroup_size,
    radius        = 0.95,
    labels        = subgroup_size,
    labeldistance = 0.9,
    colors        = subgroup_colors)

bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
kw         = dict(arrowprops=dict(arrowstyle="-"),bbox=bbox_props, zorder=0, va="center")

mypie2[0].set_visible(False)

for i, p in enumerate(mypie2):

    # print(p)
    ang = (p.theta2 - p.theta1)/2. + p.theta1

    y                   = np.sin(np.deg2rad(ang))
    x                   = np.cos(np.deg2rad(ang))

    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle     = "angle,angleA=0,angleB={}".format(ang)

    kw["arrowprops"].update({"connectionstyle": connectionstyle
    })

    # if subgroup_names[i] not in [ *ECOLI_TARGETS, 'all', 'other']: 
    ax.annotate(
        subgroup_names[i],
        xy=(x, y),
        xytext=(1.5*np.sign(x),1.4*y),
        horizontalalignment=horizontalalignment,
            **kw)


plt.setp    ( mypie , width=0.25, edgecolor='black')
plt.setp    ( mypie2, width=0.15, edgecolor='black' )

plt.margins ( 2,4 )

plt.legend(loc=(0.9, 0.1))
handles, labels = ax.get_legend_handles_labels()

plt.show()
# !------------------------------------------------------------------------------#


# ? taxids of interest
# Curate antibiotics manually for the 4 species; trna/factor categories are a little easier to curate; mrna is one solid category
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

print(f"Counted {globcount} in total for ", 300852)