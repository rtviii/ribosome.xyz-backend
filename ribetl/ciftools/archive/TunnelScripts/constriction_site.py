from itertools import chain
import os,sys
from Bio.PDB.Chain import Chain
from typing import List
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from dotenv import load_dotenv
import numpy as np
from scipy.spatial.kdtree import KDTree

def root_self(rootname:str='')->str:
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    root=os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(root)
    load_dotenv(os.path.join(root,'.env'))

if __name__=="__main__":
    root_self('ribxz')

from ciftools.Structure import fetchStructure
from ciftools.TunnelScripts.TunnelLog import (Log)
from ciftools.Neoget import _neoget

log = Log(os.getenv("TUNNEL_LOG"))

def getConstrictedProteins(pdbid:str)->List[str]:

    cypher="""
    match (n:RibosomeStructure{{rcsb_id:"{pdbid}"}})-[]-(rp:RibosomalProtein) 
    where "uL4" in rp.nomenclature or "uL22" in rp.nomenclature
    return rp.entity_poly_strand_id,rp.nomenclature, n.rcsb_id
    """.format_map({"pdbid":pdbid.upper()})
    response = _neoget(cypher)
    resplist= [*map(lambda x: ( x['rp.nomenclature'][0],x['rp.entity_poly_strand_id'] ), response)]
    uL22 = None
    uL4  = None
    for nom in resplist:

        if nom[0]=="uL22":
            uL22=nom[1]
        if nom[0]=="uL4":
            uL4=nom[1]

    return {
        "uL4" : uL4,
        "uL22": uL22
    }


def calc_constriction_site(pdbid:str):
    """Here the assumption is that the constriction site is unique and so is then the KDTree min then."""
    pdbid= pdbid.upper()

    chainL22:str                    = log.get_struct(pdbid)['uL22'].values[0]
    chainL4  :str                   = log.get_struct(pdbid)['uL4'].values[0]

    if ',' in chainL22:chainL22 = chainL22.split(',')[0]
    if ',' in chainL4:chainL4   = chainL4.split(',')[0]

    struct = fetchStructure(pdbid)[0]

    L4 :Chain = struct[chainL4]
    L22:Chain = struct[chainL22]


    l4res       = [ *L4.get_atoms() ]
    l22res      = [ *L22.get_atoms() ]

    kdtreeonl4  = KDTree([*  map(lambda x : x.get_coord(), l4res)  ])
    kdtreeonl22 = KDTree([*  map(lambda x : x.get_coord(), l22res)  ])

    nbridin22: Atom = None
    nbridin4 : Atom = None
    atom     : Atom

    distances_indexes =kdtreeonl4.query([* map(lambda x: x.get_coord(),l22res) ])
    dist = 99999999
    for x in zip( [* distances_indexes[0] ],[* distances_indexes[1] ] ):
        if x[0] < dist:
            dist     = x[0]
            nbridin4 = x[1]

    distances_indexes =kdtreeonl22.query([* map(lambda x: x.get_coord(),l4res) ])
    dist = 99999999
    for x in zip( [* distances_indexes[0] ],[* distances_indexes[1] ] ):
        if x[0] < dist:
            dist     = x[0]
            nbridin22 = x[1]


    l4atomcord  = l4res[nbridin4].get_coord()
    l22atomcord = l22res[nbridin22].get_coord()
    
    centerline=np.mean([ l4atomcord,l22atomcord ], axis=0)
    residueInL4:Residue  = l4res[nbridin4].get_parent()
    residueInL22:Residue = l22res[nbridin22].get_parent()
    

    print("""{} {} in chain {}({}) is closest to {} {} in chain{}({}).\n Centerline:{}"""
    .format(
     residueInL22.get_resname(),
     residueInL22.get_id()[1],
     chainL22,
     'uL22',
     residueInL4.get_resname(),
     residueInL4.get_id()[1],
     chainL4,
     'uL4',
     centerline))

    
    return {
        "uL22"      : residueInL22,
        "uL4"       : residueInL4,
        "centerline": centerline
    }

def added_ptc_coordinates_residues():


    for struct in log.all_structs():

        resp       = calc_constriction_site(struct)

        res22      = resp['uL22']
        res4       = resp['uL4']

        centerline = list( resp['centerline'] )

        try:
            centerline="{},{},{}".format(
            around(centerline[0],decimals=2),
            around(centerline[1],decimals=2),
            around(centerline[2], decimals=2)
            )

            log.update_struct(struct, 
            # constrictionResidueL4_id       = res4.get_id()[1],
            # constrictionResidueL22_id      = res22.get_id()[1],
            # constrictionResidueL4_resname  = res4.get_resname(),
            # constrictionResidueL22_resname = res22.get_resname(),
            constriction_coord             = str( centerline ))
        except:
            pass

        log._write()










