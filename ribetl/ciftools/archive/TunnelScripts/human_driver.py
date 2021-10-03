#!/usr/bin/env python3

import json
import os,sys
from numpy.core.numeric import NaN
import pandas as pd
from dotenv import load_dotenv

def root_self(rootname:str='')->str:

    """Returns the rootpath for the project if it's unique in the current folder tree."""
    ROOT = os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(ROOT)
    load_dotenv(os.path.join(ROOT,'.env'))

if __name__=="__main__":
    root_self('ribxz')


from ciftools.TunnelScripts.TunnelLog import (Log, TunnelWalls, get_CA_or)
from ciftools.TunnelScripts.constriction_site import ( getConstrictedProteins, calc_constriction_site,added_ptc_coordinates_residues )
from ciftools.TunnelScripts.ptc_site import STATIC_ROOT, get_ptc_residues
from ciftools.Structure import fetchStructure
from ciftools.TunnelScripts.WallsReportGeneration import InitWalls, add_nomenclature_map_to_report
from ciftools.Neoget import _neoget


log=Log(os.getenv("TUNNEL_LOG"))
# pdbid=sys.argv[1].upper()

_9606   = ["3J9M", "3J7Y", "3J92", "3JAG", "3JAH", "3JAI", "4UG0", "5A2Q", "5AJ0", "5LKS", "5LZT", "5LZU", "5LZV", "5LZW", "5LZX", "5LZY", "5LZZ", "5OOL", "5OOM", "5T2C", "6D90", "6D9J", "6EK0", "6G18", "6G5H", "6G5I", "6GAW", "6GAZ", "6GB2", "6HCF", "6I0Y", "6IP5", "6IP8", "6LQM", "6LSR", "6LSS", "6LU8", "6NU2", "6OLE", "6OLF", "6OLG", "6OLI", "6OLZ", "6OM0", "6OM7", "6P5N", "6QZP", "6R5Q", "6R6G", "6R6P", "6R7Q", "6RW4", "6RW5", "6T59", "6W6L", "6XA1", "6Y0G", "6Y2L", "6Y57", "6Z6L", "6Z6M", "6Z6N", "6ZVH", "7K5I"]    


ptcResidue = '4452'


def inited_ids():
    for struct in _9606:
        x         = getConstrictedProteins(struct)
        log.update_struct(struct,dowrite=True,**x)

def added_constriction():
    for struct in _9606:
        s22 = log.get_struct(struct)['uL22'].values[0]
        s4  = log.get_struct(struct)['uL4'].values[0]

        if NaN in [s4,s22]:
            continue
        resp   = calc_constriction_site(struct)
        resn22 = resp['uL22'].get_resname()
        resn4  = resp['uL4'].get_resname()
        resi22 = resp['uL22'].get_id()[1]
        resi4  = resp['uL4'].get_id()[1]

        center = list( map(lambda x: str(x), resp['centerline'] ) )
        center = ",".join(center)
        update={ 
        "constrictionResidueL22_resname": resn22,
        "constrictionResidueL4_resname" : resn4,
        "constrictionResidueL22_id"     : resi22,
        "constrictionResidueL4_id"      : resi4,
        "constriction_coord"            : center
        }
        log.update_struct(struct,dowrite=True,**update)

def added_ptcs():
    for struct in _9606:
        x = get_ptc_residues(fetchStructure(struct),struct,['4452'])
        if len(x)>0:
            cord = list( map(lambda x: str(x), get_CA_or( x[0] ).get_coord() ) )
            cord = ",".join(cord)
            log.update_struct(struct,dowrite=True,**{"PTC_average":str(cord)})
        else:continue

added_constriction()
