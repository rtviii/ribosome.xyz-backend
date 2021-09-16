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
from ciftools.TunnelScripts.ptc_site import STATIC_ROOT, get_ptc_residues, residues_average
from ciftools.Structure import fetchStructure
from ciftools.TunnelScripts.WallsReportGeneration import InitWalls, add_nomenclature_map_to_report
from ciftools.Neoget import _neoget


log=Log(os.getenv("TUNNEL_LOG"))
# pdbid=sys.argv[1].upper()

_300852 = ["1VY4", "1VY5", "1VY6", "1VY7", "4P6F", "4P70", "4TUA", "4TUB", "4TUC", "4TUD", "4TUE", "4U1U", "4U1V", "4U20", "4U24", "4U25", "4U26", "4U27", "4V90", "4W29", "4W2E", "4W2F", "4W2G", "4W2H", "4W2I", "4W4G", "4WF1", "4WOI", "4WPO", "4WQ1", "4WQF", "4WQR", "4WQU", "4WQY", "4WR6", "4WRA", "4WRO", "4WSD", "4WSM", "4WT1", "4WT8", "4WU1", "4WZD", "4WZO", "4Y4O", "4Y4P", "4YPB", "4YZV", "4Z3S", "4Z8C", "4ZER", "4ZSN", "5CZP", "5DFE", "5DOX", "5DOY", "5E7K", "5E81", "5EL4", "5EL5", "5EL6", "5EL7", "5F8K", "5FDU", "5FDV", "5HAU", "5HCP", "5HCQ", "5HCR", "5HD1", "5IB8", "5IBB", "5IMQ", "5J30", "5J3C", "5J4B", "5J4C", "5J4D", "5J8B", "5MDY", "5NDJ", "5NDK", "5OT7", "5UQ7", "5UQ8", "5VP2", "5VPO", "5VPP", "5W4K", "5WIS", "5WIT", "5ZLU", "6BUW", "6BZ6", "6BZ7", "6BZ8", "6CAE", "6CFJ", "6CFK", "6CFL", "6CZR", "6FKR", "6GSJ", "6GSK", "6GSL", "6GZQ", "6ND5", "6ND6", "6NDK", "6NSH", "6NTA", "6NUO", "6NWY", "6O3M", "6O97", "6OF1", "6OF6", "6OJ2", "6OPE", "6ORD", "6OTR", "6OXA", "6OXI", "6Q95", "6QNQ", "6QNR", "6UCQ", "6UO1", "6XQD", "6XQE", "7JQL", "7JQM"]


def inited_ids():
    for struct in _300852:
        x         = getConstrictedProteins(struct)
        log.update_struct(struct,**x,taxid=[300852], dowrite=True)

def added_constriction():
    for struct in _300852:
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
    for struct in _300852:
        
        x = residues_average( get_ptc_residues(fetchStructure(struct),struct,[ 2055 , 2056 , 2451 , 2452 , 2507 , 2506 ]) )
        if len(x)>0:
            cord = list( map(lambda x: str(x), x) )
            cord = ",".join(cord)
            log.update_struct(struct,dowrite=True,**{"PTC_average":str(cord)})
        else:continue


added_ptcs()