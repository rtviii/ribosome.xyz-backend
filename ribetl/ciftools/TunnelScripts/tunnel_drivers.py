#!/usr/bin/env python3
from ciftools.Chain import matchClassToStrand
from ciftools.Neoget import _neoget
import os,sys
import pandas as pd
from dotenv import load_dotenv

def root_self(rootname:str='')->str:

    """Returns the rootpath for the project if it's unique in the current folder tree."""
    ROOT = os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(ROOT)
    load_dotenv(os.path.join(ROOT,'.env'))

if __name__=="__main__":
    root_self('ribxz')


from ciftools.Structure import fetchStructure
from ciftools.TunnelScripts.TunnelLog import (Log)
from ciftools.TunnelScripts.ptc_site import (residues_average, get_ptc_residues)
from ciftools.TunnelScripts.WallsReportGeneration import InitWalls
from ciftools.TunnelScripts.WallsReportGeneration import add_nomenclature_map_to_report
TUNNELS     = os.getenv("TUNNELS")
STATIC_ROOT = os.getenv("STATIC_ROOT")
log         = Log(os.getenv('TUNNEL_LOG'))


# for struct in log.all_structs():
#     pdbid='3J7Z'
#     REPORT_PATH=os.path.join(STATIC_ROOT,pdbid,"{}_TUNNEL_REPORT.json".format(pdbid))
#     ## Generating report
#     walls =InitWalls(pdbid)
#     walls.consumeMoleDataframe(10)
#     walls.generateReport(REPORT_PATH)
#     add_nomenclature_map_to_report(pdbid,REPORT_PATH)

## Adding ptc
# for struct in log.all_structs():
#     log.update_struct(struct,
#     dowrite=True,
#     **{"PTC_average":
#         ",".join(map(lambda x:str(x), list( residues_average( get_ptc_residues(fetchStructure(struct),struct)) ) ) ) 
#         })