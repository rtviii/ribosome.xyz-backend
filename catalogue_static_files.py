#!/bin/usr/python3


import sys, os,json
import re
from pathlib import Path
from typing import Pattern




# checking for:
#! LIGANDS
#! TUNNELS REPORT
#! TUNNEL CENTERLINE


STATIC = sys.argv[1]

report = {
}
tunnels     = 0
ligands     = 0
centerlines = 0
structstotal =0
for struct in os.listdir(STATIC):
    structstotal+=1
    report[struct] = {}
    # This is the GLOB that im looking for.
    for item in ['LIGAND', 'TUNNEL_REPORT', 'CENTERLINE']:
        # the matches in the struct directories
        results  = [* map(lambda x: x.name, [*Path(os.path.join(STATIC,struct)).rglob("*"+item+"*") ]) ]
        if item == 'CENTERLINE':
            if len(results) > 0 :
                report[struct][item] = True 
                centerlines+=1
            else:
                report[struct][item] = False
        if item == 'TUNNEL_REPORT' :
            if len(results) > 0:
            # print(f"Found tunnel reports for {struct} : {results}")
                report[struct][item] = True
                tunnels+=1
            else:
                report[struct][item] = False

        if item == 'LIGAND':
            if len(results) > 0:
                def getChemid(filename:str):
                    pattenr = r'(?<=\_)(.*?)(?=\.)'
                    return re.search(pattenr,filename)[0]
                report[struct][item] = results
                ligands+=len(results)
            else:
                report[struct][item] = []

print("""
Total  Structures : {}
       Ligands    : {}
Tunnel Reports    : {}
Tunnel Centerlines: {}
""".format(
    structstotal,ligands,tunnels,centerlines
))
with open("static_files_catalogue.json", 'w') as outfile:
    json.dump(report, outfile)
