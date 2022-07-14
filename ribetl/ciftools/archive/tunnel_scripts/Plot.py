from logging import error
import os,sys,json
import numpy as np
from typing import Iterator, List
import pandas as pd
from dotenv import load_dotenv
import matplotlib
import matplotlib.pyplot  as plt
import matplotlib.patches  as patches

def root_self(rootname:str='')->str:
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    ROOT = os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(ROOT)
    load_dotenv(os.path.join(ROOT,'.env'))

root_self('ribxz')

from  ciftools.TunnelScripts.TunnelLog import (Log,get_CA_or)
from  ciftools.TunnelScripts.WallsReportGeneration import InitWalls, add_nomenclature_map_to_report


PDBID                  = sys.argv[1].upper()
STATIC_ROOT            = os.getenv("STATIC_ROOT")
STRUCT_REPORT_PATH     = os.path.join(STATIC_ROOT,PDBID,f'{PDBID}_TUNNEL_REPORT.json')

if not os.path.exists(STRUCT_REPORT_PATH):
    try:
        ## Generating report
        walls = InitWalls(PDBID)
        walls.consumeMoleDataframe(10)
        walls.generateReport(STRUCT_REPORT_PATH)
        add_nomenclature_map_to_report(PDBID,STRUCT_REPORT_PATH)
    except error:
        print(error)
        exit("Failed to generate report for {}. Does the tunnel csv exist?".format(PDBID))

log                    = Log(os.getenv('TUNNEL_LOG'))
ProteinsColorgenerator = iter( ['purple','cyan','orange',"gray",'yellow','magenta','brown','turqoise','pink'] )

with open(STRUCT_REPORT_PATH, 'rb') as infile:
    report =  json.load(infile)

walls       = InitWalls(PDBID)
tunneldf    = walls.mole_dataframe

def locate_res_in_df(rescoords:List[float],tunnel_df:pd.DataFrame): 
    respos   = np.array(rescoords)
    rowid    = -1
    currdist = 10000000

    for index,row in tunnel_df.iterrows():

        centerlinepos = np.array( ( row['X'],row['Y'], row['Z'] ) )
        dist = np.linalg.norm(centerlinepos-respos)
        if dist <  currdist:
            currdist = dist
            rowid    = index
    
    return tunnel_df.loc[rowid]

proteinsLegend = []
proteins       = [* filter(lambda x: x[1]['type']=='Protein',report['nomenclatureMap'].items()) ]
# print(proteins)


radii     = tunneldf['FreeRadius'].values
distances = tunneldf['Distance'].values

plt.plot(distances,radii, linewidth=0.7, c='gray')


residues_legend=[]
for prot in proteins:
    residuecount_legend={}
    strandName = prot[0]
    banName    = prot[1]['nomenclature'] if len( prot[1]['nomenclature'] ) >0 else "Unidentified"

    curr_color='white'
    if 'uL4' in banName :
        curr_color='green'
    if 'uL22' in banName:
        curr_color='blue'
    elif 'uL22' not in banName and 'uL4' not in banName:
        curr_color = next(ProteinsColorgenerator)

    for res in report['adjacent_strands'][strandName]:

        if res['resname']=='HOH':
            continue
        if res['resname'] not in residuecount_legend.keys(): 
            residuecount_legend[res['resname']] = 1
        else:
            residuecount_legend[res['resname']] += 1

        ####
        coordinate = res['rescoord']
        row        = locate_res_in_df(coordinate, tunneldf)
        coordinate = ( row['Distance'],row['FreeRadius'])

        if res['resname'] in ["ARG",'LYS']:
            marker = "*"
            size   = 20
        if res['resname'] in ['GLU','ASP']:
            marker = "x"
            size   = 20
        elif res['resname'] not in ['GLU','ASP','ARG','LYS']:
            size=8
            marker="v"

        resdp = { 
                "coordinate": ( row['Distance'],row['FreeRadius']),
                "marker"    : marker,
                "size"      : size,
        }
        plt.scatter(resdp['coordinate'][0],resdp['coordinate'][1], marker=resdp['marker'],s=resdp['size'], c=curr_color)


    residues_legend.append(patches.Patch(color=curr_color,label="{}({}):{}".format(
        prot[0],banName," | ".join([* map( lambda item: "#{}={}".format(item[0],item[1]),residuecount_legend.items() ) ]))))
    print(residuecount_legend)
    proteinsLegend.append(patches.Patch(color=curr_color,label=f"{strandName}({banName})"))

################### ADDING THE PTC

#ecoli            = [ 2055 , 2056 , 2451 , 2452 , 2507 , 2506 ]
#human            = [ 4452 ]

ptcresidues       = walls.get_ptc_residues([ 2055 , 2056 , 2451 , 2452 , 2507 , 2506 ])
ptc_res_locs      = [*map(lambda res: get_CA_or(res).get_coord().tolist(),ptcresidues) ]
ptc_corresp_radii = [*map(lambda respos: locate_res_in_df(respos, tunneldf), ptc_res_locs)]
ptc_corresp_radii = [*map(lambda row: ( row['Distance'],row['FreeRadius'] ),ptc_corresp_radii )]

dist              = [x[0]for x in ptc_corresp_radii]
radii             = [x[1]for x in ptc_corresp_radii]

plt.scatter(dist,radii, c='black',marker='v')

################### ADDING CONSTRICTION SITE
constrsite        = log.get_struct(PDBID)['constriction_coord'].values[0]

constrrow         = locate_res_in_df(np.array([* map(lambda x:float(x),constrsite.split(',')) ]), tunneldf)
constriction_dist = constrrow['Distance']
constriction_rad  = constrrow['FreeRadius']

plt.scatter(constriction_dist,constriction_rad,marker='^',s=120, c='red')

axs = plt.axes()
residue_summary=axs.legend(
handles=[
patches.Rectangle((0,0),1,1,fc='w',fill=False,edgecolor=None,linewidth=0,label="ARG, LYS marked with star; ASP, GLU marked with cross."),
patches.Patch(color='red',label='Contstriction;uL22/uL4 closure' ),
*residues_legend ],
fontsize=4,
loc=4)

# plt.gca().add_artist(legendproper)
plt.gca().add_artist(residue_summary)


axs.set_ylabel(ylabel="FreeRadius",fontsize=12)
axs.set_xlabel(xlabel="Distance from PTC",fontsize=12)
axs.set_title(PDBID)
figure = plt.gcf()
figure.set_size_inches(10, 4)

# plt.savefig("sample.png", dpi=100)
if sys.argv[2]=='save':
    plt.savefig('{}_radiusplot.png'.format(PDBID), dpi=1200)
if sys.argv[2]=='show':
    plt.show()

    
