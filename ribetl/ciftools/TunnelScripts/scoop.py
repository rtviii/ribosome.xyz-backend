import os,sys
from dotenv import load_dotenv

def root_self(rootname:str='')->str:
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    ROOT=os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(ROOT)
    load_dotenv(os.path.join(ROOT,'.env'))

if __name__=="__main__":
    root_self('ribxz')


from ciftools.Structure import fetchStructure
from ciftools.scripts.splitStruct import fetchChain
from typing import List
from Bio.PDB.Residue import Residue
from pymol import cmd
from ciftools.TunnelScripts.TunnelLog import Log


RADIUS      = os.getenv('SCOOP_RADIUS')
TUNNELS     = os.getenv('TUNNELS')
STATIC_ROOT = os.getenv('STATIC_ROOT')


pdbid         = sys.argv[1].upper()
log           = Log(os.getenv('TUNNEL_LOG'))
struct        = fetchStructure(pdbid)
record        = log.get_struct(pdbid)
species       = str( int( record.taxid.values[0] ) )

structpath    = os.path.join(STATIC_ROOT,pdbid,"{}.cif".format(pdbid))
scoopsavepath = os.path.join(TUNNELS,species,pdbid,'{}_{}Ascoop.pdb'.format(pdbid,RADIUS))
if not os.path.exists(os.path.dirname( scoopsavepath )):
    os.makedirs(os.path.dirname( scoopsavepath ))


cmd.load(structpath)
x = cmd.select(f'resi 2504')
cmd.select(f'br. {pdbid} w. {RADIUS} of \'sele\'')
cmd.save( scoopsavepath, 'sele' )
print("Saved to {}".format(scoopsavepath))

