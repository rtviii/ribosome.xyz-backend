
from typing import List
from Bio.PDB.Structure import Structure
import numpy as np
from Bio.PDB.Residue import Residue
import os, sys,dotenv
def root_self(rootname:str='')->str:
    """Returns the rootpath for the project if it's unique in the current folder tree."""
    root=os.path.abspath(__file__)[:os.path.abspath(__file__).find(rootname)+len(rootname)]
    sys.path.append(root)
    dotenv.load_dotenv(os.path.join(root,'.env'))

if __name__=="__main__":
    root_self('ribxz')


from ciftools.TunnelScripts.TunnelLog import (Log)
from ciftools.TunnelScripts.TunnelLog import get_CA_or
TUNNELS     = os.getenv("TUNNELS")
STATIC_ROOT = os.getenv("STATIC_ROOT")
log         = Log(os.getenv('TUNNEL_LOG'))




def get_ptc_residues(struct:Structure, pdbid:str,conserved_nucleotides:List[int])->List[Residue]:

    ECOLI_PTC_CONSERVED_NUCLEOTIDES=['2055','2451','2452','2504','2505','2506','2507']

    def belongs_to_ptc(x:Residue):
        return int(x.get_id()[1]) in conserved_nucleotides
    PTC_residues = filter(belongs_to_ptc, [*struct.get_residues()]) 
    return [* PTC_residues ]


def residues_average(residues:List[Residue], v=False):
    """Returns an average coordinate of a list of resdues by averaging over alpha carbon's or else first atom' in the iterable positions. """
    ptcres = [*map(lambda x: get_CA_or(x).get_coord() ,residues) ]
    if v:
        for res in residues:
            print(res.get_resname(), res.get_id()[1], get_CA_or(res))
    return np.average(ptcres,0)