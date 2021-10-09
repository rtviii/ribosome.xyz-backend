from pprint import pprint
from typing import List, Union
from Bio.PDB import MMCIF2Dict, MMCIFIO
from Bio.PDB import  MMCIFIO as io
from Bio.PDB import  Structure
from Bio.PDB.Superimposer import Superimposer
import sys,os



def get_dict(path:str,)->dict:
	return MMCIF2Dict.MMCIF2Dict(path)


def _get_strand_ids(mmcifdict:dict)->List[str]:
	return mmcifdict['_entity_poly.pdbx_strand_id']