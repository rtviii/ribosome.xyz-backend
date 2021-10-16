from abc import abstractproperty
from ast import arg, parse
from logging import log
import os,sys
from pprint import pprint
import json
from typing import List
from Bio import Align
from Bio.SeqRecord import SeqRecord
import dotenv
from pymol import cmd
import argparse
import glob
from ciftools.bsite_mixed import BindingSite
from Bio.Align import MultipleSeqAlignment
import numpy as np
from Bio import pairwise2
# import matplotlib as plt
import pylab as plt

sys.path.append('/home/rxz/dev/riboxyzbackend/')
dotenv.load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get("STATIC_ROOT")
nomid       = sys.argv[1]
bsites      = []

profiles = []
for  _ in os.listdir(STATIC_ROOT):
	if len(_) > 4:
		continue
	else:
		profiles = [*profiles, os.path.join(STATIC_ROOT, _, "{}.json".format(_))]

for k in profiles:
	with open(k, 'rb') as infile: 
		d = json.load(infile)
	polys  =[*d['rnas'], *d['proteins']]
	for chain in polys:
		if nomid in chain['nomenclature']:
			print("{} : {}".format(nomid,chain['entity_poly_seq_one_letter_code']))