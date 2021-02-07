#!/usr/bin/python3

import sys,os
from pymol import cmd
#from dotenv import load_dotenv

print("hey")
# Temporary file to server aligned protein 

TEMP_CHAIN='/home/rt/backend/static/_TEMP_CHAIN.cif'

handle1=sys.argv[1]
handle2=sys.argv[2]

pymol_name1=sys.argv[3]
pymol_name2=sys.argv[4]


cmd.load(handle1,pymol_name1)
cmd.load(handle2,pymol_name2)
cmd.align(pymol_name1,pymol_name2)


#print(TEMP_CHAIN)
cmd.save(TEMP_CHAIN)



