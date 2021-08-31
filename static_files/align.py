#!/usr/bin/python3

import sys,os
from pymol import cmd

# Temporary file to serve aligned protein 

handle1    =sys.argv[1]
handle2    =sys.argv[2]

pymol_name1=sys.argv[3]
pymol_name2=sys.argv[4]

TEMP_CHAIN = os.getenv("TEMP_CHAIN")

cmd.load (handle1    ,pymol_name1)
cmd.load (handle2    ,pymol_name2)
cmd.align(pymol_name1,pymol_name2)

cmd.save(TEMP_CHAIN,state=1,partial=1)



