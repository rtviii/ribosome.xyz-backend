#!/usr/bin/python3

import sys,os
from pymol import cmd
from datetime import date,datetime

# Temporary file to serve aligned protein 



os.system(f"""echo \"Last called align.py at {datetime.now()} with args :
handle1    ={sys.argv[1]}
handle2    ={sys.argv[2]}
pymol_name1={sys.argv[3]}
pymol_name2={sys.argv[4]}
\" > alignment.log""")

handle1    =sys.argv[1]
handle2    =sys.argv[2]

pymol_name1=sys.argv[3]
pymol_name2=sys.argv[4]

TEMP_CHAIN = os.getenv("TEMP_CHAIN")

cmd.load (handle1    ,pymol_name1)
cmd.load (handle2    ,pymol_name2)
cmd.super(pymol_name1,pymol_name2)
print("Aligning using SUPER")
cmd.save(TEMP_CHAIN)

print(f"Saved to {TEMP_CHAIN}")


