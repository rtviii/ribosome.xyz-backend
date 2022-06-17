import datetime
import os
import sys
from pymol import cmd



# super.ranged_super aligns chains via biopython's seqalign and returns the correct individual ranges

# os.system(f"""echo \"Last called ranged_align.py at {str( datetime.datetime.now() ).split(" ")[0]} with args :
# 	src_struct       = {sys.argv[1].upper()}
# 	tgt_struct       = {sys.argv[2].upper()}
# 	src_auth_asym_id = {sys.argv[3]}
# 	tgt_auth_asym_id = {sys.argv[4]}
# 	r1start , r1end 	 = {sys.argv[5]}
# 	r2start , r2end 	 = {sys.argv[6]}
# 	\" >> alignment{datetime.datetime.now()}.log""")

# os.system(f"""echo \"
# ranged_super.py returned 
# 	paths     : 
# 	chain1    : {c1},
# 	chain2    : {c2},
# 	source_range: [ {r1start}, {r1end} ],
# 	target_range: [ {r2start}, {r2end} ],
# 	\" >> alignment{str(datetime.datetime.now()).split(" ")[0]}.log""")

# parse args with pain
src_struct = sys.argv[1].upper()
tgt_struct = sys.argv[2].upper()
src_auth_asym_id = sys.argv[3]
tgt_auth_asym_id = sys.argv[4]
r1start, r1end = [* map(int, sys.argv[5].split("-"))]
r2start, r2end = [* map(int, sys.argv[6].split("-"))]

# ?get chain paths (in static) and correct ranged to clip out
# ?[c1, r1, c2, r2] = super.ranged_super(src_struct, src_auth_asym_id, tgt_struct, tgt_auth_asym_id, (int(rstart), int(rend)))

c1 = os.path.join(os.environ.get( "STATIC_ROOT" ), src_struct.upper(), "CHAINS", "{}_STRAND_{}.cif".format(src_struct.upper(),src_auth_asym_id))
c2 = os.path.join(os.environ.get( "STATIC_ROOT" ), tgt_struct.upper(), "CHAINS", "{}_STRAND_{}.cif".format(tgt_struct.upper(),tgt_auth_asym_id))

# log
n1 = c1.split("/")[-1].split('.')[0]
n2 = c2.split("/")[-1].split('.')[0]

# Clip chains with pymol, create snippet objects, align those and save.
cmd.load(c1)
cmd.select("resi {}-{} and m. {} ".format(r1start,r1end, n1))
cmd.create("{}_{}".format(src_struct, src_auth_asym_id), "sele")
cmd.delete(n1)

cmd.load(c2)
cmd.select("resi {}-{} and m. {} ".format(r2start,r2end, n2))
cmd.create("{}_{}".format(tgt_struct, tgt_auth_asym_id), "sele")
cmd.delete(n2)

cmd.super("{}_{}".format(src_struct, src_auth_asym_id),
          "{}_{}".format(tgt_struct, tgt_auth_asym_id))

cmd.save(os.environ.get("TEMP_CHAIN"))

