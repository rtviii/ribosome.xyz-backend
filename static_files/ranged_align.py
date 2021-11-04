import datetime
import os
import sys
from dotenv import load_dotenv
from pymol import cmd


load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
STATIC_ROOT = os.environ.get('STATIC_ROOT')
def root_self(rootname: str = ''):
	"""Returns the rootpath for the project if it's unique in the current folder tree."""
	root = os.path.abspath(__file__)[:os.path.abspath(
		__file__).find(rootname)+len(rootname)]
	sys.path.append(root)

root_self('riboxyzbackend')
import ribetl.ciftools.super as super

os.system(f"""echo \"Last called ranged_align.py at {datetime.datetime.now()} with args :
	src_struct       = {sys.argv[1].upper()}
	tgt_struct       = {sys.argv[2].upper()}
	src_auth_asym_id = {sys.argv[3]}
	tgt_auth_asym_id = {sys.argv[4]}
	rstart , rend 	 = [* map(int,{ sys.argv[5] }.split("-")) ]

	\" > alignment.log""")

src_struct       = sys.argv[1].upper()
tgt_struct       = sys.argv[2].upper()
src_auth_asym_id = sys.argv[3]
tgt_auth_asym_id = sys.argv[4]
rstart ,                 rend = [* map(int,sys.argv[5].split("-")) ]

[ c1 ,r1,c2,r2 ] = super.ranged_super( src_struct, src_auth_asym_id, tgt_struct,tgt_auth_asym_id, ( int(rstart), int(rend )))

n1 = c1.split("/")[-1].split('.')[0]
n2 = c2.split("/")[-1].split('.')[0]

cmd.load(c1)
cmd.select("resi {}-{} and m. {} ".format(*r1,n1))
cmd.create("snip1","sele")
cmd.delete(n1)

cmd.load(c2)
cmd.select("resi {}-{} and m. {} ".format(*r2,n2))
cmd.create("snip2","sele")
cmd.delete(n2)
cmd.super("snip1","snip2")
cmd.save(os.path.join(STATIC_ROOT,"RANGED_ALIGNMENT.pdb"))



