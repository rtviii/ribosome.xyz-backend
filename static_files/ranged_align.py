import argparse
import os

from dotenv import load_dotenv
from pymol import cmd
from ribetl.ciftools.super import ranged_super



if __name__ =="__main__":
	load_dotenv(dotenv_path='/home/rxz/dev/riboxyzbackend/rxz_backend/.env')
	STATIC_ROOT = os.environ.get('STATIC_ROOT')

prs = argparse.ArgumentParser()

prs.add_argument('-s', '--source_struct' , type=str, required=True)
prs.add_argument('-t', '--target_struct' , type=str, required=True)
prs.add_argument('-c', '--poly_class'    , type=str, required=True)
prs.add_argument('-r', '--residue_range' , type=str, required=True)

args = prs.parse_args()

src_struct      =            args.source_struct.upper()
tgt_struct      =            args.target_struct.upper()
poly_class      =            args.poly_class
rstart    ,rend = [* map(int,args.residue_range        .split("-")) ]

[ c1 ,r1,c2,r2 ] = ranged_super( poly_class   ,src_struct, tgt_struct,( int(rstart), int(rend )))

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
cmd.save("RANGED_ALIGNMENT.pdb")


