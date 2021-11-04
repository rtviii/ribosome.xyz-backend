from pymol import cmd
from ribetl.ciftools.super import ranged_super





@cmd.extend
def rsuper(src_struct,tgt_struct,poly_class,rstart,rend):
	cmd.delete('all')

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