import os,sys
from pymol import cmd


pdbid=sys.argv[1].upper()
cmd.load(f'/home/rxz/dev/ribetl/static/{pdbid}/{pdbid}.cif')
cmd.reset()
cmd.ray(300,300)
cmd.png(f"/home/rxz/dev/ribosome.xyz-frontend.ts/public/ray_templates/_ray_{pdbid}.png")


print('Processed {}'.format(pdbid))



