import os,sys
from pymol import cmd, util


pdbid=sys.argv[1].upper()
destination = f"/home/rxz/dev/ribosome.xyz-frontend.ts/public/ray_templates/_ray_{pdbid}.png"
# if os.path.isfile(destination):
# 	print("Exists: ", destination , ". Skipping")
# 	exit(1)
# else:
cmd.load(f'/home/rxz/dev/riboxyzbackend/ribetl/static/{pdbid}/{pdbid}.cif')
cmd.reset()
cmd.spectrum('chain')
# cmd.show('surface', '{}'.format(pdbid))
# cmd.set ( "transparency", 0.45, "{}".format(pdbid) )
cmd.ray(500,500)
cmd.png(f"/home/rxz/dev/ribosome.xyz-frontend.ts/public/ray_templates/_ray_{pdbid}.png")
print('Saved {}'.format(pdbid))



