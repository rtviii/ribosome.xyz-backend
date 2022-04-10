import os,sys
from pymol import cmd, util

cmd.load('4ug0.cif')
cmd.spectrum('chain')
cmd.show('surface')
cmd.set('transparency', 0.8)
cmd.reset()
cmd.ray(300,300)
cmd.png(f"RAY_4ug0.png")


