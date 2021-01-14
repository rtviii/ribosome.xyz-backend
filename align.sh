#!/bin/python3


import sys,os
from pymol import cmd



print(sys.argv[1])
print(sys.argv[2])

cmd.load("5MLC_STRAND_L.cif",'c1')
cmd.load("6LSR_STRAND_D.cif",'c2')
cmd.align("c1","c2")
cmd.save("/home/rxz/another.cif")
