import os, sys
import xml
from pymol import cmd 
from config import make_input_config 


# H. sappiens
# 4UG0 Chain L5
# A`4548/N1
# A`3908/N1
# C`2794/N3

# E. coli
# 4V9D Chain CA
# A 2602
# A 2062
# A 1614


sys.path.append(os.path.join(sys.path[0],'mole2'))
sys.path.append(sys.path[0])
path_moleexe   = os.path.join(sys.path[0],'mole2','mole2.exe')
path_moleconfig  = os.path.join(sys.path[0],'input.xml')
# path_pdbstruct = os.path.join(os.path.dirname(sys.path[0]), 'static','pdb-structs','4ug0.pdb')

path_pdbstruct = os.path.join(os.path.dirname(sys.path[0]), 'static','pdb-structs','radiusobject.pdb')

# print("MOLE:    ", path_moleexe)
# print("INPUT:   ", path_moleconfig)
# print("STRUCT:  ", path_pdbstruct)



if __name__ == "__main__":
    os.system('mono {} {}'.format(path_moleexe, path_moleconfig))
    # make_input_config()






