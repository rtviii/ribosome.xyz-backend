
from pymol import cmd
from pymol import stored
import os

# ===========================================================
# 3J7Y ---- h. sapiens ---- auth T
# 3J7Z ---- e. coli    ---- auth S

def toAlignUL22():
    cmd.reinitialize()
    p1 = '3j7y.pdb'
    p2 = '3j7z.pdb'
    
    cmd.load(p1) 
    cmd.select('c1', 'c. T')
    cmd.create('ch1', 'c1')
    cmd.save('chain1.pdb', 'ch1')
    cmd.delete('all')
    
    cmd.load(p2) 
    cmd.select('c2', 'c. S')
    cmd.create('ch2', 'c2')
    cmd.save('chain2.pdb', 'ch2')
    cmd.delete('all')

    cmd.reinitialize()
    cmd.load('chain1.pdb', 'c1')
    cmd.load('chain2.pdb', 'c2')   
    
    cmd.align('c1','c2')
    cmd.save('alignedUL22.pdb')
    
    os.remove('chain1.pdb')
    os.remove('chain2.pdb')
    
    
# ===========================================================
# 3J7Y ---- h. sapiens ---- auth F
# 3J7Z ---- e. coli    ---- auth E
    
def toAlignUL4():
    cmd.reinitialize()
    p1 = '3j7y.pdb'
    p2 = '3j7z.pdb'
    
    cmd.load(p1) 
    cmd.select('c1', 'c. F')
    cmd.create('ch1', 'c1')
    cmd.save('chain1.pdb', 'ch1')
    cmd.delete('all')
    
    cmd.load(p2) 
    cmd.select('c2', 'c. E')
    cmd.create('ch2', 'c2')
    cmd.save('chain2.pdb', 'ch2')
    cmd.delete('all')

    cmd.reinitialize()
    cmd.load('chain1.pdb', 'c1')
    cmd.load('chain2.pdb', 'c2')   
    
    cmd.align('c1','c2')
    cmd.save('alignedUL4.pdb')
    
    os.remove('chain1.pdb')
    os.remove('chain2.pdb')
        
    
# ===========================================================
# 4UG0 ---- h. sapiens ---- auth L7
# 4V6C ---- e. coli    ---- auth BB

def toAlign5S():
    cmd.reinitialize()
    p1 = '4ug0.cif'
    p2 = '4v6c.cif'
    
    cmd.load(p1) 
    cmd.select('c1', 'c. L7')
    cmd.create('ch1', 'c1')
    cmd.save('chain1.pdb', 'ch1')
    cmd.delete('all')
    
    cmd.load(p2) 
    cmd.select('c2', 'c. BB')
    cmd.create('ch2', 'c2')
    cmd.save('chain2.pdb', 'ch2')
    cmd.delete('all')

    cmd.reinitialize()
    cmd.load('chain1.pdb', 'c1')
    cmd.load('chain2.pdb', 'c2')   
    
    cmd.align('c1','c2')
    cmd.save('aligned5S.pdb')
    
    os.remove('chain1.pdb')
    os.remove('chain2.pdb')

    
    
cmd.extend("toAlignUL22", toAlignUL22)
cmd.extend("toAlignUL4", toAlignUL4)
cmd.extend("toAlign5S", toAlign5S)

