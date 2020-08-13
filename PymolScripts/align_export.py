import ast
from pymol import cmd, selector
from argparse import ArgumentParser

def modelChainTuple(string):
    try:
        model,chain=map(str, string.split(','))
        return model,chain
    except:
        print("Type error when reading model-chain tuples")
        raise TypeError

parser = ArgumentParser(description='Specify models.chains to align or extend pymol api',)
parser.add_argument('--model-chain-tuples', nargs='+', type=modelChainTuple)
parser.add_argument('--pmlxtnd')
args =parser.parse_args()
print(args)
"""
-fetch the structures(arg) 
-parseout the chains of interest(arg) 
-align in the same coordinate frame
-save as a pdb
"""

def chain_align_save(*args, **kwargs): 
    args = [ast.literal_eval(kvpair) for kvpair in args]
    for pair in args:
        cmd.fetch(str.lower(pair[0]))
        create_subchain_object(pair[0], pair[1])
        cmd.delete(str.lower(pair[0]))
    
    for mobile in [ '{}.{}'.format(model, chain) for (model, chain) in args ][1:]:
        cmd.super(mobile, "{}.{}".format(args[0][0], args[0][1]),reset=1,transform=1,quiet=0)

    cmd.reset()
    cmd.save('alignment.cif')
    

def create_subchain_object(pdbid, subchain):
    tempname = 'selection_{}.{}'.format(pdbid, subchain)
    cmd.select(tempname ,'m. {} and c. {}'.format( pdbid, subchain))
    cmd.create('{}.{}'.format(pdbid,subchain), tempname)
    cmd.delete(tempname)

if args.pmlxtnd:
    cmd.extend('create_subchain_object', create_subchain_object)
    cmd.extend('chain_align_save', chain_align_save);
    print("Added two commands to pymol.cmd")
else:
    mcpairs=args.model_chain_tuples
    for pair in mcpairs:
        cmd.fetch(str.lower(pair[0]))
        create_subchain_object(pair[0], pair[1])
        cmd.delete(str.lower(pair[0]))
    
    for mobile in [ '{}.{}'.format(model, chain) for (model, chain) in mcpairs ][1:]:
        cmd.super(mobile, "{}.{}".format(mcpairs[0][0], mcpairs[0][1]),reset=1,transform=1,quiet=0)

    cmd.reset()
    cmd.save('alignment.cif')