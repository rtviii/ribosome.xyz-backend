from pymol import cmd
import argparse
import os
from pandas import DataFrame


# Khanh's implementation
# xyz=cmd.get_coords(’Tunnel5’,1)
# r=[]
# cmd.iterate_state(1,’Tunnel’,’r.append(vdw)’,space=locals(),atomic=0)
# python
# from pymol import stored
# np.savetxt(’tunnel_coordinates.txt’,xyz,fmt=’\%.2f’)
# np.savetxt(’tunnel_radius.txt’,r,fmt=’\%.2f’)
# python end

if __name__ =='__main__':
    print(f"Executing {__file__} as a standalone script")

    parser = argparse.ArgumentParser(f"Argparser for {__file__}", add_help="""
    p3 mole/scoop.py -radius 50 -p ./MOLEtrials/4ug0/pymol/4ug0.cif -c L5  -r 4452 4UG0
""")
    parser.add_argument('-radius','-scoop_radius', dest='radius',help="the radius of matter to extract around a given residue")
    parser.add_argument('-o','--origin',dest='origin', help='origin of the sphere to scoop matter from.', required=False)
    parser.add_argument('-p,','--cifpath', help='Relative path to the cif file of the molecule to extract matter from.')
    parser.add_argument('-c','--chain')
    parser.add_argument('-r','--residue')
    parser.add_argument('pdbid', help='pdbid')
    args = parser.parse_args()  
    path = args.cifpath
    pdbid = args.pdbid.lower()
    origin = args.origin
    radius = args.radius
    chain = args.chain
    residue = args.residue




def scoop():
    print("got args \n",args)
    cmd.load(path)
    x = cmd.select(f'c. {chain} and resi {residue}')
    print(f'br. {pdbid} w. {radius} of \'sele\'')
    cmd.select(f'br. {pdbid} w. {radius} of \'sele\'')
    cmd.save( os.path.join(os.path.dirname(path), f'scoop_{pdbid}_{radius}.pdb'), 'sele' )


scoop()