from pymol import cmd
from Bio.PDB import Structure, FastMMCIFParser,Chain,Residue
from argparse import ArgumentParser

def cli():
    parser = ArgumentParser(description='Residue-wide homology analysis. Primary focus is the exit tunnel')
    parser.add_argument('--verbose', help='To enable atomwise logging. Not useful in any way.')
    parser.add_argument('-p','--path',help="Path to file(Relative?)")
    parser.add_argument('-t','--target-protein', help="Protein to operate inside the chosen strucutre")
    parser.add_argument('-r','--radius', help='clustering radius', type=float)
    return parser.parse_args()




def fetch_molecule():
    struct = FastMMCIFParser(QUIET=True).get_structure(structure_id='4UG0',filename='4ug0.cif')
    # Find out chain nomnenclature through neo4j 
    # In this case LP for uL22 and LC for uL4
    uL22 = struct[0].child_dict['LP']
    uL4  = struct[0].child_dict['LC']

    print(uL22)
    print(uL4)











def driver():   
    args = cli()
    fetch_molecule()

if __name__=="__main__":
    print(f"""
    ⋯⋯⋅⋅⋅⋅┈⋅⋄⋅⋅⋅⋄⋯⋅⋅⋅⋅┈⋅⋄⋅⋅⋅⋄⋅⋅⋅⋅┈⋅⋄⋅⋅⋅⋄⋯⋯⋅⋅⋅⋅┈⋅
    A template module for ribosome.xyz.
    Running {__file__} as standalone module.
    ⋯⋯⋅⋅⋅⋅┈⋅⋄⋅⋅⋅⋄⋯⋅⋅⋅⋅┈⋅⋄⋅⋅⋅⋄⋅⋅⋅⋅┈⋅⋄⋅⋅⋅⋄⋯⋯⋅⋅⋅⋅┈⋅
    \n\n
    """)
    driver()