from Bio import PDB
import pickle
import pandas as pd
import os
import numpy as np
import json
import functools
import argparse



# PATH TO CIF FILES (acquired from pdb)
directory = "./cif_models/"

def pickle_coordinates(pdbid, partition):
    if not os.path.exists('coordinate_dicts'):
        try:
            os.makedirs('coordinate_dicts')
        except:
            print("Failed to create deposition directory. Exiting")

    picklename = "./coordinate_dicts/{}.pkl".format(pdbid)
    output     = open(picklename, 'wb')
    pickle.dump(partition, output)

    print("Saved successfully to {}/coordinate_dicts/{}.pkl\n".format(os.getcwd(), pdbid))
    output.close()

    
def save_molecule_as_csv(pdbid, partition):
    print("Length of partition: ", len(partition))
    print("Parsing to .csv ...")

    df = pd.DataFrame({key: pd.Series(value) for key, value in partition.items()})

    for kvpair in partition:
        if not os.path.exists("./coordinate_csvs/"):
            try:
                os.makedirs("./coordinate_csvs/")
            except:
                print("Failed to create deposition directory(check permissions, perhaps). Exiting\n")
    outpath = './coordinate_csvs/{}.csv'.format(pdbid, pdbid)
    df.to_csv(outpath, index=True)
    print("Saved coordinates to {}.csv\n".format(pdbid))




# Usage: call with the id of the structure
# Ex. python3 extract_coordinates.py .cif_files/5jvg.cif
# Deposits a dictionary of subchains and their corresponding coordinates to /coordinate_dicts

def load_cif_from_file():

    # Parse the pdbid argument
    parser   = argparse.ArgumentParser()
    parser.add_argument('id', type=str,help='pdbid')
    parser.add_argument('format', type=str, help='csv/json/pickle')
    parser.add_argument('coord_type', type=str,help='nparray or plain list')
    parser.add_argument('precision',type=int, help="number of decimals on an atomic coordinate")
    args        = parser.parse_args()
    pdbid       = args.id
    parseformat = args.format
    coord_type  = args.coord_type
    precision   = args.precision
    filepath    = directory + pdbid + ".cif"

    # Initializing pdb parsermole
    cifparser = PDB.FastMMCIFParser(QUIET=True)
    try:
        with open(filepath) as infile:
            structure = cifparser.get_structure(pdbid, infile)
    except:
        print("Failed to open {}".format(filepath))

    print("Processing {}".format( structure))
    # Decomposing the structure by chain
    structure_atoms  = PDB.Selection.unfold_entities(structure, "A")
    structure_coords = [atom.get_coord() for atom in structure_atoms]
    allchains        = [chain for chain in structure.get_chains()]
    chainids         = [chain.get_id() for chain in allchains]

    # chain names and atom objects contained therein
    chain_atoms = [(chain.get_id(), [* chain.get_atoms()]) for chain in allchains]

    def atomArrToCord(array_of_atoms):
        return list(map(lambda atom:
                        atom.get_coord(), array_of_atoms)
                    )
    chainCoordinates = list(map((lambda chainAtomTuple: (chainAtomTuple[0], atomArrToCord(
        chainAtomTuple[1]))), chain_atoms))

    chain_coordinates_object = {}
    print("COORD TYPE {}".format(coord_type))

    for namecoordpair in chainCoordinates:
        if coord_type=='list':
            chain_coordinates_object[namecoordpair[0]] = [[ float('%.{}f'.format(precision) % coord) 
            for coord in x.tolist() ] for x in  namecoordpair[1]]
        elif coord_type=='nparray':
            chain_coordinates_object[namecoordpair[0]] = namecoordpair[1]
        else:
            print('Unknown coordinate type. Possible: [ list | nparray ]')
            print('Exiting..')
            exit

    # Sanity check
    for chain in chainids:
        if chain not in [* chain_coordinates_object.keys()]:
            print("{} is missing!".format(chain))

    if parseformat == 'csv':
        save_molecule_as_csv(pdbid, chain_coordinates_object)
    elif parseformat =='pikle': 
        pickle_coordinates(pdbid, chain_coordinates_object)
    elif parseformat == 'json':
        with open(pdbid+'.json','w') as out:
            json.dump(chain_coordinates_object, out)


load_cif_from_file()
