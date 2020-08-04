import numpy as np
from Bio import PDB
import pandas as pd
import argparse
import pandas as pd


# KMEANS on... what? residues?
# figure out how to pull prot seq, homology data for proteins
# ambient saturation field to go with


def open_struct(pdbid, filepath):
    """"Returns an open structure"""
    if filepath[-4:] == '.pdb':
        parser = PDB.PDBParser(QUIET=True)
    elif filepath[-4:] == '.cif':
        parser = PDB.FastMMCIFParser(QUIET=True)
    return parser.get_structure(pdbid, filepath)


def residue_distance(residue_one, residue_two):
    """Returns the C-alpha distance between two residues"""

    # Taking the average of each atom's position in each of the two residues
    avg1 = np.average([atom.get_coord() for atom in residue_one])
    avg2 = np.average([atom.get_coord() for atom in residue_two])
    diff_vector = avg1 - avg2
    return np.sqrt(diff_vector * diff_vector)


def chain_distance_matrix(c1, c2):
    print("Getting the dist matrix between chain {} and chain {}".format(
        c1.get_id(), c2.get_id()))
    """Returns a matrix of alpha-carbon distances between two chains"""
    dist_mat = np.zeros((len(c1), len(c2)), np.float)
    for i, residue_one in enumerate(c1):
        for j, residue_two in enumerate(c2):
            dist_mat[i, j] = residue_distance(residue_one, residue_two)
    return dist_mat


def contact_map_catalogue(pdbid, cifpath):
    model = open_struct(pdbid, cifpath)[0]
    rows = {}
    # Populate the dataframe with subchains -- keys
    for i in model:
        rows[i.get_id()] = []

    contact_map_catalogue = pd.DataFrame()
    # For each subchain in the model
    for j in model:
        # Get its distance matrix to every other subchain, includign self
        jthcolumn = [chain_distance_matrix(i, j) for i in model]
        # Add as a column
        contact_map_catalogue[j.get_id()] = jthcolumn
    return contact_map_catalogue


x = contact_map_catalogue('1xi4', './1xi4.pdb')
print(x)


# think about how to inject ribosomal sequence into the data.
# residue-wise


# # STRUCTURE --> MODEL -->  CHAINS  --> RESIDUES --> ATOMS
