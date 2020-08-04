import pymol
from pymol import cmd, selector
import os
import shutil
import numpy as np


#highlight a specific chain
def see(chain_name):
    cmd.select("chain " + chain_name)
    cmd.color("gray20", "all and (not chain " + chain_name + " )")
    cmd.color("cyan", "chain " + chain_name)
    cmd.deselect()


def nbr_prompt(fulcrum, nbr):
    return "☆" if str(nbr) == str(fulcrum) else "[" + nbr + "]"


# Print neighbors of a chain within a specified radius
def get_nbrs(chain_name, radius):
    cmd.set("specular_intensity", 0.1)
    # see(chain_name)
    cmd.select("diapason", "all within {r} of chain {c}".format(
        r=radius, c=chain_name))
    #   extrapolate chains present in the selection
    cmd.select("withinChainsOf{n}".format(n=chain_name), "bc. diapason")

    neighbours = cmd.get_chains("withinChainsOf{n}".format(n=chain_name))

    cmd.color("gray10", "all")
    cmd.color("green", "chain " + chain_name)
    for nbr in neighbours:
        if nbr != chain_name:
            cmd.color("white", "chain " + nbr)
    cmd.delete("diapason")
    cmd.delete("withinChainsOf{n}".format(n=chain_name))

    print(
        chain_name + " <--> ☆ | Present in", radius,
        "Å vicinity: "
        + " ".join(map(str, [nbr_prompt(chain_name, x) for x in neighbours])))


#Selection 
def save_coords(object_name):
    cmd.create("{}".format(object_name), object_name)
    coordinate_set = cmd.get_coordset("{}".format(object_name),1)
    print("Got ", object_name, "'s coordinates: ", coordinate_set)
    print(type(coordinate_set))

    filename_string = '{name}_coordinates.npy'.format(name=object_name)
    np.save(filename_string, coordinate_set, allow_pickle=True)
    move_to_coords(filename_string)
    print("Coordinates saved successfully under ./Coordinates/{}".format(filename_string))

def move_to_coords(filename):
    cwd = os.getcwd()
    coord_path = "Coordinates"
    # Create target directory & all intermediate directories if don't exists
    try:
        os.makedirs(coord_path)
    except OSError:
        print("Failed to move.")
        pass

    shutil.move(cwd + "\\" + filename, cwd +"\\" + coord_path + "\\" + filename)


def nbrmenu():
    print("___________________________________________________________________________________")
    print("|              SIGNATURE              |              EFFECT          \n")
    print("|get_nbrs [chain_identifier],[radius] | returns array of neighbor-chain identifiers  |\n")
    print("|see               [chain_identifier] | color-highlights a chain by identifier       |\n")
    print("|save_coords            [object_name] | saves coordinates as .npy into /Coordinates|\n")
    print("|nbrmenu                              | see this tableau again                       |")
    print("___________________________________________________________________________________\n")


cmd.extend("get_nbrs", get_nbrs)
cmd.extend("see", see)
cmd.extend("save_coords", save_coords)
cmd.extend("nbrmenu", nbrmenu)

print("Commands have been added:")
nbrmenu()
