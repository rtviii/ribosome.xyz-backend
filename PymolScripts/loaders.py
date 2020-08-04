from pymol import cmd
import os
import numpy as np


print(os.listdir(os.getcwd() + "\\Coordinates\\"))

def load_coordinates(argument):
    print(os.listdir(os.getcwd() + "/Coordinates"))
    try:
        print("Attempting to load ", os.getcwd() + "\\Coordinates\\" + argument)
        loaded_array = np.load(os.getcwd() + "\\Coordinates\\" + argument, allow_pickle=True);
        print (argument + " was loaded successfully: ", loaded_array);

    except IOError:
        print ("File does not exist;")


# def load_xyz(filename):
#     xyz = np.load(os.getcwd() + "\\Coordinates\\" + filename, allow_pickle=True);
#     xyz_iter = iter(xyz)
#     print (type(xyz_iter))
#     print (xyz_iter)
#     cmd.alter_state(1, 'empty', '(x,y,z) = xyz_iter.next()', space=locals())
    # cmd.get_coordset(filename, filename[4:len(filename) - 4], 1);

load_coordinates(raw_input("Extended filename to load from (in ./Coordinates): "))
# cmd.extend("load_xyz", load_xyz)
