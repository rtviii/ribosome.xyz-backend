# p3 mole/driver.py   --exports=t -o_points 139,154,143  133.336,172.016,150.413 -pr 4 -it 0.9 -br 2  -input ./static/pdb-structs/scoop50-4ug0.pdb --output_path ./MOLEtrials/5tunnels
# p3 mole/driver.py   --exports=t -o_points 139,154,143  133.336,172.016,150.413 -pr 4 -it 0.9 -br 2  -input ./static/pdb-structs/scoop50-4ug0.pdb --output_path ./MOLEtrials/5tunnels

import numpy as np
import pandas as pd
import glob
from pandas import DataFrame
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch
import matplotlib.patches as mpatches

import seaborn as sns; sns.set()
import argparse

# Exit/entry
#  resv 4200 P atom on exit [[159.12  185.127 152.528]]
#  resv 136 CD1 atom on entry[[105.527 156.633 160.96 ]]

str_tunnel_map = {
    "4UG0":{
        "exit":{
            "points":np.array([[159,185, 152],[105,156, 160]]),
            "residues":[  ]
        }
    },
    "5NWY":{
        "exit":{
            "points":np.array([[180,267,186],[207,197,230]]),
            "residues":['GLY/83','A/4503']
        }
    }
}

if __name__ =='__main__':
    print(f"Executing {__file__} as a standalone script")
    parser = argparse.ArgumentParser(f"Argparser for {__file__}")
    parser.add_argument('-path', '--csvpaths', help='csv path',dest='csvspath')
    parser.add_argument('struct', help='csv path')
    args = parser.parse_args()  
    csvspath = args.csvspath
    struct = args.struct


# Does it make sense to take the mean of samples with different number of measurements?
# How to extract the principal component alogn a certain axis, i.e. the one defined by the specified exits?

exits = str_tunnel_map[struct]['exit']['points']

def load_tunnel(path):
    df = pd.read_csv(path)
    x = list( df['X'].values );y = list( df['Y'].values );z = list( df['Z'].values )
    return np.array([ x,y,z ])
    
def load_tunnels(path=csvspath):
    tunnels=[]
    print(f'got path {path}')
    csvs = glob.glob(path + '/*.csv')
    for i in csvs:
        tunnels.append(load_tunnel(i))
    return tunnels
ts = load_tunnels();


def get_mean(tunnel:np.array):
    mean_x      = np.mean(tunnel[0,:])
    mean_y      = np.mean(tunnel[1,:])
    mean_z      = np.mean(tunnel[2,:])

    return np.array([[mean_x],[mean_y],[mean_z]])

    

def get_eig_pairs(tunnel):
    # get mean 
    # --> get scatter by subtracting mean from each datapoint(column) 
    # --> sum over [ outer products of these distances with themselves ]

    meanvec = get_mean(tunnel)
    scatter_matrix = np.zeros((3,3))
    for i in range(tunnel.shape[1]):
        scatter_matrix += (tunnel[:,i].reshape(3,1) - meanvec).dot((tunnel[:,i].reshape(3,1) - meanvec).T)

    eigvals, eigvecs = np.linalg.eig(scatter_matrix)
    eig_pairs = [(np.abs(eigvals[i]), eigvecs[:,i]) for i in range(len(eigvals))]

    return eig_pairs



for i,tunnel in enumerate(ts):

    print(f"\t\tTunnel {i+1}:")
    eps = get_eig_pairs(tunnel)
    eps.sort(key=lambda x: x[0], reverse=True)

    print(f"\nEigenpairs {i+1}:")
    totaleval = eps[0][0] + eps[1][0] + eps[2][0]
    for pair in eps:
        print(pair, f"\t| normalized : { round(pair[0]/totaleval, 2) }")
    print("\n")


fig = plt.figure(figsize=(7,7))
ax  = fig.add_subplot(111, projection='3d')

for ext in exits:
    ax.plot(ext[0], ext[1],ext[2], '^', markersize=10, color='purple', alpha=1)

colors = ["b", "g", "r", "c", "m", "y", "k", "w"]

for i, tunnel in enumerate(ts):
    ax.plot(tunnel[0,:], tunnel[1,:], tunnel[2,:], '-', markersize=4, color=colors[i % len(colors)], alpha=0.8, label=f'Tunnel {i}')

# ax.plot(t1[0,:], t1[1,:], t1[2,:], 'o', markersize=8, color='green', alpha=0.2)
ax.plot([ exits[0][0],exits[1][0] ],[ exits[0][1],exits[1][1] ], [ exits[0][2],exits[1][2] ])
# universal origin 
origin = [ts[0][0][0],ts[0][1][0],ts[0][2][0]]
ax.plot(origin[0],origin[1],origin[2], 'o', markersize=7, color='green', alpha=1)

plt.legend()
plt.show()



