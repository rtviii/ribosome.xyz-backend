# p3 mole/driver.py   --exports=t -o_points 139,154,143  133.336,172.016,150.413 -pr 4 -it 0.9 -br 2  -input ./static/pdb-structs/scoop50-4ug0.pdb --output_path ./MOLEtrials/5tunnels
# p3 mole/driver.py   --exports=t -o_points 139,154,143  133.336,172.016,150.413 -pr 4 -it 0.9 -br 2  -input ./static/pdb-structs/scoop50-4ug0.pdb --output_path ./MOLEtrials/5tunnels

import numpy as np
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch

import seaborn as sns; sns.set()
import argparse

# Exit/entry
#  resv 4200 P atom on exit [[159.12  185.127 152.528]]
#  resv 136 CD1 atom on entry[[105.527 156.633 160.96 ]]




seedpts = [ [[ 139 ], [ 154 ], [ 143 ]], 
[[ 133.336 ], [ 172.016 ], [ 150.413 ]]]
exits = [[[159],[185],[152]],[[105],[156],[160]]]



class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs
    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

def load_tunnel(n):
    df = pd.read_csv(f"./MOLEtrials/5tunnels/csv/tunnel_{n}.csv")
    x = list( df['X'].values );y = list( df['Y'].values );z = list( df['Z'].values )
    return np.array([ x,y,z ])
    
def load_tunnels():
    tunnels=[]
    for i in [1,2,3,4,5]:
        tunnels.append(load_tunnel(i))
    return tunnels

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

# Does it make sense to take the mean of samples with different number of measurements?
# How to extract the principal component alogn a certain axis, i.e. the one defined by the specified exits?

if __name__ =='__main__':
    print(f"Executing {__file__} as a standalone script")

    # parser = argparse.ArgumentParser(f"Argparser for {__file__}")
    # parser.add_argument('f', help='csv path')
    # args = parser.parse_args()

    # t1= load_tunnel(args.f)
    alltunnels = load_tunnels()
    [t1,t2,t3,t4,t5] = alltunnels
    
    for i,tunnel in enumerate(alltunnels):
    
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

    ax.plot(t1[0,:], t1[1,:], t1[2,:], 'o', markersize=8, color='green', alpha=0.2)
    ax.plot(t2[0,:], t2[1,:], t2[2,:], 'o', markersize=8, color='red', alpha=0.2)
    ax.plot(t3[0,:], t3[1,:], t3[2,:], 'o', markersize=8, color='yellow', alpha=0.2)
    ax.plot(t4[0,:], t4[1,:], t4[2,:], 'o', markersize=8, color='blue', alpha=0.2)
    ax.plot(t5[0,:], t5[1,:], t5[2,:], 'o', markersize=8, color='black', alpha=0.2)

    for pt in seedpts:
        ax.plot(pt[0], pt[1], pt[2], 'x', markersize=4, color='black', alpha=1)
    for ext in exits:
        ax.plot(ext[0], ext[1], ext[2], '^', markersize=10, color='purple', alpha=1)

    # Origin of all tunnels
    ax.plot([136],[161],[145], '.', markersize=12, color='red', alpha=1)
    
    plt.show()



