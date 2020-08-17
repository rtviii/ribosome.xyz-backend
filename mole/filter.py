import numpy as np
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch


import seaborn as sns; sns.set()
from sklearn.decomposition import PCA
import argparse

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



if __name__ =='__main__':
    print(f"Executing {__file__} as a standalone script")

    parser = argparse.ArgumentParser(f"Argparser for {__file__}")
    parser.add_argument('f', help='csv path')
    args = parser.parse_args()

    t1= load_tunnel(args.f)
    [t1,t2,t3,t4,t5] = load_tunnels()


    mean_x      = np.mean(t1[0,:])
    mean_y      = np.mean(t1[1,:])
    mean_z      = np.mean(t1[2,:])
    mean_vector = np.array([[mean_x],[mean_y],[mean_z]])

    scatter_matrix = np.zeros((3,3))
    for i in range(t1.shape[1]):
        scatter_matrix += (t1[:,i].reshape(3,1) - mean_vector).dot((t1[:,i].reshape(3,1) - mean_vector).T)
    print(scatter_matrix)

    eig_val_sct, eig_vec_sct = np.linalg.eig(scatter_matrix)

    eig_pairs = [(np.abs(eig_val_sct[i]), eig_vec_sct[:,i]) for i in range(len(eig_val_sct))]

    # Sort the (eigenvalue, eigenvector) tuples from high to low
    eig_pairs.sort(key=lambda x: x[0], reverse=True)
    # Visually confirm that the list is correctly sorted by decreasing eigenvalues
    for i in eig_pairs:
        print(i)

    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111, projection='3d')
# 
    ax.plot(t1[0,:], t1[1,:], t1[2,:], 'o', markersize=8, color='green', alpha=0.2)
    ax.plot(t2[0,:], t2[1,:], t2[2,:], 'o', markersize=8, color='red', alpha=0.2)
    ax.plot(t3[0,:], t3[1,:], t3[2,:], 'o', markersize=8, color='yellow', alpha=0.2)
    # ax.plot(t4[0,:], t4[1,:], t4[2,:], 'o', markersize=8, color='blue', alpha=0.2)
    # ax.plot(t5[0,:], t5[1,:], t5[2,:], 'o', markersize=8, color='black', alpha=0.2)

    ax.plot([mean_x], [mean_y], [mean_z], 'o', markersize=10, color='red', alpha=0.5)
    # for v in eig_vec_sct.T:
        # a = Arrow3D([mean_x, v[0]], [mean_y, v[1]], [mean_z, v[2]], mutation_scale=20, lw=3, arrowstyle="-|>", color="r")
        # ax.add_artist(a)
    ax.set_xlabel('x_values')
    ax.set_ylabel('y_values')
    ax.set_zlabel('z_values')

    plt.title('Eigenvectors')
    plt.show()

    print(eig_vec_cov)


