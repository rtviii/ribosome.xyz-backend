import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

if __name__ =='__main__':
    test = np.array([
        [1,2,1,2,3,-1],
        [2,3,5,6,-1,-2],
        [0,2,4,19, 24,-35]
    ])

    mean_x      = np.mean(test[0,:])
    mean_y      = np.mean(test[1,:])
    mean_z      = np.mean(test[2,:])
    mean_vector = np.array([[mean_x],[mean_y],[mean_z]])

        

    scatter_matrix = np.zeros((3,3))
    for i in range(test.shape[1]):
        scatter_matrix += (test[:,i].reshape(3,1) - mean_vector).dot((test[:,i].reshape(3,1) - mean_vector).T)
    covmat = np.cov(test)

    val_cov,vec_cov = np.linalg.eig(covmat)
    val_sct,vec_sct = np.linalg.eig(scatter_matrix)


    # Make a list of (eigenvalue, eigenvector) tuples
    eig_pairs = [(np.abs(val_sct[i]), vec_sct[:,i]) for i in range(len(val_sct))]

    # Sort the (eigenvalue, eigenvector) tuples from high to low
    eig_pairs.sort(key=lambda x: x[0], reverse=True)

    # Visually confirm that the list is correctly sorted by decreasing eigenvalues
    for i in eig_pairs:
        print(i)


    # reduced = np.hstack(( eig_pairs[0][1].reshape(3,1), eig_pairs[1][1].reshape(3,1) ))





    


    # fig = plt.figure(figsize=(7,7))
    # ax = fig.add_subplot(111, projection='3d')

    # ax.plot(test[0,:], test[1,:], test[2,:], 'o', markersize=8, color='green', alpha=0.2)
    # ax.plot([mean_x], [mean_y], [mean_z], 'o', markersize=10, color='red', alpha=0.5)
    # for v in vec_sct.T:
    #     a = Arrow3D([mean_x, v[0]], [mean_y, v[1]], [mean_z, v[2]], mutation_scale=20, lw=3, arrowstyle="-|>", color="r")
    #     ax.add_artist(a)
    # ax.set_xlabel('x_values')
    # ax.set_ylabel('y_values')
    # ax.set_zlabel('z_values')

    # plt.title('Eigenvectors')

    # plt.show()

