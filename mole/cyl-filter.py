import numpy as np
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch





def load_tunnel(n):
    df = pd.read_csv(f"./MOLEtrials/4ug0/csv/tunnel_{n}.csv")
    x = list( df['X'].values );y = list( df['Y'].values );z = list( df['Z'].values )
    return np.array([ x,y,z ])
    
def load_tunnels():
    tunnels=[]
    for i in [1,2,3,4,5]:
        tunnels.append(load_tunnel(i))
    return tunnels

  

# exits = [[[ 159 ],[ 185 ],[ 152 ]],[[ 105 ],[ 156 ],[ 160 ]]]
exits = np.array([[159,185, 152],[105,156, 160]])


point = np.array([100, 200,300])

fuclrum = exits[1]-exits[0]
fulcrum_mag = np.linalg.norm(fuclrum)
moment = np.cross(exits[0],exits[1])

def dist(pt):
    # https://math.stackexchange.com/questions/3518495/check-if-a-general-point-is-inside-a-given-cylinder
    # distance is  norm of the cross-prod of fulcrum with ( the difference of exit 1 and the point ) divided by the norm of fulcrum
    diff = pt - exits[0]
    return np.linalg.norm( np.cross(fuclrum,diff) )/np.linalg.norm(fuclrum)


alltnls = load_tunnels()

# for tunnel in alltnls:


fig = plt.figure(figsize=(7,7))
ax  = fig.add_subplot(111, projection='3d')
for ext in exits:
    ax.plot(ext[0], ext[1],ext[2], '^', markersize=10, color='purple', alpha=1)

[ t1,t2,t3,t4,t5 ] = alltnls

ax.plot(t1[0,:], t1[1,:], t1[2,:], 'o', markersize=8, color='m', alpha=0.2)
ax.plot(t2[0,:], t2[1,:], t2[2,:], 'o', markersize=8, color='c', alpha=0.2)
ax.plot(t3[0,:], t3[1,:], t3[2,:], 'o', markersize=8, color='g', alpha=0.2)
ax.plot(t4[0,:], t4[1,:], t4[2,:], 'o', markersize=8, color='k', alpha=0.2)
ax.plot(t5[0,:], t5[1,:], t5[2,:], 'o', markersize=8, color='orange', alpha=0.2)

# fulcrum

ax.plot([ exits[0][0],exits[1][0] ],[ exits[0][1],exits[1][1] ], [ exits[0][2],exits[1][2] ])

# universal origin 
ax.plot([136],[161],[145], '.', markersize=5, color='green', alpha=1)



# t = t4
for tnum,t in enumerate(alltnls):
    mx = 0
    for i in range(t.shape[1]):
        point = t[:,i]
        d = dist(point)
        # print("POINT ", i,point," dist :", d)
        mx = d if d > mx else mx;
    print("Maximum distance {tnum}:" ,mx)

plt.show()