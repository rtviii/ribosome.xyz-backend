import numpy as np
import argparse
import glob
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch




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
fuclrum     = exits[1]-exits[0]
fulcrum_mag = np.linalg.norm(fuclrum)
moment      = np.cross(exits[0],exits[1])
def dist(pt):
    # https://math.stackexchange.com/questions/3518495/check-if-a-general-point-is-inside-a-given-cylinder
    # distance is  norm of the cross-prod of fulcrum with ( the difference of exit 1 and the point ) divided by the norm of fulcrum
    diff = pt - exits[0]
    return np.linalg.norm( np.cross(fuclrum,diff) )/np.linalg.norm(fuclrum)




origin = [ts[0][0][0],ts[0][1][0],ts[0][2][0]]
print("ORIGIN", origin)



fig = plt.figure(figsize=(7,7))
ax  = fig.add_subplot(111, projection='3d')
for ext in exits:
    ax.plot(ext[0], ext[1],ext[2], '^', markersize=10, color='purple', alpha=1)

colors = ["b", "g", "r", "c", "m", "y", "k", "w"]

for i, tunnel in enumerate(ts):
    ax.plot(tunnel[0,:], tunnel[1,:], tunnel[2,:], '-', markersize=4, color=colors[i % len(colors)], alpha=0.7, label=f'Tunnel {i+1}')

# ax.plot(t1[0,:], t1[1,:], t1[2,:], 'o', markersize=8, color='green', alpha=0.2)
ax.plot([ exits[0][0],exits[1][0] ],[ exits[0][1],exits[1][1] ], [ exits[0][2],exits[1][2] ])
# universal origin 
ax.plot(origin[0],origin[1],origin[2], 'o', markersize=7, color='green', alpha=1)


final = []
for tnum,t in enumerate(ts):
    print(f"\nTunnel {tnum+1} is of shape {t.shape}")
    mx = 0
    for i in range(t.shape[1]):
        point = t[:,i]
        d = dist(point)
        mx = d if d > mx else mx;
    print(f"\nMaximum distance Tunnel{tnum+1}:" ,mx)
    print('\n')
    if mx < 25:
        final.append(f"Tunnel{tnum+1}")

print("Would keep ", final)
plt.legend()
plt.show()