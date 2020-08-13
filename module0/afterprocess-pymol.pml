xyz=cmd.get_coords(’Tunnel5’,1)
r=[]
cmd.iterate_state(1,’Tunnel’,’r.append(vdw)’,space=locals(),atomic=0)
python
from pymol import stored
np.savetxt(’tunnel_coordinates.txt’,xyz,fmt=’\%.2f’)
np.savetxt(’tunnel_radius.txt’,r,fmt=’\%.2f’)
python end