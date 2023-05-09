from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT
import numpy as np
from ase.io import Trajectory
from ase.io import read
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
plt.style.use('seaborn-poster')
d = 0.9575
t = np.pi / 180 * 104.51
e = []
water = Atoms('HHO',[(d, 0, 0),(d * np.cos(t), d * np.sin(t), 0),(0, 0, 0)],calculator=EMT())

print(water.get_total_energy())
water.set_positions([(d, 0, 0),(d * np.cos(t), d * np.sin(t), 0),(0, 0, 0)])
print(water.get_positions())
waterEMT = Trajectory('waterEMThandmk','w')
def Energy(d,t):
	water = Atoms('HHO',[(d, 0, 0),(d * np.cos(t), d * np.sin(t), 0),(0, 0, 0)],calculator=EMT())
	return water.get_total_energy()
for i in range(0,101):
	e.append([])
	for j in range (0,101):
		e[i].append(Energy(0.90+i*0.002,np.pi*(90+0.2*j)/180))
angle = np.linspace(90,110,101)
print(angle)
bondlenght = np.linspace(0.90,1.10,101)
X, Y = np.meshgrid(angle, bondlenght)
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(X, Y, e, 50, cmap='binary')
ax.set_xlabel('angle')
ax.set_ylabel('bondlength')
ax.set_zlabel('Total Engergy');
ax.view_init(60, 35)
fig
plt.show()





