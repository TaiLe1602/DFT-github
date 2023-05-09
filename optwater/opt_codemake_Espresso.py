from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.calculators.espresso import Espresso
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
pseudopotentials = {'H': 'h_pbe_v1.4.uspp.F.UPF','O': 'o_pbe_v1.2.uspp.F.UPF'}
calc = Espresso(pseudopotentials=pseudopotentials,
                tstress=True, tprnfor=True,
                ecutwfc=40, ecutrho=400, kpts=(1,1,1),
                pseudo_dir='/home/max/Documents/TestQE/all_pbe_UPF_v1.5')
water = Atoms('H2O',
              positions=[(d, 0, 0),
                         (d * np.cos(t), d * np.sin(t), 0),
                         (0, 0, 0)], cell=[8.0,8.1,8.2])

water.calc = calc
print(water.get_total_energy())
waterEs = Trajectory('waterEshandmk','w')
def Energy(d,t):
	calc = Espresso(pseudopotentials=pseudopotentials,
                tstress=True, tprnfor=True,
                ecutwfc=40, ecutrho=400, kpts=(1,1,1),
                pseudo_dir='/home/max/Documents/TestQE/all_pbe_UPF_v1.5')
	water = Atoms('H2O',
              positions=[(d, 0, 0),
                         (d * np.cos(t), d * np.sin(t), 0),
                         (0, 0, 0)], cell=[8.0,8.1,8.2])
	
	water.calc = calc
	return water.get_total_energy()
for i in range(0,21):
	e.append([]) #nhap gia tri append
	for j in range (0,21):
		e[i].append(Energy(0.95+i*0.005,np.pi*(95+0.5*j)/180))
angle = np.linspace(95,105,21)
print(angle)
bondlenght = np.linspace(0.95,1.05,21)
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





