from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT
import numpy as np
from ase.io import Trajectory
from ase.io import read, write
import ase.io as io
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
from ase.calculators.espresso import Espresso
plt.style.use('seaborn-poster')
d = 0.9575
t = np.pi / 180 * 104.51
phi = 104.51
pseudopotentials = {'H': 'h_pbe_v1.4.uspp.F.UPF','O': 'o_pbe_v1.2.uspp.F.UPF'}
calc = Espresso(pseudopotentials=pseudopotentials,
                tstress=True, tprnfor=True,
                ecutwfc=40, ecutrho=400, kpts=(1,1,1),
                pseudo_dir='/home/max/Documents/TestQE/all_pbe_UPF_v1.5')
water = Atoms('H2O',
              positions=[(d, 0, 0),
                         (d * np.cos(t), d * np.sin(t), 0),
                         (0, 0, 0)], cell=[8.0,8.1,8.2])
io.write('water_Espresso.xyz', water, format='xyz')
water.calc = calc
print(water.get_total_energy())
dyn = BFGS(water, trajectory ='water_Es.traj',restart='opt_Es.traj')
dyn.run(fmax=0.05)
atoms = read('water_Es.traj')
configs = read('water_Es.traj@0:16')
energies = [water.get_total_energy() for water in configs]
angleHOH = [water.get_angle(1,2,0) for water in configs]
bondlenght= [water.get_distance(1,2) for water in configs]
print(energies)
print(angleHOH)
print(bondlenght)
x = angleHOH
y = bondlenght
z = energies 
fig = plt.figure(figsize = (10,10))
ax = plt.axes(projection='3d')
ax.grid()

ax.scatter(x, y, z, c = 'r', s = 50)
ax.set_title('3D Scatter Plot')

# Set axes label
ax.set_xlabel('x', labelpad=20)
ax.set_ylabel('y', labelpad=20)
ax.set_zlabel('z', labelpad=20)

plt.show()



