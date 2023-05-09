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
water = Atoms('H2O',
              positions=[(d, 0, 0),
                         (d * np.cos(t), d * np.sin(t), 0),
                         (0, 0, 0)])
water.set_calculator(EMT())
print(water.get_total_energy())
dyn = BFGS(water, trajectory ='water_EMT.traj',restart='opt.traj')
dyn.run(fmax=0.00000001)
configs = read('water_EMT.traj@0:16')
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



