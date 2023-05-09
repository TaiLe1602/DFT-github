from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT
import numpy as np
from ase.io import Trajectory
from ase.io import read
from ase.build import molecule
from ase.collections import g2
d = 0.9575
t = np.pi / 180 * 104.51
water = molecule('H2O', vacuum=0)
print(g2.names)
atoms =g2['H2O']
water.set_calculator(EMT())
print(water.get_total_energy())
dyn = BFGS(water, trajectory ='water_molecule_EMT.traj')
dyn.run(fmax=0.000001)
atoms = read('water_molecule_EMT.traj')
print(atoms.get_angle(0, 1, 2))
print(atoms.get_angle(2, 0, 1))
print(atoms.get_angle(1, 2, 0))

