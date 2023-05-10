from ase.build import bulk
from ase.calculators.espresso import Espresso
from ase.constraints import UnitCellFilter
from ase.optimize import LBFGS
from ase.io import write, read

pseudopotentials = {'Na': 'na_pbe_v1.5.uspp.F.UPF',
                    'Cl': 'cl_pbe_v1.4.uspp.F.UPF'}


rocksalt = bulk('NaCl', crystalstructure='rocksalt', a=6.0)
calc = Espresso(pseudopotentials=pseudopotentials,
                tstress=True, tprnfor=True, kpts=(3, 3, 3),
                ecutwfc=40, ecutrho=200,
                pseudo_dir='/home/max/Documents/TestQE/all_pbe_UPF_v1.5')

rocksalt.set_calculator(calc)

ucf = UnitCellFilter(rocksalt)
opt = LBFGS(ucf,trajectory='opt.traj')
opt.run(fmax=0.05)
print(rocksalt.get_total_energy())

