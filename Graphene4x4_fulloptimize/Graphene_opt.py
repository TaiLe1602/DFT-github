import numpy as np
from ase import Atoms, Atom
from ase.io.trajectory import Trajectory
from ase.calculators.emt import EMT
from ase.calculators.abinit import Abinit
from ase.io import read, write
from ase.calculators.espresso import Espresso
import ase.io as io
from ase.io import Trajectory
from ase.build import bulk
from ase.units import kJ
from ase.eos import EquationOfState
import matplotlib.pyplot as plt
from ase.lattice.hexagonal import *
from ase.constraints import UnitCellFilter
from ase.optimize import LBFGS
from sympy import diff
pseudopotentials = {'C': 'c_pbe_v1.2.uspp.F.UPF'}
gra = Atoms('C2')
gra.set_scaled_positions([(0,0,0),(1/3,2/3,0)])

input_data={'system':{  'input_dft':'PBE',
                        'occupations':'smearing',
                        'degauss':0.02,
                        'vdw_corr': 'dft-d3'}, 
            'control':{'tstress':True, 
                       'tprnfor':True,
                        'calculation':'scf', 
                         'disk_io':'None'}}

def Energy(d,t):
	calc = Espresso(pseudopotentials=pseudopotentials, kpts=(12,12,1),
                input_data = input_data,
                tstress=True, tprnfor=True,
                ecutwfc=t, ecutrho=400,
                pseudo_dir='/home/max/Documents/TestQE/all_pbe_UPF_v1.5')
	gra.set_cell([d,d,10,90,90,120],scale_atoms=True)
	gra.calc=calc
	return gra.get_total_energy()

def fitfunction(x,y):
	fit = np.polyfit(x,y,3)
	f=np.poly1d(fit)
	x=np.linspace(2.4,2.6)
	return print(x[np.argmin(f(x))])
def OPT(t):
	Earray=[]
	darray = np.linspace(2.40,2.60,11)
	for i in range(0,11):
	    Earray.append(Energy(2.40+i*0.02,t))
	print(darray)
	print(Earray)
	return fitfunction(darray,Earray)

Energy(2.45,40)
OPT(40)



#traj = Trajectory('graphene.traj','w')
#for x in np.linspace(2.40,2.60,11):
#	gra.set_cell([x,x,10,90,90,120],scale_atoms=True)	
#	gra.set_calculator(calc)
#	gra.get_total_energy()
#	traj.write(gra)
#configs = read('graphene.traj@0:11')
#bl = np.linspace(2.40,2.60,11)
##energies = [gra.get_total_energy()*0.073498688455102 for gra in configs]
#print(bl)
#print(energies)
#plt.plot(bl,energies)
#plt.axis([2.4,2.6,-22.985,-22.950])
#plt.savefig("test_rasterization.pdf", dpi=150)
