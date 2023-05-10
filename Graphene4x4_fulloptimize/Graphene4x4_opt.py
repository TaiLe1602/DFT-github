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
gra.set_cell([2.46,2.46,10,90,90,120],scale_atoms=True)
gra4x4=gra*(4,4,1)
input_data={'system':{  'input_dft':'PBE',
                        'occupations':'smearing',
                        'degauss':0.02,
                        'vdw_corr': 'dft-d3'}, 
            'control':{'tstress':True, 
                       'tprnfor':True,
                        'calculation':'scf', 
                         'disk_io':'None'}}

def Energy(d,t):
	calc = Espresso(pseudopotentials=pseudopotentials, kpts=(1,1,1),
                input_data = input_data,
                tstress=True, tprnfor=True,
                ecutwfc=t, ecutrho=400,
                pseudo_dir='/workspace/stuid01/Pseudopotential')
	gra.set_cell([d,d,10,90,90,120],scale_atoms=True)
	gra4x4=gra*(4,4,1) #test gra4x4 nhung phai set cell cho gra
	gra4x4.calc=calc
	return gra4x4.get_total_energy()

def fitfunction(x,y):
	fit = np.polyfit(x,y,3)
	f=np.poly1d(fit)
	x=np.linspace(2.4,2.6)
	return x[np.argmin(f(x))]
def OPT(t):
	Earray=[]
	darray = np.linspace(2.40,2.60,11)
	for i in range(0,11):
	    Earray.append(Energy(2.40+i*0.02,t))
	print(darray)
	print(Earray)
	return fitfunction(darray,Earray)
OPT(40)
bonlength =[]
for i in range(0,6):
	bonlength.append(OPT(30+i*10))
ecut=[30,40,50,60,70,80]
print(ecut)
print(bonlength)
plt.plot(ecut,bonlength)
plt.axis([30,90,2.46,2.47])
plt.savefig("optgra4x4.pdf", dpi=150)