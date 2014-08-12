import os

# os.environ['DFTB_PREFIX'] = '/home/petragli/Store/SK-parameters/3ob-1-1/'
# os.environ['DFTB_COMMAND'] = '/home/petragli/bin/dftb+'
# os.environ['dDMC_COMMAND'] = '/home/petragli/Software/dDMC2/src/dDMC2'
os.environ['DFTB_COMMAND'] = '/home/petragli/data1/Software/bin/dftb+'
os.environ['dDMC_COMMAND'] = '/home/petragli/dDMC/src/dDMC2'
os.environ['DFTB_PREFIX'] = '/home/petragli/data2/Store/SK-parameters/3ob-1-1/'


directory = os.path.basename(os.getcwd())
method = '-'.join(directory.split('-')[1:])
print 'METHOD: #'+method.strip()+'#'


# General
molecule_file_path = './struct.xyz'
#method = 'dftbdDMC'
#method = 'D3'

#MD & GEOM
traj_output = 'traj.out'


#MD
temperature = 50 #K
time_step = 0.5 #fs
print_every = 100 #step
nstep = 200000 #step
#prop_output = '/home/petragli/MD-dftbdDMC/4T-4mon/prop.out'


#GEOM
fmax = 0.001


print method

        
if method == 'dftb_std-dDMC':
    from ase.calculators.dftbddmc import DftbdDMC
elif method == 'dftb_std':
    from ase.calculators.dftb import Dftb
elif method == 'dftb_std-D3':
    os.environ['D3'] = os.path.join(os.environ['HOME'],'Software/DFTD3/dftd3')
    os.environ['D3DF'] = 'zero'
    os.environ['H4_correction']='no'
    os.environ['D3FUNC'] = 'scc-dftb+h4'
    os.environ['H4'] = os.path.join(os.environ['HOME'],'Software/H_BONDS4/h_bonds4.DFTB')
    from ase.calculators.dftbd3h4 import DftbD3H4
elif method == 'dftb_std-D3H4':
    os.environ['D3'] = os.path.join(os.environ['HOME'],'Software/DFTD3/dftd3')
    os.environ['D3DF'] = 'zero'
    os.environ['D3FUNC'] = 'scc-dftb+h4'
    os.environ['H4'] = os.path.join(os.environ['HOME'],'Software/H_BONDS4/h_bonds4.DFTB')
    from ase.calculators.dftbd3h4 import DftbD3H4
else:
    print 'Method not existing: STOP!'
    exit(1)

from ase import Atoms
#from ase.optimize import QuasiNewton
from ase.data.molecules import molecule
from ase.io import write, read
import numpy as np


mol = read(molecule_file_path)

atoms = mol.get_chemical_symbols()

    
hubbard_derivs = {
    'c': -0.1492,
    'o': -0.1575,
    'n': -0.1535,
    'h': -0.1857,
    's': -0.11,     # From JCTC 2013 Q. Cui
    'p': -0.14      #
}
max_ang_mom = {
    'c': '"p"',
    'o': '"p"',
    'n': '"p"',
    'h': '"s"',
    's': '"d"',
    'p': '"d"'
}
spin_constants = {
    'c': '{-0.028 -0.024 -0.024 -0.022 }',
    'o': '{ -0.032 -0.028 -0.028 -0.027 }',
    'n': '{ -0.030 -0.026 -0.026 -0.025 }',
    'h': '{ -0.064 }'
}



charge = 0
spin_multiplicity = 1

dftb_parameters={}

dftb_parameters['Hamiltonian_MaxAngularMomentum_'] = ''
dftb_parameters['Hamiltonian_HubbardDerivs_'] = ''
dftb_parameters['Hamiltonian_SCC'] = 'Yes'
dftb_parameters['Hamiltonian_ThirdOrderFull'] = 'Yes'
dftb_parameters['Hamiltonian_DampXH'] = 'Yes'
dftb_parameters['Hamiltonian_DampXHExponent'] = 4.00
dftb_parameters['Hamiltonian_Charge'] = float(charge)
if spin_multiplicity > 1:
    dftb_parameters['Hamiltonian_SpinPolarisation_'] = 'Colinear'
    dftb_parameters['Hamiltonian_SpinPolarisation_UnpairedElectrons'] = int(spin_multiplicity)-1
    dftb_parameters['Hamiltonian_SpinConstants_'] = ''

for at in atoms:
#    if at.lower() not in max_ang_mom or at.lower() not in hubbard_derivs or at.lower() not in spin_constants:
    if at.lower() not in max_ang_mom or at.lower() not in hubbard_derivs:
        raise Exception('Parameters for '+str(at)+' are missing in calculator.py')
    dftb_parameters['Hamiltonian_MaxAngularMomentum_'+str(at)] = max_ang_mom[at.lower()]
    dftb_parameters['Hamiltonian_HubbardDerivs_'+str(at)] = hubbard_derivs[at.lower()]
    if spin_multiplicity > 1 :
        dftb_parameters['Hamiltonian_SpinConstants_'+str(at)] = spin_constants[at.lower()]

dftb_parameters['Hamiltonian_SCC'] = 'Yes'

ddmc_parameters = {}
ddmc_parameters['dftype'] = '3'
ddmc_parameters['param_a'] = 1.86675017658353
ddmc_parameters['param_b'] = 1.02244211019925
ddmc_parameters['param_s'] = 23.0
ddmc_parameters['readparamsflag'] = 'UP'


if method == 'dftb_std-dDMC':
    calculator = DftbdDMC(atoms=mol,label='cazzo',dftbdict=dftb_parameters, ddmcdict=ddmc_parameters)
elif method == 'dftb_std':
    calculator = Dftb(atoms=mol,label='asd')
elif method == 'dftb_std-D3':
    calculator = DftbD3H4(atoms=mol,label='cazzo',dftbdict=dftb_parameters )
elif method == 'dftb_std-D3H4':
    calculator = DftbD3H4(atoms=mol,label='cazzo',dftbdict=dftb_parameters )
else:
    print 'Method not existing: STOP!'
    exit(2)

calculator.parameters.update(dftb_parameters)

energy = mol.get_potential_energy()
anal = mol.get_forces()

#print energy
#print anal
#exit()
# # Compute the numerical Derivative!
# def numeric_force(atoms, a, i, d=0.001):
#     """Evaluate force along i'th axis on a'th atom using finite difference.

#     This will trigger two calls to get_potential_energy(), with atom a moved
#     plus/minus d in the i'th axial direction, respectively.
#     """
#     p0 = atoms.positions[a, i]
#     atoms.positions[a, i] += d
#     eplus = atoms.get_potential_energy()
#     print 
#     atoms.positions[a, i] -= 2 * d
#     eminus = atoms.get_potential_energy()
#     atoms.positions[a, i] = p0
#     return (eminus - eplus) / (2 * d)

# def numeric_forces(atoms, indices=None, axes=(0, 1, 2), d=0.001,
#                    parallel=None, name=None):
#     """Evaluate finite-difference forces on several atoms.

#     Returns an array of forces for each specified atomic index and
#     each specified axis, calculated using finite difference on each
#     atom and direction separately.  Array has same shape as if
#     returned from atoms.get_forces(); uncalculated elements are zero.

#     Calculates all forces by default."""
    
#     import numpy as np
#     from ase.parallel import world, rank, distribute_cpus
#     from ase.utils import opencew


#     if indices is None:
#         indices = range(len(atoms))
#     F_ai = np.zeros_like(atoms.positions)
#     n = len(indices) * len(axes)
#     if parallel is None:
#         atom_tasks = [atoms] * n
#         master = True
#         calc_comm = world
#     else:
#         calc_comm, tasks_comm, tasks_rank = distribute_cpus(parallel, world)
#         master = calc_comm.rank == 0
#         calculator = atoms.get_calculator()
#         calculator.set(communicator=calc_comm)
#         atom_tasks = [None] * n
#         for i in range(n):
#             if ((i - tasks_rank) % tasks_comm.size) == 0:
#                 atom_tasks[i] = atoms
#     for ia, a in enumerate(indices):
#         for ii, i in enumerate(axes):
#             atoms = atom_tasks[ia * len(axes) + ii]
#             if atoms is not None:
#                 done = 0
#                 if name:
#                     fname = '%s.%d%s.pckl' % (name, a, 'xyz'[i])
#                     fd = opencew(fname, calc_comm)
#                     if fd is None:
#                         if master:
#                             try:
#                                 F_ai[a, i] = pickle.load(open(fname))
#                                 print '# atom', a, 'xyz'[i], 'done'
#                                 done = 1
#                             except EOFError:
#                                 pass
#                         done = calc_comm.sum(done)
#                 if not done:
#                     print '# rank', rank, 'calculating atom', a, 'xyz'[i]
#                     force = numeric_force(atoms, a, i, d)
#                     if master:
#                         F_ai[a, i] = force
#                         if name:
#                             fd = open('%s.%d%s.pckl' % (name, a, 'xyz'[i]),
#                                       'w')
#                             pickle.dump(force, fd)
#                             fd.close()
#     if parallel is not None:
#         world.sum(F_ai)
#     return F_ai

# numer =  numeric_forces(mol, indices=None, axes=(0, 1, 2), d=0.001, parallel=None, name=None)


# print energy
# print numer
# print anal
# print
# print numer/anal

# Geometry optimization
#from ase.optimize import BFGS
from ase.optimize.bfgslinesearch import BFGSLineSearch
from ase.io.trajectory import PickleTrajectory

#dyn = BFGS(mol, trajectory=traj_output)
dyn = BFGSLineSearch(mol, trajectory=traj_output)
# traj = PickleTrajectory(traj_output, 'a', mol)
# dyn.attach(traj.write, interval = 1)
dyn.run(fmax=fmax)

# # MD NVE
# from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
# from ase.md.langevin import Langevin
# from ase import units
# from ase.md import MDLogger
# from ase.io.trajectory import PickleTrajectory

# useAsap = False

# MaxwellBoltzmannDistribution(mol,temperature*units.kB)

# dyn = Langevin(mol, time_step*units.fs, temperature*units.kB, 0.01 )

# #Function to print the potential, kinetic and total energy
# def printenergy(a=mol):
#     epot = a.get_potential_energy() / len(a)
#     ekin = a.get_kinetic_energy() / len(a)
#     print ("Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  Etot = %.3feV" %
#            (epot, ekin, ekin/(1.5*units.kB), epot+ekin))


# traj = PickleTrajectory(traj_output,"a",mol)
# dyn.attach(MDLogger(dyn,mol,prop_output),interval=print_every)
# dyn.attach(traj.write, interval = print_every)


# dyn.run(nstep)
# traj.close()
