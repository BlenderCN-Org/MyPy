import os
import sys
import check
from errors_handling import *


class ASEInterface():
    """This class allow to run computation using the ASE package in an easy way."""

    def __init__(self, job = 'singlePoint'):
        """Define the kind of job you would like to do.

           Possible options:
             - singlePoint
             - geomOpt
             - MD
             - compareForces
        """

        self.definedParams = {}
        self.definedParams['job'] = job


    def setMachine(self, machine):
        """Set the necessary enviromental variables depending on the machine.

        Set which computer you are using:
         - machine <string> {workstation|lcmdlc1|lcmdlc2}
        """

        machines = ['workstation','lcmdlc1','lcmdlc2']
        if not machine in machines:
            raise InputError(machine, 'The chosen machine does not exist')
        self.definedParams['machine'] = machine


    def setGeneralParams(self, method, initial_structure, charge=0, spinm=1, workdir=os.getcwd()):
        """Set the always necessary parameters.

           Set in order the following parameters:
            - method   <string>   {dftb_std|dftb_std-D3|dftb_std-dDMC|dftb_std-D3H4}
            - xyz file <string>   Path of the xyz file
            - charge   <integer>  Charge of the molecules
            - spinm    <integer>  Spin Multiplicity
            - workdir  <string>   Path of where to save everithing (default = .)
        """

        methods = ['dftb_std', 'dftb_std-D3', 'dftb_std-dDMC', 'dftb_std-D3H4', 'dftb_mio11']

        # Check if file exists and if method is a valid method
        if  method not in methods:
            raise InputError(method,'The chosen method does not exist!')
        check.iffile(initial_structure)
        check.ifdir(workdir)
        check.ifinteger(charge, 'The charge has to be an integer')
        check.ifinteger(spinm, 'The spin multiplicity has to be an integer')

        self.definedParams['method'] = str(method)
        self.definedParams['initial_structure'] = str(initial_structure)
        self.definedParams['workdir'] = str(workdir)
        self.definedParams['charge'] = int(charge)
        self.definedParams['spinm'] = int(spinm)


    def setMDParams(self,thermostat,time_step,nstep,output_prefix,
                    nprint,init_temperature,temperature):
        """Set the parameters needed to run an MD computation.

        Set in order the following parameters:
         - thermostat       <string>   {Langevin|None} [If None an NVE simulation is performed]
         - time_step        <float>
         - nstep            <integer>                  [Number of step]
         - output_prefix    <string>                   [Prefix for the output files]
         - nprint           <integer>                  [Print properties and xyz every nprint steps]
         - init_temperature <float>                    [Initial temperature]
         - temperature      <float>                    [Temperature - placeholder if thermostat is None]
        """

        # Check that everithing is right
        thermostats = ['Langevin','None']
        if not thermostat in thermostats:
            raise InputError(thermostat,'The chosen thermostat does not exist!')
        check.iffloat(time_step, 'The time step has to be a float!')
        check.iffloat(temperature, 'The temperature has to be a float!')
        check.iffloat(init_temperature, 'The init_temperature has to be a float!')
        check.ifinteger(nstep, 'The nstep has to be an integer!')
        check.ifstring(output_prefix, 'The output_prefix has to be a string!')
        check.ifinteger(nprint, 'The nprint has to be an integer!')

        self.definedParams['thermostat']        = str(thermostat)
        self.definedParams['time_step']         = float(time_step)
        self.definedParams['temperature']       = float(temperature)
        self.definedParams['init_temperature']  = float(init_temperature)
        self.definedParams['nstep']             = int(nstep)
        self.definedParams['output_prefix']     = str(output_prefix)
        self.definedParams['nprint']            = int(nprint)


    def setGeomOptParams(self, output_prefix, fmax=0.001, driver='BFGSLineSearch'):
    # def setGeomOptParams(self, output_prefix, fmax=0.001, driver='BFGS'):
        """Set the parameters needed to run a Geometry optimization.

        Set in order the following parameters:
         - fmax          <float>   [Threshold: maximum force intensity allowded]
         - output_prefix <string>  [See MD]
        """

        check.iffloat(fmax, 'The fmax has to be a float!')
        check.ifstring(output_prefix, 'The output_prefix has to be a string!')

        self.definedParams['fmax'] = float(fmax)
        self.definedParams['output_prefix'] = str(output_prefix)
        self.definedParams['driver'] = str(driver)

    def __xyzReader(self):
        """Read columns other then the 4th. If those columns are present they have a meaning:
            P -> position of the respective atom is constrained (if the atom is not constrainde you MUST specify 'f' (free))
            float -> is the charge of the atom. All the charge MUST be specified!"""
        import numpy as np
        rewritexyz = False
        xyz = np.loadtxt(self.definedParams['initial_structure'],skiprows=2,usecols=[1,2,3])
        atoms = np.loadtxt(self.definedParams['initial_structure'], skiprows=2,usecols=[0],dtype=str)
        lastcolumns = []
        try:
            i = 4
            while True:
                lastcolumns.append(np.loadtxt(self.definedParams['initial_structure'], skiprows=2,usecols=[i], dtype=str))
                i +=1
                rewritexyz = True
        except:
            pass

        self.__constraint = (np.zeros_like(atoms,  dtype=int) == 1).tolist()

        for col in lastcolumns:
            for tp in [float, int,  str]:
                if isinstance(col[0], tp):
                    col = np.array(col, dtype=tp)
                    break
            if np.issubdtype(col.dtype, int):
                pass
            elif np.issubdtype(col.dtype,  str):
                self.__constraint = (np.array(map(lambda x: x.upper(), col)) == 'P').tolist()
            elif np.issubdtype(col.dtype, float):
                self.__charges = col

        if rewritexyz:
            self.definedParams['initial_structure'] = 'generated_from_ase-'+self.definedParams['initial_structure']
            f = open(self.definedParams['initial_structure'], 'w')
            f.write('%i\n\n' % len(atoms))
            for atom, coord in zip(atoms, xyz):
                f.write('%s   %12.6f   %12.6f   %12.6f\n' % (atom,  coord[0], coord[1], coord[2]))
            f.close


    def __applyConstraints(self):
        from ase.constraints import FixAtoms
<<<<<<< HEAD
        if any(self.__constraint):
            print 'Position Constrained: '+' '.join(map(lambda x: str(x), self.__constraint))
        self.__fixedAtoms = FixAtoms(mask = self.__constraint)
=======
        if not any(self.__constraint):
            print 'Position Constrained: '+' '.join(map(lambda x: str(x), self.__constraint))
            self.__fixedAtoms = FixAtoms(mask = self.__constraint)
>>>>>>> 024e04470e2f616ac4b72704507ff58d3bc902b6

    def getParams(self):
        """Print the value of all the parameters

        If this output is written to a file then can be reread using
        the readParams method.
        """

        msg = '#%20s --> %s\n' % tuple(['Params','Value'])
        for k,v in self.definedParams.items():
            msg += '%20s --> %s\n' % tuple([k,v])

        return msg

    def readParams(self, filepath):
        """Read the parameters saved with getParams.

        After readed all the parameters you can directly start the
        job with a call to job_start. Offcourse you need to specify
        all the needed parameters.

        WARN! Using with caution: there are no control on the variables.
        """

        check.iffile(filepath)

        with open(filepath,'r') as f:
            for line in f.readlines():
                if line[0] != '#':
                    k,v = line.split('-->')
                    self.definedParams[k.strip()] = v.strip()


    def getEnvironment(self):
        """Print the contents of all the enviromental variables used by ASE in this context."""

        self.__setEnvironment()

        env = ['DFTB_COMMAND','dDMC_COMMAND','DFTB_PREFIX','D3','D3DF','H4_correction','D3FUNC','H4']

        msg = '%20s --> %20s\n' % tuple(['#Variable','Value'])
        for v in env:
            try:
                msg += '%20s --> %20s\n' % tuple([v,os.environ[v]])
            except:
                msg += '%20s --> %20s\n' % tuple([v,'IS MISSING!!'])

        return msg


    def __setEnvironment(self):

        machine = self.definedParams['machine']
        method = self.definedParams['method']


        if machine == 'workstation':
            os.environ['DFTB_COMMAND'] = '/home/petragli/data1/Software/bin/dftb+'
            os.environ['dDMC_COMMAND'] = '/home/petragli/dDMC/src/dDMC2'
            os.environ['DFTB_PREFIX'] = '/home/petragli/data2/Store/SK-parameters/3ob-1-1/'
            if method == 'dftb_mio11':
                os.environ['DFTB_PREFIX'] = '/home/petragli/data2/Store/SK-parameters/mio-1-1-NHorg/'
        elif machine == 'lcmdlc1':
            os.environ['DFTB_PREFIX'] = '/lcmd-data/petragli/Store/SK-parameters/3ob-1-1/'
            os.environ['DFTB_COMMAND'] = '/lcmd-data/petragli/Software/bin/dftb+'
            os.environ['dDMC_COMMAND'] = '/lcmd-data/petragli/Software/dDMC/src/dDMC2'
            if method == 'dftb_mio11':
                raise ImplementationError(method, 'Not implemented for this machine')
        elif machine == 'lcmdlc2':
            os.environ['DFTB_PREFIX'] = '/home/petragli/Store/SK-parameters/3ob-1-1/'
            os.environ['DFTB_COMMAND'] = '/home/petragli/bin/dftb+'
            os.environ['dDMC_COMMAND'] = '/home/petragli/Software/dDMC/src/dDMC2'
            if method == 'dftb_mio11':
                raise ImplementationError(method, 'Not implemented for this machine')
        else:
            raise ImplementationError(machine, 'Enviromental variable for this machine are not implemented')


        if method == 'dftb_std-dDMC':
            pass
        elif method == 'dftb_std':
            pass
        elif method == 'dftb_mio11':
            pass
        elif method == 'dftb_std-D3':
            os.environ['D3'] = os.path.join(os.environ['HOME'],'Software/DFTD3/dftd3')
            os.environ['D3DF'] = 'zero'
            os.environ['H4_correction']='no'
            os.environ['D3FUNC'] = 'scc-dftb+h4'
            os.environ['H4'] = os.path.join(os.environ['HOME'],'Software/H_BONDS4/h_bonds4.DFTB')
        elif method == 'dftb_std-D3H4':
            os.environ['D3'] = os.path.join(os.environ['HOME'],'Software/DFTD3/dftd3')
            os.environ['D3DF'] = 'zero'
            os.environ['D3FUNC'] = 'scc-dftb+h4'
            os.environ['H4'] = os.path.join(os.environ['HOME'],'Software/H_BONDS4/h_bonds4.DFTB')
        else:
            raise ImplementationError(method,'Enviromental variables for this method are not implemented')


    def __ASEInit(self):
        # Commented because signed from the editor
        #from ase import Atoms
        #from ase.data.molecules import molecule
        from ase.io import read
        #import numpy as np

        self.mol = read(self.definedParams['initial_structure'])
        self.mol.set_constraint(self.__fixedAtoms)


    def __setDFTBParams(self):
        """Set all the parameters necessary to run teh DFTB computation."""

        # 3OB parameters
        hubbard_derivs = {
            'c': -0.1492,
            'o': -0.1575,
            'n': -0.1535,
            'h': -0.1857,
            's': -0.11,     # From JCTC 2013 Q. Cui
            'p': -0.14      #
        }
        if self.definedParams['method'] == 'dftb_mio11':
            # Mio11 parameters
            hubbard_derivs = {
                'c': -0.1492,
                'o': -0.1575,
                'n': -0.1535,
                'h': -0.07,
                's': -0.069,     # From JCTC 2013 Q. Cui
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

        atoms = self.mol.get_chemical_symbols()

        charge = self.definedParams['charge']
        spin_multiplicity = self.definedParams['spinm']

        self.__dftb_parameters={}

        self.__dftb_parameters['Hamiltonian_MaxAngularMomentum_'] = ''
        self.__dftb_parameters['Hamiltonian_HubbardDerivs_'] = ''
        self.__dftb_parameters['Hamiltonian_SCC'] = 'Yes'
        self.__dftb_parameters['Hamiltonian_ThirdOrderFull'] = 'Yes'
        self.__dftb_parameters['Hamiltonian_DampXH'] = 'Yes'
        self.__dftb_parameters['Hamiltonian_DampXHExponent'] = 4.00 # 3OB
        if self.definedParams['method'] == 'dftb_mio11':
            self.__dftb_parameters['Hamiltonian_DampXHExponent'] = 4.95 # MIO11
        self.__dftb_parameters['Hamiltonian_Charge'] = float(charge)
        if spin_multiplicity > 1:
            self.__dftb_parameters['Hamiltonian_SpinPolarisation_'] = 'Colinear'
            self.__dftb_parameters['Hamiltonian_SpinPolarisation_UnpairedElectrons'] = int(spin_multiplicity)-1
            self.__dftb_parameters['Hamiltonian_SpinConstants_'] = ''

        for at in atoms:
            if at.lower() not in max_ang_mom or at.lower() not in hubbard_derivs:
                raise Exception('DFTB parameters for '+str(at)+' are missing')
            self.__dftb_parameters['Hamiltonian_MaxAngularMomentum_'+str(at)] = max_ang_mom[at.lower()]
            self.__dftb_parameters['Hamiltonian_HubbardDerivs_'+str(at)] = hubbard_derivs[at.lower()]
            if spin_multiplicity > 1 :
                if at.lower() not in spin_constants.keys():
                    raise Exception('DFTB parameters for '+str(at)+' are missing')
                self.__dftb_parameters['Hamiltonian_SpinConstants_'+str(at)] = spin_constants[at.lower()]

        self.__dftb_parameters['Hamiltonian_SCC'] = 'Yes'


    def __setdDMCParams(self):
        """Set all the parameters to run dDMC computation."""

        self.__ddmc_parameters = {}
        self.__ddmc_parameters['dftype'] = '3'
        self.__ddmc_parameters['param_a'] = 1.86675017658353
        self.__ddmc_parameters['param_b'] = 1.02244211019925
        self.__ddmc_parameters['param_s'] = 23.0
        self.__ddmc_parameters['readparamsflag'] = 'UP'


    def __setASECalculator(self):
        """Initialize the calculator and set the parameters"""

        method = self.definedParams['method']

        if method == 'dftb_std-dDMC':
            self.__setdDMCParams()
            from ase.calculators.dftbddmc import DftbdDMC
            self.calculator = DftbdDMC(atoms=self.mol,label='cazzo',dftbdict=self.__dftb_parameters, ddmcdict=self.__ddmc_parameters)
        elif method == 'dftb_std' or method == 'dftb_mio11':
            from ase.calculators.dftb import Dftb
            self.calculator = Dftb(atoms=self.mol,label='asd')
        elif method == 'dftb_std-D3':
            from ase.calculators.dftbd3h4 import DftbD3H4
            self.calculator = DftbD3H4(atoms=self.mol,label='cazzo',dftbdict=self.__dftb_parameters )
        elif method == 'dftb_std-D3H4':
            from ase.calculators.dftbd3h4 import DftbD3H4
            self.calculator = DftbD3H4(atoms=self.mol,label='cazzo',dftbdict=self.__dftb_parameters )
        else:
            raise ImplementationError(method,'Calculator for this method is not implemented')

        self.calculator.parameters.update(self.__dftb_parameters)


    # Compute the numerical Derivative!
    def __numeric_force(self,atoms, a, i, d=0.001):
        """Evaluate force along i'th axis on a'th atom using finite difference.

        This will trigger two calls to get_potential_energy(), with atom a moved
        plus/minus d in the i'th axial direction, respectively.
        """
        p0 = atoms.positions[a, i]
        atoms.positions[a, i] += d
        eplus = atoms.get_potential_energy()
        atoms.positions[a, i] -= 2 * d
        eminus = atoms.get_potential_energy()
        atoms.positions[a, i] = p0
        return (eminus - eplus) / (2 * d)

    def numeric_forces(self,atoms, axes=(0, 1, 2), d=0.001,
                       parallel=None, name=None):
        """Evaluate finite-difference forces on several atoms.

        Returns an array of forces for each specified atomic index and
        each specified axis, calculated using finite difference on each
        atom and direction separately.  Array has same shape as if
        returned from atoms.get_forces(); uncalculated elements are zero.

        Calculates all forces by default."""

        import numpy as np
        from ase.parallel import world, distribute_cpus
        from ase.utils import opencew

        indices = range(len(atoms))
        F_ai = np.zeros_like(atoms.positions)
        n = len(indices) * len(axes)
        total_calculation = len(indices)*len(axes)

        if parallel is None:
            atom_tasks = [atoms] * n
            master = True
            calc_comm = world
        else:
            calc_comm, tasks_comm, tasks_rank = distribute_cpus(parallel, world)
            master = calc_comm.rank == 0
            calculator = atoms.get_calculator()
            calculator.set(communicator=calc_comm)
            atom_tasks = [None] * n
            for i in range(n):
                if ((i - tasks_rank) % tasks_comm.size) == 0:
                    atom_tasks[i] = atoms
        counter = 0
        for ia, a in enumerate(indices):
            for ii, i in enumerate(axes):
                atoms = atom_tasks[ia * len(axes) + ii]
                if atoms is not None:
                    done = 0
                    if name:
                        fname = '%s.%d%s.pckl' % (name, a, 'xyz'[i])
                        fd = opencew(fname, calc_comm)
                        if fd is None:
                            if master:
                                try:
                                    F_ai[a, i] = pickle.load(open(fname))
                                    print '# atom', a, 'xyz'[i], 'done'
                                    done = 1
                                except EOFError:
                                    pass
                            done = calc_comm.sum(done)
                    if not done:
                        # print '# rank', rank, 'calculating atom', a, 'xyz'[i]
                        force = self.__numeric_force(atoms, a, i, d)
                        if master:
                            F_ai[a, i] = force
                            if name:
                                fd = open('%s.%d%s.pckl' % (name, a, 'xyz'[i]),
                                          'w')
                                pickle.dump(force, fd)
                                fd.close()
                counter += 1
                fan = ['-', '\\', '|', '/']
                sys.stderr.write("\r[%-100s]%3.1f%% %1s" % tuple([int(float(counter)/float(total_calculation)*100)*'=',float(counter)/float(total_calculation)*100,fan[counter%4]]))
                sys.stderr.flush()
        sys.stderr.write('\n')
        if parallel is not None:
            world.sum(F_ai)
        return F_ai


    # Define the job function:
    def singlePoint(self):
        """Return the single point energy."""
        energy = self.mol.get_potential_energy()
        return energy


    def compareForces(self):
        """Compare forces.

        Return the forces computed numerically
        and analitically with the same method."""
        energy = self.mol.get_potential_energy()
        anal = self.mol.get_forces()
        numer =  self.numeric_forces(self.mol, axes=(0, 1, 2), d=0.001,
                       parallel=None, name=None)

        return anal, numer


    def geomOpt(self):
        """Geometry Optimization"""
        from ase.optimize import BFGS
        from ase.optimize.bfgslinesearch import BFGSLineSearch
        # from ase.optimize.oldqn import GoodOldQuasiNewton
        #from ase.io.trajectory import PickleTrajectory
        drivers = {
            'BFGSLineSearch' : BFGSLineSearch,
            'BFGS'           : BFGS
        }

        # dyn = GoodOldQuasiNewton(self.mol,
        #                 trajectory=os.path.join(self.definedParams['workdir'],
        #                 self.definedParams['output_prefix'])+'.traj')

        dyn = drivers[self.definedParams['driver']](self.mol,
                        trajectory=os.path.join(self.definedParams['workdir'],
                        self.definedParams['output_prefix'])+'.traj')
        dyn.run(fmax=self.definedParams['fmax'])
        return 'done'


    def MD(self):
        """Molecular Dynamic"""
        from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
        from ase import units
        from ase.md import MDLogger
        from ase.io.trajectory import PickleTrajectory
        from ase.md.langevin import Langevin
        from ase.md.verlet import VelocityVerlet

        dyndrivers = {
            'Langevin': Langevin,
            'None': VelocityVerlet,
        }

        useAsap = False

        mol = self.mol
        temperature = self.definedParams['temperature']
        init_temperature = self.definedParams['init_temperature']
        time_step = self.definedParams['time_step']
        nstep = self.definedParams['nstep']
        nprint = self.definedParams['nprint']
        thermostat = self.definedParams['thermostat']
        prop_file = os.path.join(self.definedParams['workdir'],
                                 self.definedParams['output_prefix']+'.out')
        traj_file = os.path.join(self.definedParams['workdir'],
                                 self.definedParams['output_prefix']+'.traj')

        MaxwellBoltzmannDistribution(mol,init_temperature*units.kB)

        if thermostat == 'None':
            dyn = VelocityVerlet(mol, time_step*units.fs)
        elif thermostat == 'Langevin':
            dyn = Langevin(mol, time_step*units.fs, temperature*units.kB, 0.01 )
        else:
            raise ImplementationError(method,'Thermostat is not implemented in the MD function')

        #Function to print the potential, kinetic and total energy
        traj = PickleTrajectory(traj_file,"a",mol)
        dyn.attach(MDLogger(dyn,mol,prop_file),interval=nprint)
        dyn.attach(traj.write, interval = nprint)

        dyn.run(nstep)
        traj.close()


    # Run the computation:
    def __initCalc(self):
        """Set everething!."""
        self.__xyzReader()
        self.__applyConstraints()
        self.__setEnvironment()
        self.__ASEInit()
        self.__setDFTBParams()
        self.__setASECalculator()

    def runJob(self):
        """Interface to start the job."""

        self.__initCalc()

        jobs = {
            'singlePoint'  : self.singlePoint,
            'geomOpt'      : self.geomOpt,
            'compareForces': self.compareForces,
            'MD'           : self.MD
        }

        job = self.definedParams['job']
        return jobs[job]()



if __name__ == '__main__':

    asd = ASEInterface()

    asd.setMachine('workstation')
    asd.setGeneralParams('dftb_std-dDMC','c6.xyz')
    # asd.setGeomOptParams('testout',)
    #    asd.setMDParams('Langevin',1.0,200,'testout_MD',1,300.,300.)
    print asd.getParams()
    print asd.runJob()
