#!/usr/bin/python
"""Using the ASEInterface to run MD with different 
timestep computing the energy drift.

   The only thing to set is the timesteps array and the 
   method you want to test."""

import os
from MyPy.ASE.ASEInterface import ASEInterface
from MyPy.PRINT.advancement_bar import percentage_with_fan

def main():
    timesteps = [0.5, 1.0, 1.5, 2.0, 2.1, 2.2, 2.3, 2.4 ]   # Time steps to tests in fs
    method = 'dftb_std-dDMC'    # Method to test
    virtual_time = 5000          # Virtual time in fs
    nprint_time = 2             # Interval between two printed output in fs
    drift_on_time = 1000             # Time on which the drift is computed in fs
    
    
    ase = ASEInterface(job='MD')
    ase.setMachine('workstation')
    ase.setGeneralParams(method,'c6.xyz')

    counter = 0
    m = []
    bar = percentage_with_fan(len(timesteps))
    for timestep in timesteps:
        outf = open('out_ts_'+str(timestep)+'.out','w')
        nstep = int(float(virtual_time)/timestep)
        nprint = int(float(nprint_time)/float(timestep))
        drift_on = int(float(drift_on_time)/float(timestep))
        ase.setMDParams('None',timestep,nstep,'check_ts_'+str(timestep),1,300,300)
        outf.write(ase.getParams())
        if not os.path.exists('check_ts_'+str(timestep)+'.out'):
            try:
                ase.runJob()
            except:
                pass
        m.append(float(linearRegression('check_ts_'+str(timestep)+'.out',drift_on)))
        counter +=1
        bar.bar(counter)

    print
    for t,d in zip(timesteps,m):
        print '%6.2f  %6.2f' % tuple([t,d])
        

def linearRegression(filep,drift_on):
    import numpy as np

    filec = open(filep,'r').readlines()
    
    time = []
    energy = []
    for line in filec:
        if not line[0] == 'T':
            t, e = line.split()[:2]
            time.append(float(t))
            energy.append(float(e))

    time = np.array(time[-drift_on:])
    A = np.array([time[-drift_on:],np.ones(len(time[-drift_on:]))])
    energy = np.array(energy[-drift_on:])

    w = np.linalg.lstsq(A.T,energy)[0]
    
    # from pylab import plot,show
    # line = w[0]*time + float(w[1])
    # plot(np.array(time),line,'r-',np.array(time),energy,'o')
    # show()

    return w[0]

if __name__ == '__main__':
    main()
