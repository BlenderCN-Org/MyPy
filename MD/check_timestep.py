#!/usr/bin/python
"""Using the ASEInterface to run MD with different 
timestep computing the energy drift.

   The only thing to set is the timesteps array and the 
   method you want to test."""

import os
from MyPy.ASE.ASEInterface import ASEInterface
from MyPy.PRINT.advancement_bar import percentage_with_fan

def main():
    timesteps = [0.5, 1.0, 1.5, 2.0, 2.2, 2.4, 2.6 ]   # Time steps to tests in fs
    method = 'dftb_std'    # Method to test
    virtual_time = 10000          # Virtual time in fs
    nprint_time = 2             # Interval between two printed output in fs
    drift_on_time = 10000             # Time on which the drift is computed in fs
    
    
    ase = ASEInterface(job='MD')
    ase.setMachine('workstation')
    ase.setGeneralParams(method,'c6.xyz')

    counter = 0
    m = []
    maxtime = []
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
        reg = linearRegression('check_ts_'+str(timestep)+'.out',drift_on)
        m.append(float(reg[0]))
        maxtime.append(float(reg[1]))
        counter +=1
        bar.bar(counter)

    print
    for t,d,mt in zip(timesteps, m, maxtime):
        print '%6.2f  %6.2f  %6.2f' % tuple([t,d,mt])
        

def linearRegression(filep,drift_on):
    import numpy as np

    filec = open(filep,'r').readlines()
    
    time = []
    energy = []
    for line in filec:
        if not line[0] == 'T' and line.find('nan') == -1:
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

    return w[0], time[-1]

if __name__ == '__main__':
    main()
