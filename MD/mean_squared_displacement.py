#!/usr/bin/python

from xyz_trajsplit import mdxyz
import numpy as np
import numpy.linalg as LA
import sys

# Diffusion coefficient of water: 2.299E-9 m^2/s


filep = sys.argv[1]
timestep = float(sys.argv[2])
linear_regression_points = 500

y = []
xyz = mdxyz(filep)
xyz.read_xyz()

all_atoms='1-'+str(xyz.natom)

f0 = np.array(xyz.atom_xyz(all_atoms,frame=1))
for frame in xrange(1,xyz.frames):
    f1 = np.array(xyz.atom_xyz(all_atoms,frame=frame))
    
    mean_sqrt_displ = np.sum((f1-f0)**2)/xyz.natom
    y.append(mean_sqrt_displ)
#    print '%12.6g    %12.6g' % ((frame-1)*timestep,mean_sqrt_displ)



x = np.arange((xyz.frames-linear_regression_points)*timestep,xyz.frames*timestep,timestep)
yr = y[-1*linear_regression_points:]

A = np.vstack([x,np.ones(len(x))]).T
#slope, intercept,r_value,p_value,std_err =  LA.lstsq(A,y)
params,r_value =  LA.lstsq(A,yr)[:2]

print 'Diffusion: ',params[0]/6 * (1E-10)**2 *(1E15)

import pylab as pl

x = np.arange(timestep,xyz.frames*timestep,timestep)
line = params[0]*x+params[1]
pl.plot(x,y,'o',x,line,'r-')
pl.show()
