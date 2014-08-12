#!/usr/bin/python

from xyz_trajsplit import mdxyz
import sys

filep = sys.argv[1]

baric = []

xyz = mdxyz(filep)
xyz.read_xyz()

print xyz.frames

for frame in xrange(1,xyz.frames+1):
    print frame
    baric.append(xyz.atoms_baricenter('1-'+str(xyz.natom),frame=frame).tolist())
print baric
    
xyz.add_atom('He',baric)

for frame in xrange(1,xyz.frames+1):
    xyz.write_xyz('1-'+str(xyz.natom),'asd2.xyz',append=True,frame=frame)
    
