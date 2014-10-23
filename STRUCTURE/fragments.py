#!/usr/bin/python

import os, sys
import numpy as np

labels = {'h': 1, 'c': 6, 'n': 7, 'o': 8, 'f': 9, 's': 16, 'cl': 17, 'fe': 26}
factor = 1.05
#                0.31
radius=[None,0.35,0.28,1.28,0.96,0.84,0.76,0.71,0.66,
        0.57,0.58,1.66,1.41,1.21,1.11,1.07,1.05,
        1.02,1.06,2.03,1.76,1.70,1.60,1.53,1.39,
        1.61,1.52,1.50,1.24,1.32,1.22,1.22,1.20,
        1.19,1.20,1.20,1.16,2.20,1.95,1.90,1.75,
        1.64,1.54,1.47,1.46,1.42,1.39,1.45,1.44,
        1.42,1.39,1.39,1.38,1.39,1.40,2.44,2.15,
        2.07,2.04,2.03,2.01,1.99,1.98,1.98,1.96,
        1.94,1.92,1.92,1.89,1.90,1.87,1.87,1.75,
        1.70,1.62,1.51,1.44,1.41,1.36,1.36,1.32,
        1.45,1.46,1.48,1.40,1.50,1.50,2.60,2.21,
        2.15,2.06,2.00,1.96,1.90,1.87,1.80,1.69]



def readxyz(fname):

    atomnos = []
    atomcoords = []
    linecount = 0
    for line in open(fname).readlines():
        linecount += 1
        words = line.split()
        if linecount == 2 : cstring = line
        if len(words) == 4 and linecount > 2:
            atomnos.append(labels[words[0].lower()])
            atomcoords.append(map(float, words[1:]))
    return atomnos, atomcoords,cstring

def writexyz(fname, nxyzs, cstring):

    rever_labels={}
    for it in labels.items():
        rever_labels[it[1]]=it[0]

    f = open(fname, 'w')

    f.write('%i\n  %s\n'%(len(nxyzs), cstring.rstrip()))
    for nxyz in nxyzs:
        f.write("%2s%12.6f%12.6f%12.6f\n"%(rever_labels[nxyz[0]],nxyz[1],nxyz[2],nxyz[3]))
    f.close()

    return

def intersect(list1, list2):

    result=[]
    for x in list1:

        if x in list2:
            result.append(x)

    return result

def union(list1, list2):

    result = list1
    for x in list2:

        if x not in result:
            result.append(x)

    return result

def neighbours(atomnos, atomcoords):

    voisins = {}
    for at1 in range(len(atomnos)):
        voisins[at1] = []
        xyz1 = np.array(atomcoords[at1])
        rad1 = radius[atomnos[at1]]
        voisin = {}
        for at2 in range(len(atomnos)):
            if at2 == at1 :
                continue

            xyz2 = np.array(atomcoords[at2])
            rad2 = radius[atomnos[at2]]
            d = np.sqrt(np.dot(xyz1 - xyz2, xyz1 - xyz2))
            if  d < factor * (rad1 + rad2):
                voisins[at1].append(at2)

    return voisins

def molsplit(voisins):
    cluster={}
    cl = 0
    keys = voisins.keys()

    while len(keys) > 0:
        at1 = keys[0]
        if at1 in voisins.keys():
            cluster[cl] = [at1]
            union(cluster[cl],voisins[at1])
            for at2 in cluster[cl]:
                union(cluster[cl],voisins[at2])
                del voisins[at2]
        keys = voisins.keys()
        cl += 1

    return cluster
    
    

def molecular_mass(atom_indexes,atomnos):
    mass = 0.
    for at in atom_indexes:
        mass += atomnos[at]
    return mass


def baricenter(atom_indexes,atomnos,atomcoords):
    atomcoords_a = []
    for at in atomcoords:
        atomcoords_a.append(np.array(at))
        
    b = np.zeros(3)
    for at in atom_indexes:
        b += atomcoords_a[at]*atomnos[at]
    return b/molecular_mass(atom_indexes,atomnos)


if __name__ == '__main__':

    fname = sys.argv[1]
    bname = os.path.splitext(fname)[0]

    atomnos, atomcoords, cstring = readxyz(sys.argv[1])

    voisins = neighbours(atomnos, atomcoords)
#    print voisins

    monomers = molsplit(voisins)
    
    for m in monomers.keys():
        monomer = []
        for atom in monomers[m]:
            monomer.append([atomnos[atom]]+atomcoords[atom])
        fname = '%s_mono%i.xyz'%(bname, m+1)
#        writexyz(fname, monomer, cstring)

    print "Mon 1 molecular mass: ",molecular_mass(monomers[0],atomnos)
    print "Mon 2 molecular mass: ",molecular_mass(monomers[1],atomnos)

    print "Mon 1: ", sorted(monomers[0])
    print "Mon 2: ", sorted(monomers[1])
    
    mon1_b = baricenter(monomers[0],atomnos,atomcoords)
    mon2_b = baricenter(monomers[1],atomnos,atomcoords)
    print 'Monomer Distance: %12.4f' % np.linalg.norm(mon1_b - mon2_b)
    #print type(baricenter(range(len(atomnos)),atomnos,atomcoords))
    
    
