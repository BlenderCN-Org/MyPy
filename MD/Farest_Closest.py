#!/usr/bin/python

def analyzing(frame,return_dict):
    Farest = 0.0
    Closest = 1E15
    print 'Frame n.',frame
    for atom1 in xrange(1,xyz.natom):
#        print 'Atom1', atom1
        for atom2 in xrange(atom1+1,xyz.natom):
#           print frame,atom1,atom2
            a1 = np.array(xyz.atom_xyz(str(atom1),frame)[0])
#            print a1
            a2 = np.array(xyz.atom_xyz(str(atom2),frame)[0])
            dist = LA.norm(a1-a2)

            if dist > Farest:
                Farest = dist
                Farest_atoms = [atom1,atom2]
#                print 'Farest upgraded: ',Farest
            if dist < Closest:
                Closest = dist
                Closest_atoms = [atom1,atom2]
#                print 'Closest upgraded: ',Closest
    return_dict[frame] = (Farest,Farest_atoms,Closest,Closest_atoms)



if __name__ == '__main__':

    from xyz_trajsplit import mdxyz
    import numpy as np
    import numpy.linalg as LA
    import sys
    from multiprocessing import Process, Manager

    filep = sys.argv[1]


    xyz = mdxyz(filep)
    xyz.read_xyz()

    Farest = 0.0
    Closest = 1E15


    manager = Manager()
    return_dict = manager.dict()
    proc = []
    analyzed_frames = []
    
    if xyz.frames > 16:
        starting_frame = xyz.frames - 16
    else:
        starting_frame = 0

    for frame in xrange(starting_frame,xyz.frames):
        analyzed_frames.append(frame)
        p = Process(target=analyzing,args=(frame,return_dict))
        p.start()
        proc.append(p)

        for p in proc:
            p.join()


        for frame in analyzed_frames:
            if return_dict[frame][0] > Farest:
                Farest = return_dict[frame][0]
                Farest_atoms = return_dict[frame][1]
            if return_dict[frame][2] < Closest:
                Closest = return_dict[frame][2]
                Closest_atoms = return_dict[frame][3]


    
    print Farest,  Farest_atoms
    print Closest, Closest_atoms
