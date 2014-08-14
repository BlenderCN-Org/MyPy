#!/usr/bin/python

import numpy as np
import os

class meanPlane:
    
    def __init__(self,points,guess_plane):
        
        self.points = self._checkType(points)
        self.guess_plane = self._checkType(guess_plane)
        self.other_points = []
        

    def _checkType(self,l):
        testtype_nparr = type(np.array([1.,2.,3.]))
        testtype_list = type([1.,2.,3.])

        if type(l) != testtype_nparr:
            if type(l) != testtype_list:
                raise Exception("Wrong input format")
            else:
                return np.array(l)
        else:
            return l

    def _residuals(self,p,y):
        distance= np.abs(np.dot(y,p[0:3])+p[3]) / np.sqrt(np.dot(p[0:3],p[0:3]))
        return distance

    def optimize_plane(self):
        from scipy.optimize import leastsq as lstsq
        self.plsq = lstsq(self._residuals, self.guess_plane, args=(self.points),full_output=1)
        return self.plsq

    
    def point_from_plane(self,y):
        p = self.plsq[0]
        self.other_points.append(y)
#        print p
        distance= np.abs(np.dot(y,p[0:3])+p[3]) / np.sqrt(np.dot(p[0:3],p[0:3]))
        return distance

    def show_graph(self):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        a,b,c,d = self.plsq[0]
        points = np.transpose(self.points)
        
        xx,yy = np.meshgrid(np.arange(-5,5,.1),np.arange(-5,5,.1))
        if c == 0 : c = 1E-10
        zz = (-a*xx-b*yy-d)/c

        plot3d = plt.figure().gca(projection='3d')
#        plot3d.plot_surface(xx,yy,zz,color='yellow',alpha=0.2)
        points = np.transpose(self.points)
        others = np.transpose(self.other_points)
        plot3d.scatter(points[0],points[1],points[2], color='red')
        plot3d.scatter(others[0],others[1],others[2], color='green')
        plot3d.set_zlim(-10.,10.)
        plt.show()



class mdxyz():
    
    def __init__(self,filepath):
        self.filepath = filepath
        self.readed = False
        self.pi = np.arccos(-1.0)
        self.read_xyz()
        
    def read_xyz(self):
        import re
        
        self.at_xyz = []
        self.frames = 0
        self.at_type = []

        with open(self.filepath,'r') as xyzf:
            for line in xyzf:
                if re.match('^\s*\d+\s*$',line):
                    self.natom = int(line.strip())
                    self.frames += 1
                    continue
                if re.match('^\s*\w{1,2}(?:\s+[+-]?\d+\.?\d+[EDed]?[-+]?\d{0,2}){3}\s*$',line):
#                    print '#'+line.strip()+'#'
                    if len(self.at_type) < self.natom:
                        self.at_type.append(line.split()[0])
                    self.at_xyz.append(map(float,line.split()[1:4]))
                    continue
            
    def add_atom(self,at_type,at_xyz):
        import re
        match = re.match('^[a-zA-Z][a-zA-Z]$',at_type)
        
        if type(at_type) != str or type(at_xyz) != list:
            print type(at_type),type(at_xyz)
            print "The atom type has to be a string and the atom coordinates has to be a list"
            import sys
            sys.exit()

        # for n in at_xyz:
        #     if type(n) != float:
        #         print "The atom coordinates have to be float!"
        #         import sys
        #         sys.exit()
                

        if not match:
            print "The atom type has to be of two letters"
            import sys
            sys.exit()
        
        self.at_type.insert(0,at_type)

        for i in xrange(self.frames):
            try:
                xyz = at_xyz(i)
            except:
                xyz = at_xyz[-1]

            self.at_xyz.insert((i*self.natom)+i,xyz)

        self.natom += 1
        

    def frame(self,n):
        """Return a list of all the atoms and position for a specified frame"""
        fr = []
        start_index = int((n-1)*self.natom)
        end_index = int(n*self.natom)
        frame_type = self.at_type
        frame_xyz  = self.at_xyz[start_index:end_index]
        for t,xyz in zip (frame_type,frame_xyz):
             fr.append([t,xyz])
        return fr

    
    def _checklist(self,string):
        import re

        match = re.match('^(\d+)\-(\d+)$',string)
        if match:
            start_index = int(match.group(1))-1
            end_index = int(match.group(2))
        else:
            start_index = int(string)-1
            end_index = int(string)
        
        return start_index,end_index



    def atom_xyz(self,n,frame=1):
        import re
        
        atoms = []
        at_list = n.split(',')
        
        for elem in at_list:
            
            start_id, end_id = self._checklist(elem)
            start_id += self.natom*(frame-1)
            end_id += self.natom*(frame-1)
#            atoms += atoms_in_frame[start_id:end_id]
            atoms += self.at_xyz[start_id:end_id]

        return atoms


    def atom_type(self,n):
        import re
        
        atoms = []
        at_list = n.split(',')
        
        for elem in at_list:
            start_id, end_id = self._checklist(elem)
            atoms += self.at_type[start_id:end_id]

        return atoms


    def atom(self,n,frame=1):

        xyz = self.atom_xyz(n,frame=frame)
        typ = self.atom_type(n)

        return map(list,zip(typ,xyz))


    def write_xyz(self,n,filepath,frame=1,append=True):
        
        import MyPy.DATA.atomdata
        atdata = MyPy.DATA.atomdata.atomdata()

        op_mode = 'a'
        if not append or not os.path.exists(filepath): op_mode = 'w'
        atoms = self.atom(n,frame=frame)
#        print atoms
        with open(filepath,op_mode) as fxyz:
            fxyz.write(' %10i \n\n' % len(atoms))
            for at in atoms:
                try:
                    int(at[0])
                    at[0] = atdata.atnum_r[int(at[0])]
                except:
                    pass
                fxyz.write('%4s   %12.6f  %12.6f  %12.6f\n' % (at[0],at[1][0],at[1][1],at[1][2]))

    def atoms_centroid(self,n,frame=1):
        
        xyz = self.atom_xyz(n,frame=frame)
        
        natoms = len(xyz)
        return np.sum(np.transpose(np.array(xyz)),axis=1)/natoms


    def atoms_baricenter(self,n,frame=1):
        from MyPy.DATA.atomdata import atomdata
        data = atomdata()

        xyz = self.atom_xyz(n,frame=frame)
        at_weight = map(lambda x: data.get_weight(x),self.atom_type(n))
#        print np.array(at_weight).T
#        print np.array(xyz).T
#        print np.array(at_weight) * np.array(xyz).T

        baricenter = np.sum(np.array(at_weight) * np.array(xyz).T,axis=1)/np.sum(np.array(at_weight))
        
        return baricenter


    def atoms_distance(self,a1,a2,frame=1):
        
        xyz = np.array(self.atom_xyz(str(a1)+','+str(a2),frame=frame))

        return np.linalg.norm(xyz[0,:] - xyz[1,:])


    def atoms_angle(self,a1,a2,a3,frame=1):
        
        
        d12 = self.atoms_distance(a1,a2,frame=frame)
        d13 = self.atoms_distance(a1,a3,frame=frame)
        d23 = self.atoms_distance(a3,a2,frame=frame)

        return np.arccos((d12**2 + d23**2 - d13**2)/(2 * d12 * d23)) * (180.0/self.pi)


    def atoms_angle_4(self,a1,a2,b1,b2,frame=1):

        
        atoms_list = str(a1)+','+str(a2)+','+str(b1)+','+str(b2)
        xyz = np.array(self.atom_xyz(atoms_list,frame=frame))

        u = xyz[0,:]-xyz[1,:]
        v = xyz[2,:]-xyz[3,:]
        
        cos_alpha = np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v))
        return np.arccos(cos_alpha)*180/self.pi
                       

    
    def atoms_died(self,a1,a2,a3,a4,frame=1):
        
        atoms_list = str(a1)+','+str(a2)+','+str(a3)+','+str(a4)
        xyz = np.array(self.atom_xyz(atoms_list,frame=frame))

        A1 = xyz[0,:]
        X1 = xyz[1,:]
        X2 = xyz[2,:]
        A2 = xyz[3,:]

        
        Ua = np.cross((X1-A1),(X2-A1))
        Ub = np.cross((X1-A2),(X2-A2))

        phi = np.arccos(np.dot(Ua,Ub)/(np.linalg.norm(Ua)*np.linalg.norm(Ub)))

        return phi * 180.0/self.pi

    # def trajectory(self,n,frames='all'):
        
    #     traj = []

    #     if frames == 'all':
    #         frames = '1-'+str(self.frames)

    #     n = n.split(',')
    #     frames = frame.split(',')

    #     for frame in frames:

    #         start_id,end_id = self._checklist(frame)
    #         for f in range(start_id,end_id):
                
                
                
            
####
if __name__ == '__main__':        
            

    import sys
    
    xyz = mdxyz('asd.xyz')
    xyz.read_xyz()
    print xyz.atom('0',frame=1)

    exit()


    # points = [[1.,4.,5.],[4.,4.,2.],[2.,-4.,1.],[-5,-5.,0.]]
    # gplane = [0,0,2,0]
    # print len(points)

    # plane = meanPlane(points,gplane)
    # final = plane.optimize_plane()
    # print final
    # plane.show_graph()

    system = '4T-4mon'


    if system == '4T-4mon':
        ## 4T-4mon
        xyz = mdxyz('./asd.xyz')

        xyz.read_xyz()
        dihed         = []
        c_alpha       = []
        c_alphavcore1 = []
        c_alphavcore2 = []
        pipi_dist     = []
        c_alphavpept  = []
        for f in xrange(0,xyz.frames):
        #    dihed.append(xyz.atoms_died(197,171,337,363,frame=f))
    #        dihed.append(xyz.atoms_died(533,507,337,353,frame=f))
    #        dihed.append(xyz.atoms_died(533,507,337,353,frame=f))
            dihed.append(xyz.atoms_died(507,533,363,337,frame=f))
            c_alphavpept.append(xyz.atoms_angle(206,372,427,frame=f))
            c_alpha.append(xyz.atoms_distance(372,206,frame=f))
            c_alpha.append(xyz.atoms_distance(383,217,frame=f))
            c_alphavcore1.append(xyz.atoms_angle_4(383,217,197,171,frame=f))
            c_alphavcore2.append(xyz.atoms_angle_4(363,337,372,206,frame=f))
            #pipi
            plane = meanPlane(xyz.atom_xyz('171-197',frame=f),[1,1,1,1])
            other_core = xyz.atom_xyz('337-363',frame=f)
            plane.other_points = other_core
            centroid = xyz.atoms_centroid('337-363',frame=f)
            final = plane.optimize_plane()
            pipi_dist.append(plane.point_from_plane(centroid))

        c_alpha       = np.array(c_alpha)
        pipi_dist     = np.array(pipi_dist)
        dihed         = np.array(dihed)
        c_alphavcore1 = np.array(c_alphavcore1)
        c_alphavcore2 = np.array(c_alphavcore2)
        c_alphavpept  = np.array(c_alphavpept)




    if system == '4T-6mon':
        ## 4T-4mon
        xyz = mdxyz('./asd.xyz')

        xyz.read_xyz()
        dihed         = []
        c_alpha       = []
        c_alphavcore1 = []
        c_alphavcore2 = []
        pipi_dist     = []
        c_alphavpept  = []
        for f in xrange(80,xyz.frames):
        #    dihed.append(xyz.atoms_died(197,171,337,363,frame=f))
            dihed.append(xyz.atoms_died(507,533,367,341,frame=f))
        #    dihed.append(xyz.atoms_died(507,533,363,337,frame=f))
            c_alphavpept.append(xyz.atoms_angle(387,553,567,frame=f))
            c_alpha.append(xyz.atoms_distance(387,553,frame=f))
            c_alpha.append(xyz.atoms_distance(376,542,frame=f))
            c_alphavcore1.append(xyz.atoms_angle_4(341,367,387,401,frame=f))
            c_alphavcore2.append(xyz.atoms_angle_4(533,507,542,597,frame=f))
            #pipi
            plane = meanPlane(xyz.atom_xyz('341-368',frame=f),[1,1,1,1])
    #        other_core = xyz.atom_xyz('534-507',frame=f)
    #        plane.other_points = other_core
            centroid = xyz.atoms_centroid('507-534',frame=f)
            final = plane.optimize_plane()
            pipi_dist.append(plane.point_from_plane(centroid))

        c_alpha       = np.array(c_alpha)
        pipi_dist     = np.array(pipi_dist)
        dihed         = np.array(dihed)
        c_alphavcore1 = np.array(c_alphavcore1)
        c_alphavcore2 = np.array(c_alphavcore2)
        c_alphavpept  = np.array(c_alphavpept)




    if system == '4T-4mon-PM6D3':
        ## 4T-4mon
        xyz = mdxyz('./asd.xyz')

        xyz.read_xyz()
        dihed         = []
        c_alpha       = []
        c_alphavcore1 = []
        c_alphavcore2 = []
        pipi_dist     = []
        c_alphavpept  = []
        for f in xrange(80,xyz.frames):
        #    dihed.append(xyz.atoms_died(197,171,337,363,frame=f))
            dihed.append(xyz.atoms_died(507,533,367,341,frame=f))
        #    dihed.append(xyz.atoms_died(507,533,363,337,frame=f))
            c_alphavpept.append(xyz.atoms_angle(387,553,567,frame=f))
            c_alpha.append(xyz.atoms_distance(387,553,frame=f))
            c_alpha.append(xyz.atoms_distance(376,542,frame=f))
            c_alphavcore1.append(xyz.atoms_angle_4(341,367,387,401,frame=f))
            c_alphavcore2.append(xyz.atoms_angle_4(533,507,542,597,frame=f))
            #pipi
            plane = meanPlane(xyz.atom_xyz('341-368',frame=f),[1,1,1,1])
    #        other_core = xyz.atom_xyz('534-507',frame=f)
    #        plane.other_points = other_core
            centroid = xyz.atoms_centroid('507-534',frame=f)
            final = plane.optimize_plane()
            pipi_dist.append(plane.point_from_plane(centroid))

        c_alpha       = np.array(c_alpha)
        pipi_dist     = np.array(pipi_dist)
        dihed         = np.array(dihed)
        c_alphavcore1 = np.array(c_alphavcore1)
        c_alphavcore2 = np.array(c_alphavcore2)
        c_alphavpept  = np.array(c_alphavpept)



    print xyz.frames
    print 'c_alpha dist: ',np.average(c_alpha),' +-',np.std(c_alpha)
    print 'pipi dist: ',np.average(pipi_dist),' +-',np.std(pipi_dist)
    print 'c_alpha ^ core: ', 180 - np.average(np.concatenate((c_alphavcore1,c_alphavcore2))),' +-',np.std(np.concatenate((c_alphavcore1,c_alphavcore2)))
    print 'dihed: ',np.average(dihed),' +-',np.std(dihed)
    print 'c_alpha ^ pept: ',np.average(c_alphavpept),' +-',np.std(c_alphavpept)



