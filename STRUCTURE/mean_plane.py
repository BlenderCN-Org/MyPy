#!/usr/bin/env python

import numpy as np
import sys

sys.path.append('/home/petragli/MyCodes')


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


if __name__ == '__main__':
    print 'Write some working example'
    print 'For now refer to the xyz_trajsplit.py comments'
    sys.exit()
