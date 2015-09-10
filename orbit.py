import numpy as np
import sys

class Orbit:
    G = 6.67259e-8 #Gravitational constant (cm^3/g/s^2)

    def __init__(self,p):
        '''
        Creates an Orbit object

        :param p:
        array of orbital parameters. [t0, a, e, i, W, w, m1, m2]
        
        ##### all angles are in radians and all other units are cgs (ARG WHY?!?!) #####
        '''
        
        #pull out orbital parameters
        self.t0 = p[0]
        self.a = p[1]
        self.e = p[2]
        self.i = p[3]
        self.W = p[4]
        self.w = p[5]
        self.m1 = p[6]
        self.m2 = p[7]
        self.P = np.sqrt((4.*self.a**3*np.pi**2)/(Orbit.G*(self.m1+self.m2)))

    def get_params(self):
        '''
        Returns orbital parameters
        '''
        return [self.t0,self.a,self.e,self.i,self.W,self.w,self.m1,self.m2,self.P]

    def calcCoord(self, t):
        '''
        Calculates and returns the cartesian and radial coordinates.

        :param t:
        Single time or array of times to compute the coordinates of.
        '''

        M = 2*np.pi*(t-self.t0)/self.P % (2.*np.pi)

        #set up transcendental equation and derivative to find eccentric anomaly
        def g(E):
            return E-self.e*np.sin(E)-M
        def gp(E):
            return 1-self.e*np.cos(E)

        E = NRmethod(g, gp, 1000, M)

        f = 2.*np.arctan(np.sqrt((1.+self.e)/(1.-self.e))*np.tan(E/2.))

        r=self.a*(1.-self.e*np.cos(E))

        #coordinates in sky plane
        X = r*(np.cos(self.W)*np.cos(self.w+f)-np.sin(self.W)*np.sin(self.w+f)*np.cos(self.i))
        Y = r*(np.sin(self.W)*np.cos(self.w+f)+np.cos(self.W)*np.sin(self.w+f)*np.cos(self.i))
        Z = r*np.sin(self.w+f)*np.sin(self.i)

        return [X,Y,Z,r,f,E,M]

    def calcObs(self, t, coords=None):
        ''' Calculates observables from orbital parameters. If orbital parameters are not provided, this routine calculates them from the given times. 
        returns array of observables and coordinates (see calcCoord) [sep, PA, R1, PA1, R2, PA2, RV, coords].
        Rx and PAx are the seperation with respect to the barycenter.

        :param t:
        Single time or array of time to compute the coordinates of.

        :param coords:
        Coordinates (in same format as the output of calcCoord) of the system at the relavent times. If this parameter is provided, the given coordinates are used rather than recalculated.
        '''

        if coords is None:
            X, Y, Z, r, f, E, M = self.calcCoord(t)
        else:
            X, Y, Z, r, f, E, M = coords

        #position angle east of north
        PA = np.arctan2(Y,X) % (2.*np.pi)
        sep = np.sqrt(X**2+Y**2)

        R1 = (self.m2/(self.m1+self.m2))*sep
        R2 = (self.m1/(self.m1+self.m2))*sep

        PA1 = (PA + np.pi) % (2.*np.pi)
        PA2 = PA

        RVoM = (2*np.pi*self.a*np.sin(self.i))/((self.m1+self.m2)*self.P*np.sqrt(1-self.e**2))*(np.cos(self.w+f)+self.e*np.cos(self.w))

        if coords is None:
            return [sep, PA, R1, PA1, R2, PA2, RVoM*self.m2, RVoM*self.m1, X, Y, Z, r, f, E, M]
        else:
            return np.append([sep, PA, R1, PA1, R2, PA2, RVoM*self.m2, RVoM*self.m1],coords)
        
def NRmethod(f, fp, n, x0):
    ''' Newton Raphson root finding method.
    returns position of root.
    :param f:
    Function to find the root of.

    :param fp:
    First derivative of f.

    :param n:
    Number of iterations to run.

    :param x0:
    Initial guess.
    '''

    xi = x0


    for i in xrange(n):
        xip1 = xi - f(xi)/fp(xi)
        xi=xip1
        percent = float(i) / n
        hashes = '#' * int(round(percent * 20))
        spaces = ' ' * (20 - len(hashes))
        sys.stdout.write("\rPercent: [{0}] {1}%".format(hashes + spaces, int(round(percent * 100))))
        sys.stdout.flush()

    if (f(xip1) > 0.1).any():
        print "### WARNING: found root to be ", f(xip1)
    else: print ""
    return xip1

