import numpy as np

G = 6.67259e-8 #Gravitational constant (cm^3/g/s^2)

def calcObs(p):
    ''' Calculates observables from orbital parameters.
    returns array of observables [sep, PA, R1, PA1, R2, PA2, RV].
    Rx and PAx are the seperation (in the same units as a) with respect to the barycenter.
    :param p:
    array of orbital parameters. [t, a, e, i, W, w, m1, m2]
    assuming t0 = 0
    
    ##### all angles are in radians #####
    '''

    #pull out orbital parameters
    t = p[0]
    a = p[1]
    e = p[2]
    i = p[3]
    W = p[4]
    w = p[5]
    m1 = p[6]
    m2 = p[7]

    P = np.sqrt((4.*a**3*np.pi**2)/(G*(m1+m2)))

    M = 2*np.pi*t/P % (2.*np.pi)

    #set up transcendental equation and derivative to find eccentric anomaly
    def g(E):
        return E-e*np.sin(E)-M
    def gp(E):
        return 1-e*np.cos(E)

    E = NRmethod(g, gp, 1000, M)
    r=a*(1.-e*np.cos(E))

    f = 2.*np.arctan(np.sqrt((1.+e)/(1.-e))*np.tan(E/2.))

    #coordinates in sky plane
    X = r*(np.cos(W)*np.cos(w+f)-np.sin(W)*np.sin(w+f)*np.cos(i))
    Y = r*(np.sin(W)*np.cos(w+f)+np.cos(W)*np.sin(w+f)*np.cos(i))
    Z = r*np.sin(w+f)*np.sin(i)

    #position angle west of north
    PA = np.arctan2(Y,X) % (2.*np.pi)
    sep = np.sqrt(X**2+Y**2)

    R1 = (m2/(m1+m2))*sep
    R2 = (m1/(m1+m2))*sep

    PA1 = (PA + np.pi) % (2.*np.pi)
    PA2 = PA

    RV = (m2*2*np.pi*a*np.sin(i))/((m1+m2)*P*np.sqrt(1-e**2))*(np.cos(w+f)+e*np.cos(w))

    return [sep, PA, R1, PA1, R2, PA2, RV]
    
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

    if (f(xip1) > 0.1).any():
        print "### WARNING: found root to be ", f(xip1)
    return xip1

