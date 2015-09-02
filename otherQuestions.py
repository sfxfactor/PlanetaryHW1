import numpy as np
import matplotlib.pyplot as plt
from orbit import *
from scipy.signal import argrelextrema

AU = 1.49597871e13 #cm
Msun = 1.989e33 #g
Mjup = 1.89813e30 #g
Day = 86164.0916 #seconds in a day

#####Q2

# HD 80606 b orbit parameters from http://exoplanet.eu/catalog/hd_80606_b/#6832
# updated Aug. 16, 2014

t0 = 2454424.857*Day #JD
a = 0.449 * AU
e = 0.93366
i = np.radians(89.285)
W = np.radians(160.98) # or -19.02 http://adsabs.harvard.edu/abs/2014AAS...22341102W
w = np.radians(300.651)
m1 = 0.98 * Msun
m2 = 3.94 * Mjup

t1 = 2457235.5 #JD of Sep 1 2015 
t2 = 2457387.5 #JD of Dec 31 2015

#calculate range of time in seconds
JDs = np.arange(t1,t2,0.001)
ts = JDs*Day

p=[t0,a,e,i,W,w,m1,m2]
HD80606b=Orbit(p)
obs=HD80606b.calcObs(ts)
Rvs=obs[6]
#Rvs=np.append(Rvs,Rv)

plt.clf()
plt.plot(JDs ,Rvs/100.)
plt.xlabel("Julain date [days]")
plt.ylabel("Radial Velocity [m/s]")
plt.title("HD 80606 b")
plt.savefig("Q2.pdf")
print "Extreme neg vel on: ", JDs[argrelextrema(Rvs,np.less)]
print "Extreme pos vel on: ", JDs[argrelextrema(Rvs,np.greater)]

plt.clf()
plt.plot(JDs ,obs[0]/AU)
plt.xlabel("Julain date [days]")
plt.ylabel("Projected seperation [AU]")
plt.title("HD 80606 b")
print "minimum projected seperation on JD "+str(JDs[argrelextrema(obs[0],np.less)])
plt.savefig("Q3.pdf")




