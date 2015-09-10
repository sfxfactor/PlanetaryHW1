import numpy as np
import matplotlib.pyplot as plt
from orbit import *
from scipy.signal import argrelextrema

AU = 1.49597871e13 #cm
Msun = 1.989e33 #g
Mjup = 1.89813e30 #g
Day = 86164.0916 #seconds in a day
Rjup = 7.1492e9 #radius of jupiter in cm
Rs = 6.9599e10 #radius of sun in cm

#####Q2

# HD 80606 b orbit parameters

t0 = 2454424.852*Day #JD
a = 0.449 * AU
e = 0.9332
i = np.radians(89.32)
W = np.radians(160.98) # or -19.02 http://adsabs.harvard.edu/abs/2014AAS...22341102W
w = np.radians(300.80)
m1 = 0.97 * Msun
m2 = 3.94 * Mjup

t1 = 2457235.5 #JD of Sep 1 2015 
t2 = 2457387.5 #JD of Dec 31 2015

#calculate range of time in seconds
JDs = np.arange(t1,t2,0.001)
ts = JDs*Day

#generate model
p=[t0,a,e,i,W,w,m1,m2]
HD80606b=Orbit(p)
obs=HD80606b.calcObs(ts)
Rvs=obs[6]
#Rvs=np.append(Rvs,Rv)

#plot RV
plt.clf()
plt.plot(JDs ,Rvs/100.)
plt.xlabel("Julian date [days]")
plt.ylabel("Radial Velocity [m/s]")
plt.title("HD 80606 b")
plt.savefig("Q2.pdf")

#write file with important values
f = open('results.out','w')
f.write("Extreme neg vel on: {}\n".format(JDs[argrelextrema(Rvs,np.less)]))
f.write("Extreme pos vel on: {}\n".format(JDs[argrelextrema(Rvs,np.greater)]))

'''
plt.clf()
plt.plot(JDs ,obs[0]/AU)
plt.xlabel("Julain date [days]")
plt.ylabel("Projected seperation [AU]")
plt.title("HD 80606 b")
plt.savefig("Q3.pdf")
'''

##### Q3 #####
f.write("minimum projected seperation on JD {}\n".format(JDs[argrelextrema(obs[0],np.less)]))

rpl = 0.98*Rjup
rs = 0.978*Rs
f.write("1/4 contact occurs on JD {}\n".format(JDs[argrelextrema(np.abs(obs[0]-rpl-rs),np.less)]))
f.write("2/3 Contact occurs on JD {}\n".format(JDs[argrelextrema(np.abs(obs[0]+rpl-rs),np.less)]))

f.close()


