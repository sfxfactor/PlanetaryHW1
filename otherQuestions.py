import numpy as np
import matplotlib.pyplot as plt
import orbits as orb

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

t1 = 2457266. #JD of Sep 1 2015 
t2 = 2457387. #JD of Dec 31 2015

#calculate range of time in seconds
JDs = np.arange(t1,t2,0.001)
ts = JDs*Day

p=[t0,a,e,i,W,w,m1,m2]
obs=orb.calcObs(p,ts)
Rvs=obs[6]
#Rvs=np.append(Rvs,Rv)

plt.clf()
plt.plot(JDs ,Rvs/100.)
plt.xlabel("Julain date [days]")
plt.ylabel("Radial Velocity [m/s]")
plt.title("HD 80606 b")
plt.savefig("Q2.pdf")

plt.clf()
plt.plot(JDs ,obs[0]/AU)
plt.xlabel("Julain date [days]")
plt.ylabel("Projected seperation [AU]")
plt.title("HD 80606 b")
print "primary transit on JD "+str(JDs[np.where(obs[6]>0.)][np.where(obs[0]==np.min(obs[0]))])
print "second transit on JD "+str(JDs[np.where(obs[6]<0.)][np.where(obs[0]==np.min(obs[0]))])
plt.savefig("Q3.pdf")




