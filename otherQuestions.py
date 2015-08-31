import numpy as np
import orbits as orb

AU = 1.49597871e13 #cm
Msun = 1.989e33 #g
Mjup = 1.89813e30 #g
Day = 86164.0916 #seconds in a day

# HD 80606 b orbit parameters from http://exoplanet.eu/catalog/hd_80606_b/#6832
# updated Aug. 16, 2014

Tperi = 2454424.857 #JD
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
ts = (np.arange(t1,t2)-Tperi)*Day
print ts

Rvs=np.array([])

for t in ts:
    p=[t,a,e,i,W,w,m1,m2]
    obs=orb.calcObs(p)
    Rv=obs[6]
    Rvs=np.append(Rvs,Rv)
print Rvs
