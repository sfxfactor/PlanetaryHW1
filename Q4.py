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

#generate model
p=[t0,a,e,i,W,w,m1,m2]
HD80606b=Orbit(p)

##### Q4 #####

pi=0.01713 #arcsec \pm 0.00577
D=1./pi #distance in PC
gp = 6.7e-6 #gaia precision in arcsec
err = gp*np.ones(100)

dt = 5.*365.256*Day #time window in seconds
samples = np.random.random_sample(size=100)*dt 
samples = np.sort(samples) #100 randomly sampled times
alltimes = np.linspace(0.,1.,num=1000)*dt + 2456863.5*Day #finer resolution for true values
alltimesJDs = alltimes/Day
epochs = samples + 2456863.5*Day #shift begining of observations to July 25 2014 Gaia's first science day
JDs = epochs/Day

#calculate observables
Gobs=HD80606b.calcObs(epochs)
obs=HD80606b.calcObs(alltimes)
GRcm, GPAcm = [Gobs[2], Gobs[3]]
Rcm, PAcm = [obs[2], obs[3]]
GDRA, GDDEC = (GRcm/(AU*D))*[np.sin(GPAcm),np.cos(GPAcm)]
DRA, DDEC = (Rcm/(AU*D))*[np.sin(PAcm),np.cos(PAcm)]

# Gaussian random noise
GRAnoise = np.random.normal(loc=0.0, scale=gp, size=100)
GDECnoise = np.random.normal(loc=0.0, scale=gp, size=100)
GDRA += GRAnoise
GDDEC += GDECnoise

#plot reflex motion
plt.clf()
f, axarr = plt.subplots(2, sharex=True)
axarr[0].errorbar(JDs,GDRA*1e6,yerr=err*1e6,fmt='k.')
axarr[0].plot(alltimesJDs,DRA*1e6)
axarr[0].set_title('Reflex motion')
axarr[0].set_ylabel(r"$\Delta\alpha$ [$\mu as$]")
axarr[1].errorbar(JDs,GDDEC*1e6,yerr=err*1e6,fmt='k.')
axarr[1].plot(alltimesJDs,DDEC*1e6)
axarr[1].set_ylabel(r"$\Delta\delta$ [$\mu as$]")
plt.xlabel("Julian date [days]")
plt.savefig("1time.pdf")

plt.clf()
plt.errorbar(GDRA*1e6,GDDEC*1e6,xerr=err*1e6,yerr=err*1e6,fmt='k.')
plt.plot(DRA*1e6,DDEC*1e6)
plt.title('Reflex motion')
plt.xlabel(r"$\Delta\alpha$ [$\mu as$]")
plt.ylabel(r"$\Delta\delta$ [$\mu as$]")
plt.savefig('1radec.pdf')

#J2000. coordinates
a = np.radians((9+22./60.+37.568/3600.)*15.) #RA
d = np.radians(50+36./60.+13.48/3600.) #DEC
e = np.radians(23.44) #obliquity

#calculate longitude of the sun for gaia observations
n = JDs - 2451545. #days since J2000.0
L = (280.46 + 0.9856474*n) % 360. #mean longitude of sun
g = np.radians(357.528 + 0.9856003*n) % (2.*np.pi) #mean anomaly of sun
theta = np.radians(L + 1.915*np.sin(g) + 0.020*np.sin(2.*g)) #ecliptic longitude of sun
#same calculation for full resolution curve
an = alltimesJDs - 2451545.
aL = (280.46 + 0.9856474*an) % 360.
ag = np.radians(357.528 + 0.9856003*an) % (2.*np.pi)
atheta = np.radians(aL + 1.915*np.sin(ag) + 0.020*np.sin(2.*ag))

#calculate parallactic motion
GDRA = GDRA - pi*(np.cos(theta)*np.sin(a)-np.sin(theta)*np.cos(e)*np.cos(a))*np.sin(d)
GDDEC = GDDEC - pi*(np.cos(e)*np.sin(a)*np.sin(d)-np.sin(e)*np.cos(d))-pi*np.cos(a)*np.sin(d)*np.cos(theta)

DRA = DRA - pi*(np.cos(atheta)*np.sin(a)-np.sin(atheta)*np.cos(e)*np.cos(a))*np.sin(d)
DDEC = DDEC - pi*(np.cos(e)*np.sin(a)*np.sin(d)-np.sin(e)*np.cos(d))-pi*np.cos(a)*np.sin(d)*np.cos(atheta)

plt.clf()
f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(JDs,GDRA*1000.,'k.')
axarr[0].plot(alltimesJDs,DRA*1000.)
axarr[0].set_ylabel(r"$\Delta\alpha$ [mas]")
axarr[0].set_title('Reflex + parallax motion')
axarr[1].plot(JDs,GDDEC*1000.,'k.')
axarr[1].plot(alltimesJDs,DDEC*1000.)
axarr[1].set_ylabel(r"$\Delta\delta$ [mas]")
plt.xlabel("Julian date [days]")
plt.savefig("2time.pdf")

plt.clf()
plt.plot(GDRA*1000.,GDDEC*1000.,'k.')
plt.plot(DRA*1000.,DDEC*1000.)
plt.title('Reflex + parallax motion')
plt.xlabel(r"$\Delta\alpha$ [mas]")
plt.ylabel(r"$\Delta\delta$ [mas]")
plt.savefig('2radec.pdf')

#add propper motion from hipparcos
GDRA = GDRA + 0.04698*n/365.24
GDDEC = GDDEC + 0.00692*n/365.24
DRA = DRA + 0.04698*an/365.24
DDEC = DDEC + 0.00692*an/365.24

plt.clf()
f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(JDs,GDRA*1000.,'k.')
axarr[0].plot(alltimesJDs,DRA*1000.)
axarr[0].set_ylabel(r"$\Delta\alpha$ [mas]")
axarr[0].set_title('Reflex + parallax + proper motion')
axarr[1].plot(JDs,GDDEC*1000.,'k.')
axarr[1].plot(alltimesJDs,DDEC*1000.)
axarr[1].set_ylabel(r"$\Delta\delta$ [mas]")
plt.xlabel("Julian date [days]")
plt.savefig("3time.pdf")

plt.clf()
plt.plot(GDRA*1000.,GDDEC*1000.,'k.')
plt.plot(DRA*1000.,DDEC*1000.)
plt.title('Reflex + parallax + proper motion')
plt.xlabel(r"$\Delta\alpha$ [mas]")
plt.ylabel(r"$\Delta\delta$ [mas]")
plt.savefig('3radec.pdf')
