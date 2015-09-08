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

# HD 80606 b orbit parameters from http://exoplanet.eu/catalog/hd_80606_b/#6832
# updated Aug. 16, 2014

t0 = 2454424.852*Day #JD
a = 0.449 * AU
e = 0.9332
i = np.radians(89.32)
W = np.radians(160.98) # or -19.02 http://adsabs.harvard.edu/abs/2014AAS...22341102W
w = np.radians(300.80)
m1 = 0.97 * Msun
m2 = 3.94 * Mjup

p=[t0,a,e,i,W,w,m1,m2]
HD80606b=Orbit(p)

##### Q4 #####

pi=0.01713 #arcsec \pm 0.00577
D=1./pi #distance in PC
gp = 6.7e-6 #gaia precision in arcsec
err = gp*np.ones(100)

dt = 5.*365.256*Day
samples = np.random.random_sample(size=100)*dt 
samples = np.sort(samples)
epochs = samples + 2456863.5*Day #July 25 2014 Gaia's first science day
JDs = epochs/Day

Gobs=HD80606b.calcObs(epochs)
Rcm, PAcm = [Gobs[2], Gobs[3]]
DRA, DDEC = (Rcm/(AU*D))*[-np.cos(PAcm),np.sin(PAcm)]

## Gaussian random noise
RAnoise = np.random.normal(loc=0.0, scale=gp, size=100)
DECnoise = np.random.normal(loc=0.0, scale=gp, size=100)
DRA += RAnoise
DDEC += DECnoise

plt.clf()
f, axarr = plt.subplots(2, sharex=True)
axarr[0].errorbar(JDs,DRA*1e6,yerr=err*1e6,fmt='k.')
axarr[0].set_title('Reflex motion')
axarr[0].set_ylabel(r"$\Delta\alpha$ [$\mu as$]")
axarr[1].errorbar(JDs,DDEC*1e6,yerr=err*1e6,fmt='k.')
axarr[1].set_ylabel(r"$\Delta\delta$ [$\mu as$]")
plt.xlabel("Julian date [days]")
plt.savefig("1time.pdf")

plt.clf()
plt.errorbar(DRA*1e6,DDEC*1e6,xerr=err*1e6,yerr=err*1e6,fmt='k.')
plt.title('Reflex motion')
plt.xlabel(r"$\Delta\alpha$ [$\mu as$]")
plt.ylabel(r"$\Delta\delta$ [$\mu as$]")
plt.savefig('1radec.pdf')

a = np.radians((9+22./60.+37.568/3600.)*15.)
d = np.radians(50+36./60.+13.48/3600.)
e = np.radians(23.44)
n = JDs - 2451545.
L = (280.46 + 0.9856474*n) % 360.
g = np.radians(357.528 + 0.9856003*n) % (2.*np.pi)
theta = np.radians(L + 1.915*np.sin(g) + 0.020*np.sin(2.*g))

DRA = DRA - pi*(np.cos(theta)*np.sin(a)-np.sin(theta)*np.cos(e)*np.cos(a))*np.sin(d)
DDEC = DDEC - pi*(np.cos(e)*np.sin(a)*np.sin(d)-np.sin(e)*np.cos(d))-pi*np.cos(a)*np.sin(d)*np.cos(theta)

plt.clf()
f, axarr = plt.subplots(2, sharex=True)
axarr[0].errorbar(JDs,DRA*1000.,yerr=err*1000.,fmt='k.')
axarr[0].set_ylabel(r"$\Delta\alpha$ [mas]")
axarr[0].set_title('Reflex + parallax motion')
axarr[1].errorbar(JDs,DDEC*1000.,yerr=err*1000.,fmt='k.')
axarr[1].set_ylabel(r"$\Delta\delta$ [mas]")
plt.xlabel("Julian date [days]")
plt.savefig("2time.pdf")

plt.clf()
plt.errorbar(DRA*1000.,DDEC*1000.,xerr=err*1000.,yerr=err*1000.,fmt='k.')
plt.title('Reflex + parallax motion')
plt.xlabel(r"$\Delta\alpha$ [mas]")
plt.ylabel(r"$\Delta\delta$ [mas]")
plt.savefig('2radec.pdf')

DRA = DRA + 0.04698*n/365.24
DDEC = DDEC + 0.00692*n/365.24

plt.clf()
f, axarr = plt.subplots(2, sharex=True)
axarr[0].errorbar(JDs,DRA*1000.,yerr=err*1000.,fmt='k.')
axarr[0].set_ylabel(r"$\Delta\alpha$ [mas]")
axarr[0].set_title('Reflex + parallax + proper motion')
axarr[1].errorbar(JDs,DDEC*1000.,yerr=err*1000.,fmt='k.')
axarr[1].set_ylabel(r"$\Delta\delta$ [mas]")
plt.xlabel("Julian date [days]")
plt.savefig("3time.pdf")

plt.clf()
plt.errorbar(DRA*1000.,DDEC*1000.,xerr=err*1000.,yerr=err*1000.,fmt='k.')
plt.title('Reflex + parallax + proper motion')
plt.xlabel(r"$\Delta\alpha$ [mas]")
plt.ylabel(r"$\Delta\delta$ [mas]")
plt.savefig('3radec.pdf')
