from scipy import special
import numpy as np
import astropy.constants as const
import astropy.units as units
from math import *
import matplotlib.pyplot as plt
from matplotlib import rc
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rc('text', usetex=True)
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18

temp=5000.
elpress=1e2
c=2.99792*1e10
kerg = 1.380658e-16
kev = 8.61734e-5
elmass = 9.10939e-28
h = 6.62607e-27 #[ergs]
keVT = kev*temp
kergT = kerg*temp
eldens = elpress/kergT


def planck(temp,wav):
	B=(2*h*c**2*wav**-5)*(1./(np.exp(h*c/(wav*kerg*temp))-1.))
	return B

#----------3.1
"""
wav = np.arange(1000,20801,200)
b = np.zeros(wav.shape)

plt.xlabel(r'wavelength $\lambda / \AA$', size=18)
plt.ylabel(r'Planck function', size=18)

#plt.yscale('log')
#plt.xscale('log')
plt.xlim(1e3,2.1e4)
for T in range(8000,5000-1,-200):
	b[:] = planck(T, wav[:]*1e-8)
	plt.plot(wav,b,'-')	
	plt.hold('on')
plt.savefig("planckcurvesnonlog")
plt.show()
"""
#---------3.2 Radiation through an isothermal layer
"""
B = 2.
tau = np.arange(0.01,10.01, 0.01)
intensity = np.zeros(tau.shape)
for I0 in range(4,-1,-1):
	intensity[:] = I0 * np.exp(-tau[:]) + B*(1-np.exp(-tau[:]))
	plt.plot(tau, intensity, label = 'intensity I0 = ' + str(I0))

#plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'optical depth $\tau$', size=18)
plt.ylabel(r'intensity', size=18)
plt.legend(fontsize=16, loc='lower right')

plt.savefig("opticaldepth")
plt.show()
"""
# -------- 3.3 spectral lines from a solar reversing layer -------
# This function computes the voigt profile in python.
# one needs to introduce the two parameters for the voigt function.
# contains faddeeva function with which voight profile can evaluated
# spot-checked with idl voigt routine
def voigt(gamma,x):
	z = (x+1j*gamma)
	V = special.wofz(z).real 
	return V
"""
u = np.arange(-10,10.1,0.1)
a = np.array([0.001,0.01,0.1,1])
vau = np.zeros((a.shape[0],u.shape[0]))

for i in range(4):
	vau[i,:] = voigt(a[i],u[:])
	plt.plot(u[:],vau[i,:], label = 'a = ' + np.str(a[i]))

plt.ylim(0,1)
plt.xlim(-10,10)
plt.legend(fontsize=16)
plt.ylabel('voigt profile', size=18)
plt.savefig("voigt")
plt.show()


for i in range(4):
	vau[i,:] = voigt(a[i],u[:])
	plt.plot(u[:],vau[i,:], label = 'a = ' + np.str(a[i]))

plt.yscale('log')
plt.legend(fontsize=16, loc = 8)
plt.xlabel('u', size=18)
plt.ylabel('logarithmic voigt profile', size=18)
plt.savefig("voigtlog")
plt.show()
"""
#---------Schuster-Schwarzchild line profile-------------

Ts = 4200. # solar surface temperature
Tl = 5700. # solar T-min temperature = 'reversing layer'
a = 0.1 # damping parameter
wav = 5000.0e-8 # wavelength in cm
tau0 = 1. # reversing layer thickness at line center
u = np.arange(-10,10.1,0.1)
intensity = np.zeros(u.shape)

for i in range(201):
	tau = tau0 * voigt(a, u[i])
	intensity[i] = planck(Ts,wav) * np.exp(-tau) + planck(Tl,wav)*(1.-np.exp(-tau))
"""
plt.plot(u,intensity)
plt.show()
"""


logtau0 = np.arange(-2,2.1,0.5)

for itau in range(9):
	for i in range(201):
		tau = 10.**(logtau0[itau]) * voigt(a, u[i])
		intensity[i] = planck(Ts,wav) * np.exp(-tau) + planck(Tl,wav)*(1.-np.exp(-tau))
	plt.plot(u,intensity, label = r'$\log{(\tau_0)} = $' + np.str(logtau0[itau]))

plt.legend(loc=3, fontsize=16)
plt.savefig("lineprofilesE")
plt.show()




for iwav in range(1,4):
	wav = (iwav**2+1.)*1.0e-5 # wav = 2000, 5000, 10000 angstrom

for itau in range(8):
	for i in range(201):
		tau = 10.**(logtau0[itau]) * voigt(a,u[i])
		intensity[i] = planck(Ts,wav) * exp(-tau) + planck(Tl,wav)*(1.-exp(-tau))
	intensity = intensity / intensity[0]
	plt.plot(u,intensity[:],  linewidth=1.)
plt.savefig("lineprofilesscaledE")
plt.show()

#------ 3.4 The equivalent width of spectral lines

def profile(a,tau0,u):
	Ts = 4200.#5700.
	Tl = 5700.
	wav = 2000.0e-8
	intensity = np.zeros(u.size)
	usize = u.size
	for i in range(usize):
		tau = tau0 * voigt(a, u[i])
		intensity[i] = planck(Ts,wav)*np.exp(-tau) + planck(Tl,wav)*(1.-np.exp(-tau))
	return intensity
# Checking the profile
u = np.arange(-200,200.4,0.4)
a = 0.1
tau0 = 1.0e2

intensity = profile(a,tau0,u)

plt.plot(u,intensity)
plt.show()


# relative
reldepth = (intensity[0]-intensity)/intensity[0]
plt.plot(u,reldepth)
plt.show()
eqw = sum(reldepth)*0.4
#print eqw

#-----   3.5 The curve of growth

tau0 = np.logspace(-2, 4, 61)
eqw = np.zeros(tau0.size)
for i in range(61):
	intensity = profile(a,tau0[i],u)
	reldepth = (intensity[0] - intensity) / intensity[0]
	eqw[i] = sum(reldepth)*0.4

plt.plot(tau0,abs(eqw))
plt.xlabel(r'$\tau_0$', size=18)
plt.ylabel(r'abs(equivalent width) |$W_{\lambda}|$', size=18)
plt.xscale('log')
plt.yscale('log')
plt.savefig("curveofgrowthE")
plt.show()


