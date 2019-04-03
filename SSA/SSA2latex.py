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


kerg = 1.380658e-16
keV = 8.61734e-5
elmass = 9.10939e-28
h = 6.62607e-27 #[ergs]
chiion = np.array([7, 16, 31, 51])
chiion = np.array([6, 12, 51, 67]) 

def partfunc_E(temp):
	u = np.zeros(4)
	for r in range(4):
		for s in range(chiion[r]):
			u[r] = u[r] + np.exp(-s / keV / temp)
	return u 

def boltz_E(temp, r, s):

	u = partfunc_E(temp)
	relnrs = 1. / u[r - 1] * np.exp(-(s - 1) / (keV* temp))
	return relnrs

"""
CORRECT!
for s in range(1,11): # now the loop starts at 1 and finishes at 10
	print boltz_E(5000., 1., s) #DO FOR 3 TEMPS
"""

def saha_E(temp, elpress, ionstage):

	keVT = keV * temp
	kergT = kerg * temp
	eldens = elpress / kergT
	
	u = partfunc_E(temp)
	u = np.append(u, 2) 
	sahaconst = (2. * np.pi * elmass * kergT / (h**2))**1.5 * 2. / eldens
	nstage = np.zeros(5)
	nstage[0] = 1.
	for r in range(4):
		nstage[r + 1] = nstage[r] * sahaconst * u[r + 1] / u[r] * np.exp(-chiion[r] / keVT)
	ntotal = np.sum(nstage)
	nstagerel = nstage / ntotal
	return nstagerel[ionstage - 1]

def sahabolt_E(temp, elpress, ion, level):
	return saha_E(temp, elpress, ion) * boltz_E(temp, ion, level)



def sahabolt_H(temp,elpress,level):

	keVT = keV*temp
	kergT = kerg*temp
	eldens = elpress/kergT

	nrlevels = 100
	g = np.zeros((2,nrlevels))
	chiexc = np.zeros((2,nrlevels))

	for s in range(nrlevels):
		g[0,s] = 2.*(s+1.)**2.
		chiexc[0,s] = 13.598*(1.-1./(s+1.)**2.)	
	g[1,0] = 1.
	chiexc[1,0] = 0. 

	# partition functions
	u = np.zeros([2])
	for s in range(nrlevels):
		u[0] = u[0] + g[0,s]*exp(-chiexc[0,s]/keVT)
	u[1] = g[1,0]
	
	# Saha
	sahaconst = (2*np.pi*elmass*kergT /(h*h))**(1.5)*2./eldens
	nstage = np.zeros(2)
	nstage[0] = 1.
	nstage[1] = nstage[0] * sahaconst * u[1]/u[0] * np.exp(-13.598/keVT)
	ntotal = np.sum(nstage)

	# Boltzmann
	nlevel = nstage[0]*g[0,level-1]/u[0]*np.exp(-chiexc[0,level-1]/keVT)
	nlevelrel = nlevel/ntotal


	"""
	print "-----------"
	for s in range(6):
		print s+1, g[0,s], chiexc[0,s], g[0,s]*np.exp(-chiexc[0,s]/keVT)

	print "-----------"
	for s in range(0,nrlevels,10):
		print s+1, g[0,s], chiexc[0,s], g[0,s]*np.exp(-chiexc[0,s]/keVT)
	"""
	return nlevelrel


"""
temp = np.arange(1000,20001,100)
CaH = np.zeros(temp.shape)
Caabund = 2.0e-6

for i in range(0,191):
	NCa = sahabolt_E(temp[i],1e2,2,1) # is equal to sahabolt_Ca
	NH = sahabolt_H(temp[i],1e2,2)
	CaH[i] = NCa*Caabund/NH


plt.plot(temp,CaH, label=r'strength ratio Ca$^+$K / H$\alpha$')
plt.yscale('log')

plt.xlabel(r'temperature $T / K$', size=14)
plt.ylabel(r'Ca II K / H$\alpha$', size=14)
plt.legend(fontsize=14)
plt.show()

print 'Ca/H ratio at 5000 K = ', CaH[np.argwhere(temp==5000)][0][0]
"""
#-----------------2.9 Solar Ca+ K versus Ha-temperature sensitivity----
"""
temp = np.arange(2000,12001,100)
dNCadT = np.zeros(temp.shape)
dNHdT = np.zeros(temp.shape)
dT = 1.
for i in range(101):
	NCa = sahabolt_E(temp[i],1e2,2,1)
	NCa2 = sahabolt_E(temp[i]-dT,1e2,2,1)
	dNCadT[i] = (NCa - NCa2)/(dT*NCa)
	NH = sahabolt_H(temp[i],1e2,2)
	NH2 = sahabolt_H(temp[i]-dT,1e2,2)
	dNHdT[i] = (NH-NH2)/(dT*NH)


NCa = np.zeros(temp.shape)
NH = np.zeros(temp.shape)

for i in range(101):
	NCa[i] = sahabolt_E(temp[i],1e2,2,1)
	NH[i] = sahabolt_H(temp[i],1e2,2)
plt.figure()
plt.plot(temp,np.absolute(dNHdT), label=r'H')
plt.plot(temp,np.absolute(dNCadT), label=r'Ca$^+$K')


plt.plot(temp,NH/np.amax(NH), ls='--',  label = 'rel. pop. H')
plt.plot(temp,NCa/np.amax(NCa), ls='--', label = r'rel. pop. Ca$^+$')

plt.yscale('log')
plt.ylim(1e-9,1)
plt.xlabel(r'temperature $T/K$', size=18)
plt.ylabel(r"$\left| \left( \Delta n(r,s) / \Delta T \right) /  n(r,s) \right|$", size=18)

plt.legend(loc=4, fontsize=18)
plt.show()
"""
# -------------------- 2.10 Hot stars versus Cool Stars ------------
"""
for T in np.arange(2e3,2e4+1,2e3):
	print T, sahabolt_H(T,1e2,1)

temp = np.arange(1e3,2e4+1,1e2)
nH = np.zeros(temp.shape)
for i in range(191):
	nH[i] = sahabolt_H(temp[i],1e2,1)

plt.plot(temp,nH)
plt.xlabel('temperature $T/K$', size=14)
plt.ylabel('neutral hydrogen fraction', size=14)
plt.show()
"""

