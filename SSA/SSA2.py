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
temp=5000
elpress=1e3
k = const.k_B
kerg = k.cgs.value
keV = k.to(units.eV/units.K).value	
h = const.h.cgs.value
keVT = keV * temp
kergT = kerg * temp
elmass =  const.m_e.cgs.value 
eldens = elpress / kergT
chiion = np.array([7, 16, 31, 51 ])
# energy levels and weights for hydrogen
nrlevels = 100 # reasonable partition function cut-off value
g = np.zeros((2,nrlevels)) # declarations weights (too many for proton)	3
chiexc = np.zeros((2,nrlevels)) # declaration excitation energies (idem)

#---------2.4 Saha-Boltzmann populations of scadeenium-------------------------
def partfunc_E(temp):
	chiion = np.array([7, 16, 31, 51])# Schadee ionization energies into numpy array
	k = 8.61734e-5	# Boltzmann constant in eV/deg
	u = np.zeros(4)	# declare a 4 zero-element array
	for r in range(4):
		for s in range(chiion[r]):
			u[r] = u[r] + np.exp(-s / k / temp)
	return u # returns all the values of u array

#print partfunc_E(temp)
#print partfunc_E(20000)

def boltz_E(temp, r, s):
	u = partfunc_E(temp)
	KeV = 8.61734e-5 # This constant does need to be defined here again if it was before
	relnrs = 1. / u[r - 1] * np.exp(-(s - 1) / (KeV * temp))
	return relnrs
"""
for s in range(1,11):
	print boltz_E(temp,1, s)


Temp = 5000
0.90181507784
0.0885447147611
0.00869376295073
0.000853597128269
8.38104353107e-05
8.22892771584e-06
8.07957280039e-07 --
7.93292867443e-08
7.78894613717e-09
7.64757688081e-10

Temp = 10000
0.686858855862
0.21522369225
0.0674392377861
0.0211317385442
0.00662152166249
0.00207481978045
0.000650134114296 --
0.000203716183234
6.38334189807e-05
2.00018737543e-05

Temp = 20000
0.447942229527
0.250745585061
0.140360395343
0.078569840327
0.0439812084735
0.0246194556427
0.0137812856259  --
0.00771437988955
0.00431829501949
0.00241726128896

"""


def saha_E(temp, elpress, ionstage):
	k = const.k_B
	kerg = k.cgs.value
	keV = k.to(units.eV/units.K).value	
	h = const.h.cgs.value
	keVT = keV * temp
	kergT = kerg * temp
	elmass =  const.m_e.cgs.value 
	eldens = elpress / kergT
	chiion = np.array([7, 16, 31, 51 ])

	u = partfunc_E(temp)
	u = np.append(u, 2)
	# With this command we are adding a new element to the array
	sahaconst = (2. * np.pi * elmass * kergT / (h**2))**1.5 * 2. / eldens	
	nstage = np.zeros(5)
	nstage[0] = 1.	# We set the first element of the array to a value 1
	for r in range(4):
		nstage[r + 1] = nstage[r] * sahaconst * u[r + 1] / u[r] * np.exp(-chiion[r] / keVT)
	ntotal = np.sum(nstage)
	nstagerel = nstage / ntotal
	return nstagerel[ionstage - 1]

"""
for r in range(1,6):
	print saha_E(20000,1e3,r)
for r in range(1,6):
	print saha_E(20000,1e1,r)

2.72791311909e-10
0.000180285342838
0.632013153927
0.367804840818
1.7196401184e-06

7.28809862609e-16
4.81663932127e-08
0.0168853405432
0.982655179254
0.000459432036035
"""


# -----------2.5 PAYNE CURVES FOR SHADEENIUM-------------------------------

def sahabolt_E(T,elpress,ion,level):
	return saha_E(T,elpress,ion)*boltz_E(T,ion,level)
"""
for s in range(1,6):
	print sahabolt_E(5000,1e3,1,s)
for s in range(1,6):
	print sahabolt_E(20000,1e3,1,s)
for s in range(1,6):
	print sahabolt_E(10000,1e3,2,s)
for s in range(1,6):
	print sahabolt_E(20000,1e3,4,s)

0.817095991836
0.0802265711757
0.0078770460104
0.000773407774266
7.59370434685e-05

1.22194748452e-10
6.84012171041e-11
3.82890963858e-11
2.14331698193e-11
1.19976915588e-11

0.648955432457
0.203346849338
0.0637176901026
0.0199656107052
0.00625612149764

0.161917934159
0.0906371501867
0.0507361524627
0.0284006851652
0.0158979126067
"""
"""
What causes the steep flanks on the left and right side of each peak?
What happens for T\/0 T/\inf?
"""

"""
temp = np.arange(0,30001,1000)
pop = np.zeros((5,31))
for T in np.arange(1,31):
	for r in np.arange(1,5):
		pop[r,T] = sahabolt_E(temp[T],131,r,1) #Change s for energy level

colors = ['-b', '-y', '-g', '-r']
labellst = ['ground stage', 'first ion stage', 'second ion stage', 'third ion stage']

plt.figure(0)
#ground-state plot
for i in range(1,5):
	plt.plot(temp,pop[i,:],colors[i-1], label=labellst[i-1])

#-------------------------------------------------------------------------------
for T in np.arange(1,31):
	for r in np.arange(1,5):
		pop[r,T] = sahabolt_E(temp[T],131,r,2) #Change s for energy level


#First excited plot
for i in range(1,5):
	plt.plot(temp,pop[i,:],colors[i-1])
#-------------------------------------------------------------------------------
for T in np.arange(1,31):
	for r in np.arange(1,5):
		pop[r,T] = sahabolt_E(temp[T],131,r,4) #Change s for energy level


#4th plot
for i in range(1,5):
	plt.plot(temp,pop[i,:],colors[i-1])
#-------------------------------------------------------------------------------


plt.yscale('log')
plt.ylim([1e-3,1.1])
plt.legend(loc='lower right')
plt.xlabel('temperature',size=14)
plt.ylabel('population',size=14)
plt.savefig('sahabolt_E_variation')

plt.show()
"""


# ------------------2.7 Saha-Boltzmann populations of hydrogen------------------

def sahabolt_H(temp,elpress,level):

	k = const.k_B
	kerg = k.cgs.value
	keV = k.to(units.eV/units.K).value	
	h = const.h.cgs.value
	keVT = keV * temp
	kergT = kerg * temp
	elmass =  const.m_e.cgs.value 
	eldens = elpress / kergT
	chiion = np.array([7, 16, 31, 51 ])

	for s in range(nrlevels):
		g[0,s] = 2.*(s+1.)**2. # statistical weights
		chiexc[0,s] = 13.598*(1.-1./(s+1.)**2.) # excitation weights
	
	g[1,0] = 1. # statistical weights free proton
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
	ntotal = np.sum(nstage)	# sum both stages = total hydrogen density
	
	# Boltzmann
	nlevel = nstage[0]*g[0,level-1]/u[0]*np.exp(-chiexc[0,level-1]/keVT)
	nlevelrel = nlevel/ntotal # fraction of total hydrogen density
		
	return nlevelrel

print "sahabolt_H(6000,1e2,1):", sahabolt_H(6000,1e2,1)

print "-----------"
for s in range(6):
	print s+1, g[0,s], chiexc[0,s], g[0,s]*np.exp(-chiexc[0,s]/keVT)

print "-----------"
for s in range(0,nrlevels,10):
	print s+1, g[0,s], chiexc[0,s], g[0,s]*np.exp(-chiexc[0,s]/keVT)

"""
sahabolt_H(6000,1e2,1): 0.999999981472
-----------
1 2.0 0.0 2.0
2 8.0 10.1985 4.20197414626e-10
3 18.0 12.0871111111 1.18031923499e-11
4 32.0 12.748125 4.52485028902e-12
5 50.0 13.05408 3.475642845e-12
6 72.0 13.2202777778 3.4031228753e-12
-----------
1 2.0 0.0 2.0
11 242.0 13.4856198347 6.1788471307e-12
21 882.0 13.5671655329 1.86365798838e-11
31 1922.0 13.5838501561 3.90691220224e-11
41 3362.0 13.5899107674 6.73859180545e-11
51 5202.0 13.5927720108 1.03575677222e-10
61 7442.0 13.594345606 1.47635563167e-10
71 10082.0 13.5953025193 1.99564590917e-10
81 13122.0 13.5959274501 2.59362344892e-10
91 16562.0 13.5963579278 3.27028625016e-10

"""

#-------2.8 Solar Ca+ K versus Ha: line strength---------------------

temp = np.arange(1000,20001,100)
CaH = np.zeros(temp.shape)
Caabund = 2.0e-6
for i in range(0,191):
	NCa = sahabolt_E(temp[i],1e2,2,1)# is equal to sahabolt_Ca
	NH = sahabolt_H(temp[i],1e2,2)
	CaH[i] = NCa*Caabund/NH

plt.plot(temp,CaH, label=r'strength ratio Ca$^+$K / H$\alpha$')
plt.yscale('log')
plt.xlabel(r'temperature $T / K$', size=14)
plt.ylabel(r'Ca II K / H$\alpha$', size=14)
plt.legend(fontsize=14)

plt.savefig('strengthratioCa+K')

plt.show()

print 'Ca/H ratio at 5000 K = ', CaH[np.argwhere(temp==5000)][0][0]

"""
#-----------------2.9 Solar Ca+ K versus Ha-temperature sensitivity----

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

plt.figure()
plt.plot(temp,np.absolute(dNHdT), label=r'H')
plt.plot(temp,np.absolute(dNCadT), label=r'Ca$^+$K')

plt.yscale('log')
#plt.ylim(1e-9,1)
plt.xlabel(r'temperature $T/K$', size=14)
plt.ylabel(r"$\left| \left( \Delta n(r,s) / \Delta T \right) /  n(r,s) \right|$", size=20)
plt.legend(loc=4, fontsize=12)
NCa = np.zeros(temp.shape)
NH = np.zeros(temp.shape)

for i in range(101):
	NCa[i] = sahabolt_E(temp[i],1e2,2,1)
	NH[i] = sahabolt_H(temp[i],1e2,2)

plt.plot(temp,NH/np.amax(NH), ls='--',  label = 'rel. pop. H')
plt.plot(temp,NCa/np.amax(NCa), ls='--', label = r'rel. pop. Ca$^+$')
plt.show()

# -------------------- 2.10 Hot stars versus Cool Stars ------------

for T in np.arange(2e3,2e4+1,2e3):
	print T, sahabolt_H(T,1e2,1)

temp = np.arange(1e3,2e4+1,1e2)
nH = np.zeros(temp.shape)
for i in range(191):
	nH[i] = sahabolt_H(temp[i],1e2,1)

plt.plot(temp,nH)
plt.xlabel('temperature $T/K$', size=14)
plt.ylabel('neutral hydrogen fraction', size=14)

"""
