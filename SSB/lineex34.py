import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import special
h, tau5, colm, temp, vturb, nhyd, nprot, nel, ptot, pgasptot, dens = np.loadtxt("falc.dat",usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)
sigma, obs = np.loadtxt("int_nad.dat",usecols=(0,2), unpack=True)
wav1 = 1./sigma*1e8 #Convert to um
wav1 = 0.99972683*wav1+0.0107-196.25/wav1 #Converting to air wavelength

chiion = [5.139, 47.29, 71.64]  #Sodium ionization energies [eV]
# s=1 : ground state, angular momentum l=0
# s=2 : Na D1 upper level l=1
# s=3 : Na D2 upper level l=1
e = 4.803204e-10 #electron charge [statcoulomb]
hp = 6.62607e-27 # Planck constant (erg s)
c = 2.99792e10 # light speed [cm/s]
k=1.380658e-16	# Boltzmann constant [erg/K]
def partfunc_Na(temp):
	# partition functions Na
	# input: temp (K)
	# output: float array(3) = partition functions U1,U2,U3
	u=np.zeros(3)
	# partition function Na I: follow Appendix D of Gray 1992
	# log(U1(T)) = c0 + c1*log(theta) + c2*log(theta)^2 +
	#              c3*log(theta)^3 + c4 log(theta)^4
	# with theta=5040./T
	theta=5040./temp
	# partition function Na I : Appendix D of Gray (1992)
	c0=0.30955
	c1=-0.17778
	c2=1.10594
	c3=-2.42847
	c4=1.70721
	logU1 = (c0 + c1*np.log10(theta) + c2*np.log10(theta)**2 +c3*np.log10(theta)**3 + c4*np.log10(theta)**4)
	u[0]=10**logU1

	# partition function Na II and Na III: approximate by the
	# statistical weights of the ion ground states
	u[1]=1 # from Allen 1976
	u[2]=6 # from Allen 1976
	return u

def boltz_Na(temp, r, s):
	keV = 8.61734e-5
	u = partfunc_Na(temp)
	erg2eV=1/1.60219e-12 # erg to eV conversion
	E_n=np.zeros(3) # energy level: E_n[0]=0 : ground state
	E_n[1]=hp*c/5895.94e-8*erg2eV # Na D1: 2.10285 eV
	E_n[2]=hp*c/5889.97e-8*erg2eV # Na D2: 2.10498 eV

	g = [2.,2.,4.]
	relnrs = g[s-1] / u[r-1]*np.exp(-(E_n[s-1])/ (keV* temp))
	return relnrs

#THIS GIVES PLOT IN INSTRUCTIONS
distribution = np.zeros(len(temp))

distribution1 = np.zeros(len(temp))

distribution2 = np.zeros(len(temp))
for i in range(len(temp)):
	distribution[i] = boltz_Na(temp[i],1,1)
	distribution1[i] = boltz_Na(temp[i],1,2)
	distribution2[i] = boltz_Na(temp[i],1,3)
plt.plot(h, distribution)
plt.plot(h, distribution1)
plt.plot(h, distribution2)
plt.legend(['s=1 ground state', 's=2 Na D1', 's=3 Na D2'])
plt.xlim([-100,2000])
plt.show()

def saha_Na(temp, eldens, r):
	kerg = 1.380658e-16
	keV = 8.61734e-5
	elmass = 9.10939e-28
	hp = 6.62607e-27 #[ergs]
	keVT = keV * temp
	kergT = kerg * temp	
	chiion = [5.139, 47.29, 71.64] 
	u = partfunc_Na(temp)
	u = np.append(u, 2) 	# With this command we are adding a new element to the array
	sahaconst = (2. * np.pi * elmass * kergT / (hp**2))**1.5 * 2. / eldens
	nstage = np.zeros(4)
	nstage[0] = 1. 	# We set the first element of the array to a value 1
	for i in range(len(chiion)):
		nstage[i + 1] = nstage[i] * sahaconst * u[i + 1] / u[i] * np.exp(-chiion[i] / keVT)
	ntotal = np.sum(nstage)
	nstagerel = nstage / ntotal
	return nstagerel[r - 1]

#THIS SHOWS THE RIGHT PLOT
distribution = np.zeros((len(temp),2))
for j in range(2):
	for i in range(len(temp)):
		distribution[i,j] =  saha_Na(temp[i],nel[i],j+1)

plt.semilogy(h, distribution[:,0])
plt.semilogy(h, distribution[:,1],'--')
plt.xlim([-100,2000])
plt.ylim([1e-4,1e1])
plt.legend(['s=1','s=2'])
plt.show()


def sahabolt_Na(temp, eldens, r, s):
	return saha_Na(temp, eldens, r) * boltz_Na(temp, r, s)
"""
distribution = np.zeros(len(temp))
for i in range(len(temp)):
	distribution[i] = sahabolt_Na(temp[i],nel[i],1,2)
plt.plot(h,distribution)
plt.show()
"""
def voigt(A,x):
	z = (x+1j*A)
	V = special.wofz(z).real
	return V



def gammavdw_NaD(temp, pgas, s):
	# Van der Waals broadening for Na D1 and Na D2
	# s=2 : Na D1
	# s=3 : Na D2
	# using classical recipe by Unsold
	# following recipe in SSB
	rsq_u = rsq_NaD(s)
	rsq_l = rsq_NaD(1)
	# lower level D1 and D2 lines is ground state s=1
	loggvdw=6.33 + 0.4*np.log10(rsq_u - rsq_l)+ np.log10(pgas) - 0.7*np.log10(temp)
	return 10**loggvdw

def rsq_NaD(s):
	# compute mean square radius of level s of Na D1 and Na D2 transitions
	#  -> needed for van der Waals broadening in SSB
	# s=1 : ground state, angular momentum l=0
	# s=2 : Na D1 upper level l=1
	# s=3 : Na D2 upper level l=1
	hp=6.62607e-27 # Planck constant (erg s)
	c=2.99792e10 # light speed [cm/s]
	erg2eV=1/1.60219e-12 # erg to eV conversion
	E_ionization = 5.139 # [eV] ionization energy
	E_n=np.zeros(3) # energy level: E_n[0]=0 : ground state
	E_n[1]=hp*c/5895.94e-8*erg2eV # Na D1: 2.10285 eV
	E_n[2]=hp*c/5889.97e-8*erg2eV # Na D2: 2.10498 eV
	Z=1.  # ionization stage, neutral Na: Na I
	Rydberg=13.6 # [eV] Rydberg constant
	l=[0.,1.,1.] # angular quantum number
	nstar_sq = Rydberg*Z**2 / (E_ionization - E_n[s-1])
	rsq=nstar_sq / 2. / Z**2*(5*nstar_sq + 1 - 3*l[s-1]*(l[s-1] + 1))
	return rsq


def planck(temp,wav):
	h = 6.626076e-27 # Planck constant [erg s]
	blambda = 2*h*c**2/(wav**5*(np.exp(h*c/(wav*k*temp))-1)) #Planck function

	return blambda


def exthmin(wav,temp,eldens):

	# H-minus extinction, from Gray 1992
	# input:
	#  wav = wavelength [Angstrom] (float or float array)
	#  temp = temperature [K]
	#  eldens = electron density [electrons cm-3]
	# output:
	#  H-minus bf+ff extinction [cm^2 per neutral hydrogen atom]
	#  assuming LTE ionization H/H-min
	# physics constants in cgs (all cm)
	
	h=6.626076e-27	# Planck constant [erg s]
	
	# other parameters
	theta=5040./temp
	elpress=eldens*k*temp

	# evaluate H-min bound-free per H-min ion ? Gray (8.11)
	# his alpha = my sigma in NGSB/AFYC (per particle without stimulated)
	sigmabf = (1.99654 -1.18267E-5*wav +2.64243E-6*wav**2-4.40524E-10*wav**3 +3.23992E-14*wav**4-1.39568E-18*wav**5 +2.78701E-23*wav**6)
	sigmabf *= 1e-18 # cm^2 per H-min ion
	if np.size(wav) > 1:
		sigmabf[np.where(wav>16444)] = 0 # H-min ionization limit at lambda = 1.6444 micron
	elif ( np.size(wav) == 1):
		if wav> 16444:
			sigmabf = 0

	# convert into bound-free per neutral H atom assuming Saha = Gray p135
	# units: cm2 per neutral H atom in whatever level (whole stage)
	graysaha=4.158E-10*elpress*theta**2.5*10.**(0.754*theta)# Gray (8.12)
	kappabf=sigmabf*graysaha # per neutral H atom
	kappabf=kappabf*(1.-np.exp(-h*c/(wav*1E-8*k*temp)))# correct stimulated

	# check Gray's Saha-Boltzmann with AFYC (edition 1999) p168
	# logratio=-0.1761-np.log10(elpress)+np.log10(2.)+2.5*np.log10(temp)-theta*0.754
	# print 'Hmin/H ratio=',1/(10.**logratio) # OK, same as Gray factor SB

	# evaluate H-min free-free including stimulated emission = Gray p136
	lwav = np.log10(wav)
	f0 = - 2.2763 - 1.6850*lwav + 0.76661*lwav**2 - 0.0533464*lwav**3
	f1 =   15.2827 - 9.2846*lwav + 1.99381*lwav**2 - 0.142631*lwav**3
	f2 = - 197.789 + 190.266*lwav - 67.9775*lwav**2 + 10.6913*lwav**3 - 0.625151*lwav**4
	ltheta = np.log10(theta)
	kappaff = 1e-26*elpress*10**(f0+f1*ltheta+f2*ltheta**2) #Gray(8.13)
	return kappabf+kappaff

def dopplerwidth(wav,temp,v_t):
	'''Takes in central wavelength in cm,temperature in K, v_t in cm/s, 
	and m in grams and returns dopplerwidth in cm.'''	
	m = 22.99*1.6606e-24 #g
	kerg = 1.380658e-16
	return wav/c*np.sqrt(2.*kerg*temp/m + v_t*v_t)

def NaD1_ext(wav,temp,eldens,nhyd,vmicro,pgas):
	echarge=4.803204e-10          # electron charge [statcoulomb] (ESU: electrostatic unit, cgs)
	kerg = 1.380658e-16	
	s = 2
	r = 1
	AN = 1.8e-6 #Na abundance
	nna = nhyd*AN	
	flu = 0.318 #Oscillator strength for NaI D 1 (0.631 for NaI D 2)
	me = 9.1093897e-28 #g
	nlte = sahabolt_Na(temp,eldens,r,1) #Ground state
	gamma = gammavdw_NaD(temp,pgas,s)

	lambda0 = 5895.94e-8	#lambda0 [cm] is center wavelength of Na D lines
	lambdaD = dopplerwidth(wav,temp,vmicro)
	#lambdaD = lambda0/c*np.sqrt(2*k*temp/m+vmicro**2)

	v_voigt = (wav-lambda0)/lambdaD
	a_voigt = wav**2*gamma/(4*np.pi*c*lambdaD)
	voigt_NaD = voigt(a_voigt, v_voigt) / lambdaD
	
	#print "gamma",gamma
	#print "SahaBoltz:", nlte
	#print "wav:", wav
	#print "voigt_NaD", voigt_NaD

	#voigt_NaD REALLY big!
	NaD1extinction = np.sqrt(np.pi)*echarge**2/(me*c)*wav**2/c*nlte*nna*flu*voigt_NaD*(1-np.exp(-hp*c/(wav*kerg*temp)))
	return NaD1extinction

sigma_Thomson= 6.648e-25
pgas = ptot*pgasptot#Calculating gas pressure

#----ADOPTED WAVELENGTH SPECTER----------
#nw=1000
#offs=np.linspace(-2e-8,2e-8,num=nw)  # brute force with very dense wavelength grid
# bit more intelligent: dense sampling in core, sparser in wings and even sparser far wings
l0=np.linspace(.01,.25,num=25)
l1=np.linspace(l0[-1]+.01,1,num=10)
l2=np.linspace(l1[-1]+.1,2,num=10)
offs=np.concatenate((np.flipud(-l2),np.flipud(-l1),np.flipud(-l0),[0],l0,l1,l2))*1e-8

lambda0 = 5895.94e-8	
wav=offs+lambda0
nw=len(offs)
int_calc = np.zeros(nw)
nh=len(h)
for w in range(nw):
    wl=wav[w]
    ext=np.zeros(nh)
    tau=np.zeros(nh)
    integrand=np.zeros(nh)
    intt=0.
    for i in range(nh):
        cext = exthmin(wl*1e8, temp[i], nel[i])*(nhyd[i]-nprot[i]) + sigma_Thomson*nel[i]
        lext = NaD1_ext(wl, temp[i], nel[i], nhyd[i], vturb[i]*1e5, pgas[i])
        ext[i]  = cext + lext
        tau[i] = tau[i-1] + 0.5 * (ext[i] + ext[i-1]) * (h[i-1]-h[i])*1E5
        integrand[i] = planck(temp[i],wl)*np.exp(-tau[i])
        intt += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
    int_calc[w]=intt



# ---------------- New color scheme -----------------
# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.) 
#-------------------------------------------------------

fig = plt.figure()
#--------------------Configuration------------------------------
plt.grid("on")
# Remove the plot frame lines. They are unnecessary chartjunk.    
ax = plt.subplot(111)    
ax.spines["top"].set_visible(False)    
ax.spines["bottom"].set_visible(False)    
ax.spines["right"].set_visible(False)    
ax.spines["left"].set_visible(False)   
plt.xticks(fontsize=14)    
plt.yticks(fontsize=14)    
plt.grid(b=True, which='major')
plt.grid(b=True, which='minor', alpha=0.2)
plt.minorticks_on()
# Remove the tick marks; they are unnecessary with the tick lines we just plotted.    
plt.tick_params(axis="both", which="both", bottom="off", top="off",    
                labelbottom="on", left="off", right="off", labelleft="on")  


plt.plot(wav*1e8, int_calc/int_calc[0],color=tableau20[4])
plt.plot(wav1[2500:3500], obs[2500:3500]/max(obs),'--',color=tableau20[2])
plt.xlim(5894,5+5893)
plt.title("Na I D1 in LTE in FALC")
plt.xlabel("wavelength")
plt.ylabel('Normalized intensity')
#-------------------------------------------------------------
# if you want/need to save the plot in some format, you can use
# (bbox and pad make the figure to be tighten to the plot-box)
fig.savefig('extinct.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()

