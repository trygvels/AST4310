import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
h, tau5, colm, temp, vturb, nhyd, nprot, nel, ptot, pgasptot, dens = np.loadtxt("falc.dat",usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)
wav, F,Fcont,I,Icont= np.loadtxt("solspect.dat",usecols=(0,1,2,3,4), unpack=True)

print wav

def planck(temp,wav):
	k = 1.38065e-16
	h = 6.626076e-27 # Planck constant [erg s]
	c = 2.997929e10  # Velocity of light [cm/s]

	blambda = 2*h*c**2/(wav**5*(np.exp(h*c/(wav*k*temp))-1))

	return blambda
def brightTemp(wav,I):
	k = 1.38065e-16
	h = 6.626076e-27 # Planck constant [erg s]
	c = 2.997929e10  # Velocity of light [micro m/s]

	Tb = h*c/(wav*k*np.log(2*h*c**2/(wav**5*I)+1))
	return Tb
		

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
	k=1.380658e-16	# Boltzmann constant [erg/K]
	h=6.626076e-27	# Planck constant [erg s]
	c=2.997929e10	# velocity of light [cm/s]
	
	# other parameters
	theta=5040./temp
	elpress=eldens*k*temp

	# evaluate H-min bound-free per H-min ion ? Gray (8.11)
	# his alpha = my sigma in NGSB/AFYC (per particle without stimulated)
	sigmabf = (1.99654 -1.18267E-5*wav +2.64243E-6*wav**2-4.40524E-10*wav**3 			  +3.23992E-14*wav**4-1.39568E-18*wav**5 +2.78701E-23*wav**6)
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
fig = plt.figure()

sigma_Thomson= 6.648e-25
# SSB 2.4 page 16 emergent intensity, contribution function and mean height of formation in FALC sigma_Thomson = 6.648E-25 # Thomson cross-section [cm^2]
mu = np.linspace(0.1,1.,10)
print wav
contfuncMu = np.zeros((len(mu),len(wav)))
contfuncCalc=np.zeros(len(wav))
for j in range(len(wav)):
	
	for w in range(len(mu)):
		ext = np.zeros(len(tau5))
		tau = np.zeros(len(tau5))
		integrand = np.zeros(len(tau5))
		contfunc = np.zeros(len(tau5))
		intt = 0.0
		hint = 0.0
		for i in range(1, len(tau5)):
			ext[i] = (exthmin(wav[j]*1e4, temp[i], nel[i])*(nhyd[i]-nprot[i])+ sigma_Thomson*nel[i])
			tau[i] = tau[i-1] + 0.5*(ext[i] + ext[i-1])*(h[i-1]-h[i])*1e5
			integrand[i] = planck(temp[i],wav[j]*1e-4)*np.exp(-tau[i]/mu[w])
			intt += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])/mu[w]
			hint += h[i]*0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])/mu[w]
			contfunc[i] = integrand[i]*ext[i]
		contfuncMu[w,j]=intt
		# note : exthmin has wavelength in [Angstrom], planck in [cm]
		#mean = hint / intt
		#tau5[69]=1
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

sintheta = np.sqrt(1-mu**2)
wavs=[0.2,0.3,0.4,0.6,0.7,0.9,1.,2.,5.]
data = np.zeros((10,10))
#for j in range(10):
#	for i in range(10):
#		data[i,j] = contfuncMu[:,np.where(wav==wavs[j-1])][i]

#for j in range(10):
#	plt.plot(sintheta,data[:,j]/data[0,0],color=tableau20[j])

plt.grid(b=True, which='major')
plt.grid(b=True, which='minor', alpha=0.2)
plt.minorticks_on()

plt.xlabel(r"$sin(\theta)=r/R")
plt.ylabel(r"Intensity [$10^{14} erg s^{-1} cm^{-2} ster^{-1} cm^{-1}$]")
plt.title('Limb darkening in FALC')
plt.legend([r'$\mu = 0.1$', r'$\mu = 0.2$', r'$\mu = 0.3$', r'$\mu = 0.4$', r'$\mu = 0.5$', r'$\mu = 0.6$', r'$\mu = 0.7$', r'$\mu = 0.8$', r'$\mu = 0.9$', r'$\mu = 1.0$'])

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
# Remove the tick marks; they are unnecessary with the tick lines we just plotted.    
plt.tick_params(axis="both", which="both", bottom="off", top="off",    
                labelbottom="on", left="off", right="off", labelleft="on")  
#-------------------------------------------------------------
# if you want/need to save the plot in some format, you can use
# (bbox and pad make the figure to be tighten to the plot-box)
fig.savefig('limbdarkRatio.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()








