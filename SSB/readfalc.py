import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc



# reading falc.dat
h, tau5, colm, temp, vturb, nhyd, nprot, nel, ptot, pgasptot, dens = np.loadtxt("falc.dat",usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)

# ALL CONSTANTS HERE
m_h = 1.67352e-24
m_he = 3.97*m_h
m_e = 9.10939e-28
m_p = 1.67262e-24
n_he = 0.1*nhyd
pgas = ptot-dens*vturb**2/2
k = 1.38065e-16
n_phot = 20*temp**3
# plotting

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
#Plot
fig = plt.figure()
print "Photon density at deepest model location: %g" % n_phot[np.argmin(h)]
print "Hydrogen density at deepest model location: %g" % nhyd[np.argmin(h)]
print "Photon density at highest model location: %g" % (20*5770**3/(2*np.pi))
print "Hydrogen density at highest model location: %g" % nhyd[np.argmax(h)]

plt.semilogy(h, nhyd, color=tableau20[4])
plt.semilogy(h,n_phot, color=tableau20[6])

#Changing fonts
plt.rc('text', usetex=True)
plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
plt.legend([r"$n_{H}$",r"$n_{Photon}$"], loc=3)
#plt.xlim(-500,2.5e3)
#plt.ylim(0,10000)
plt.xlabel(r"heigth h [km]")
plt.ylabel(r"Number density [$cm^3$]")
plt.title(r"Photon number density (FALC)")
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
fig.savefig('photondens.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()
