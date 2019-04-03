import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

wav, F,Fcont,I,Icont= np.loadtxt("solspect.dat",usecols=(0,1,2,3,4), unpack=True)
#Constants in cgs
k = 1.38065e-16
h = 6.626076e-27 # Planck constant [erg s]
c = 2.997929e10  # Velocity of light [micro m/s]

def planck(temp,wav):
	blambda = 2*h*c**2/(wav**5*(np.exp(h*c/(wav*k*temp))-1))

	return blambda
def brightTemp(wav,I):
	var1= h*c/(wav*k)
	var2= 2*h*c**2/(I*wav**5)*1e-4
	return var1/np.log(var2+1)
		
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

#--------------- Plotting -----------------------
fig = plt.figure()
plt.plot(wav, brightTemp(wav*1e-4,Icont*1e10), color=tableau20[4])

#Changing fonts
plt.rc('text', usetex=True)
plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})

#plt.legend([r'$B_{\lambda}(T=6300)$',r'$B_{\lambda}(T=5777)$'],loc=1)

#plt.xlim(0,2)
#plt.ylim(0,10000)
plt.xlabel(r"wavelength $\lambda$ [$\mu m$]")
plt.ylabel(r"$T_b$ [K]")
plt.title(r"Brightness temperature (solspect.dat)")



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
fig.savefig('brightTemp.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()
