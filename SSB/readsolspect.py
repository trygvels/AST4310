import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.ticker as mtick
# reading solspect.dat
wav, F,Fcont,I,Icont= np.loadtxt("solspect.dat",usecols=(0,1,2,3,4), unpack=True)



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
c=3e14
factor = wav**2/c*1e10
print 'max(Ic) = ',np.max(Icont*factor),'at',wav[np.where(Icont*factor == np.max(Icont*factor))]
#Plot
fig = plt.figure()
plt.plot(wav,F*factor, color=tableau20[5])
plt.plot(wav,Fcont*factor, color=tableau20[4])
plt.plot(wav,I*factor, color=tableau20[7])
plt.plot(wav,Icont*factor, color=tableau20[6])

#Changing fonts
plt.rc('text', usetex=True)
plt.rc('font', **{'family':'sans-serif','sans-serif':['Helvetica']})
plt.legend([r"$F_{cont}^{\nu}$",r"$F_{cont}^{\nu}$ smoothed",r"$I_{cont}^{\nu}$",r"$F_{cont}^{\nu}$ smoothed"], loc=1)
#plt.xlim(-500,2.5e3)
#plt.ylim(0,10000)
plt.xlabel(r"wavelength $\lambda$ [$\mu m$]")
plt.ylabel(r"intensity, flux [$erg cm^{-2} s^{-1} ster^{-1} Hz^{-1}$]")
plt.title(r"solar continuum radiation (solspect.dat)")

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
ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%g'))
# Remove the tick marks; they are unnecessary with the tick lines we just plotted.    
plt.tick_params(axis="both", which="both", bottom="off", top="off",    
                labelbottom="on", left="off", right="off", labelleft="on")  
#-------------------------------------------------------------

# if you want/need to save the plot in some format, you can use
# (bbox and pad make the figure to be tighten to the plot-box)
fig.savefig('spectraldistFreq.pdf', bbox_inches='tight',pad_inches=0.106)
plt.show()
