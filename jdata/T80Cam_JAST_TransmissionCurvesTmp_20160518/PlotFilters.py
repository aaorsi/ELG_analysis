import numpy as np
from astropy.io import ascii
import pylab as pl
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0 = 73.0, Om0 = 0.25)
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

###############################################################################
#                                FUNCTIONS                                    #
###############################################################################
def PlotFiltros():
    tabla=ascii.read('r_band_trans.dat')
    rSloan=tabla['L']
    transmisionr=tabla['t']
    
    tabla=ascii.read('i_band_trans.dat')
    iSloan=tabla['L']
    transmisioni=tabla['t']
    
    tabla=ascii.read('z_band_trans.dat')
    zSloan=tabla['L']
    transmisionz=tabla['t']
    
    tabla=ascii.read('g_band_trans.dat')
    gSloan=tabla['L']
    transmisiong=tabla['t']
    
    plt.plot(rSloan,transmisionr,color='black',linewidth=0.7,linestyle='--',label='$\mathrm{r_{SDSS}}$')
    plt.plot(iSloan,transmisioni,color='blue',linewidth=0.7,linestyle='--',label='$\mathrm{i_{SDSS}}$')
    plt.plot(zSloan,transmisionz,color='green',linewidth=0.7,linestyle='--',label='$\mathrm{z_{SDSS}}$')
    plt.plot(gSloan,transmisiong,color='red',linewidth=0.7,linestyle='--',label='$\mathrm{g_{SDSS}}$')
    plt.title("JPlus filters",fontsize=20)
    for i in range(12):
        longOnda=[]
        transmision=[]
        c=['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'gray', 'lightcoral', 'maroon', 'lightpink', 'lime', 'sienna']
        na=['F348', 'F378', 'F395','F410', 'F430', '$\mathrm{g_{sdss}}$', 'F515', '$\mathrm{r_{sdss}}$', 'F660', '$\mathrm{i_{sdss}}$', 'F861', '$\mathrm{z_{sdss}}$']
        centre=[3485,3785,3950,4100,4300,4750,5150,6250,6600,7725,8610,9150]
        nfiltro=i+1
        name1='jp'
        name2='_trans.dat'
        namefile=name1 + str(nfiltro) + name2
        file1 = open(namefile,'r')
        next(file1)
        for line in file1:
            aux= line.split()
            longOnda.append(aux[0])
            transmision.append(aux[1])
        
        plt.plot(longOnda,transmision,linewidth=1.5, color=c[i],label=na[i])
        n=str(centre[i])
        #plt.plot(centre[i], 0.01 ,'ro')
        #plt.text(centre[i], 0.87 , n, fontsize=8)
        #plt.axvline(x=centre[i], ymin=0, ymax=1, c="black",linewidth=0.5)
        #plt.plot(centre[i],0,'ro',label='n')
        plt.legend(loc="upper right",prop={'size':12.5})
        plt.xlabel('$\mathrm{\lambda}[\AA]$',fontsize=18)
        plt.ylabel('Transmision',fontsize=18)
        
###############################################################################
#                               END FUNCTIONS                                 #
###############################################################################


PlotFiltros()




