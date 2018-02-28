import numpy as np
import math as math
import sys
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.colors as colors
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0 = 73.0, Om0 = 0.25)
import deepdish as dd
from scipy.interpolate import interp1d
import os.path


c = 2.99792458 * 10 **18 #Angst/s

Lines_Added = ['Lyalpha', 'Hbeta', 'Halpha','OII_3727','OII_3729','OIII_5007','OIII_4959',\
               'OI_6300','NII_6548','NII_6584','SII_6717','SII_6731','NeIII_3870']
Lines_lambda_0 = np.array([1216.0, 4861.0, 6563.0, 3727.0, 3729.0, 5007.0, 4959.0\
                    ,6300.0, 6548.0, 6584.0, 6717.0, 6731.0, 3870.0]) # I remove the CII_158um, NII_205um lines

plotdir = '/global/users/dizquierdo/david/Code_python/plots/'
#Lines_emission = ['Halpha', 'OIII+Hbeta','OII']

def redshif_lines(x):
        ''' Devuelve el redshift interval of the line inside f660'''
        red = {
        'Ha': [0.,0.0268170044187],
        'OI': [0.0298412698413,0.0696825396825],
        'Hb': [0.334704793252, 0.386340259206],
        'OIII5007': [0.29578589974, 0.345915717995],
        'OIII4959': [0.308328291994, 0.35894333535],
        'OII3727': [0.739876642532, 0.808156694392]
        }.get(x,'Error')
        if red == 'Error':
                print 'Da fuck?!? Line not found', x
                sys.exit()
        return red



def import_JPLUS_filters(filter2import):
        name = 'jp' + str(filter2import+1) + '_trans.dat'
        direc = '/home/CEFCA/aaorsi/work/elg_jplus/jdata/T80Cam_JAST_TransmissionCurvesTmp_20160518//' + name
        wl, trans, dummy1, dummy2, dummy3 = np.loadtxt(direc, unpack = True, skiprows = 1, comments = '#')
        return wl , trans


def lambdaPivot(nFilters):
        lPivot=np.zeros(nFilters)
        for i in np.arange(0,nFilters,1):
                longi, trans = import_JPLUS_filters(i)
                LS=np.zeros(len(longi))
                S_L=np.zeros(len(longi))
                for j in range(len(LS)):
                        LS[j] = (trans[j]*longi[j])
                        S_L[j]= (trans[j]/longi[j])
                lPivot[i]=(np.trapz(LS)/np.trapz(S_L))**0.5
        return lPivot


def alpha(Filter):
        Long, Trans = import_JPLUS_filters(Filter)
        Long = np.array(Long)
        Trans = np.array(Trans)
        suma = 0
        norm = 0
        for i in range(len(Long)):
                suma = suma + (Long[i]**2)*Trans[i]
                norm = norm + (Long[i] * Trans[i])
        alpha = suma/norm

        return alpha


def beta(Filter,longi):
        Long, Trans = import_JPLUS_filters(Filter)

        Long = np.array(Long)
        Trans = np.array(Trans)
        norm = 0
        for i in range(len(Long)):
                norm = norm + (Long[i] * Trans[i])
        aux = np.where(Long > longi)
        pos = aux[0][0]
        num = longi * Trans[pos+1]
        beta = num/norm
        return beta


def norm(Filter):
        longi, trans = import_JPLUS_filters(Filter)
        return np.trapz(trans/longi, x = longi, dx = longi[1] - longi[0])


def name_to_num_function():
        names_filters = ['uJAVA' , 'J0378', 'J0395', 'J0410', 'J0430', 'gSDSS',
                         'J0515', 'rSDSS', 'J0660', 'iSDSS', 'J0861', 'zSDSS']
        name_to_num = {}
        for nombre in names_filters:
                name_to_num[nombre] = names_filters.index(nombre)

        return name_to_num


def Filters_centers():
        names_filters = ['uJAVA' , 'J0378', 'J0395', 'J0410', 'J0430', 'gSDSS',
                         'J0515', 'rSDSS', 'J0660', 'iSDSS', 'J0861', 'zSDSS']
        fil_cent = [3485 , 3785, 3950, 4100, 4300, 4803,
                         5150, 6254, 6600, 7668, 8610, 9114]
        name_to_center = {}
        for nombre in names_filters:
                name_to_center[nombre] = fil_cent[names_filters.index(nombre)]
        return name_to_center



def app_mag(data, filname):
# Compute apparent magnitudes from lightcone input
# data: lightcone data
# filname: string with filter name

  name_to_num = name_to_num_function()
  nn_line_filter = name_to_num[filname]
  distance = ((data['pos'][:,0]**2 + data['pos'][:,1]**2 + data['pos'][:,2]**2)**0.5) / cosmo.h
  dl = (1+data['redshift'])*(distance)

  
  mag = data['ObsMagDust'][:,nn_line_filter]\
    + 5*(np.log10(dl[:])+5) # apparent 
  return mag


name_to_num = name_to_num_function()
centers = Filters_centers()
nFilters = 12

def tfm(ggC, line_filter = 'J0660', broad_withline= 'rSDSS', broad_noline='iSDSS'):
  
  lcentre = centers[line_filter]
  nn_line_filter = name_to_num[line_filter]
  nn_broad_withline = name_to_num[broad_withline]
  nn_broad_noline = name_to_num[broad_noline]
  lPiv = lambdaPivot(nFilters)
  alpha660 = alpha(nn_line_filter)
  alphai = alpha(nn_broad_noline)
  alphaR=  alpha(nn_broad_withline)
  beta660 = beta(nn_line_filter,lcentre)
  betaR = beta(nn_broad_withline,lcentre)

  mag_line            = app_mag(ggC, line_filter)
  mag_broad_withline  = app_mag(ggC, broad_withline)
  mag_broad_noline    = app_mag(ggC, broad_noline)

  Fr   = (10**(-0.4*(mag_broad_withline + 48.6)))*c/(lPiv[nn_broad_withline]**2) # erg/s cm2 A
  Fi   = (10**(-0.4*(mag_broad_noline + 48.6)))*c/(lPiv[nn_broad_noline]**2) # erg/s cm2 A
  F660 = (10**(-0.4*(mag_line + 48.6)))*c/(lPiv[nn_line_filter]**2) # erg/s cm2 A

  Fline = ((Fr-Fi)-((alphaR-alphai)/(alpha660-alphai)*(F660-Fi)))/((beta660*((alphai-alphaR)/(alpha660-alphai)))+betaR)
  M = (F660 - Fi - (beta660*Fline))/(alpha660-alphai)
  N = Fi - alphai*((F660-Fi-(beta660*Fline))/(alpha660-alphai))
  Fcont = (M * lcentre) + N # erg / s cm2 A
  Fcont = (lPiv[nn_line_filter]**2/c) * Fcont  # erg / s cm2 Hz
  
  mag_Cont_3FM = -2.5*np.log10(Fcont) - 48.6
  
  dm3FM = mag_Cont_3FM - mag_line

  return dm3FM

def gen_3fm(linemag, broad_line, broad_noline, LineFilterName='J0660', Broad_LineName='rSDSS',Broad_NoLineName='iSDSS'):
  lPiv      = lambdaPivot(nFilters)
  lcentre   = centers[LineFilterName]
  alpha660  = alpha(name_to_num[LineFilterName])
  alphai    = alpha(name_to_num[Broad_NoLineName])
  alphaR    = alpha(name_to_num[Broad_LineName])
  beta660   = beta(name_to_num[LineFilterName],lcentre)
  betaR     = beta(name_to_num[Broad_LineName],lcentre)

  Fr   = (10**(-0.4*(broad_line + 48.6)))*c/(lPiv[name_to_num[Broad_LineName]]**2) # erg/s cm2 A
  Fi   = (10**(-0.4*(broad_noline + 48.6)))*c/(lPiv[name_to_num[Broad_NoLineName]]**2) # erg/s cm2 A
  F660 = (10**(-0.4*(linemag + 48.6)))*c/(lPiv[name_to_num[LineFilterName]]**2) # erg/s cm2 A

  Fline = ((Fr-Fi)-((alphaR-alphai)/(alpha660-alphai)*(F660-Fi)))/((beta660*((alphai-alphaR)/(alpha660-alphai)))+betaR)
  M = (F660 - Fi - (beta660*Fline))/(alpha660-alphai)
  N = Fi - alphai*((F660-Fi-(beta660*Fline))/(alpha660-alphai))
  Fcont = (M * lcentre) + N # erg / s cm2 A
  Fcont *= (lPiv[name_to_num[LineFilterName]]**2/c) # erg / s cm2 Hz
  
  mag_Cont_3FM = -2.5*np.log10(Fcont) - 48.6
  
  dm3FM = mag_Cont_3FM - linemag

  return dm3FM


def Three_FM(ggC,ggL,ancho_sin_contaminar = 'iSDSS'):
  ########################################################## 
  # This function given a galaxy with lines ggL and a galaxy
  # without lines and with an auxiliar filter retrieves the
  # estimated cont and the true one
  ########################################################## 

  name_to_num = name_to_num_function()
  centers = Filters_centers()
  estrecho = 'J0660'
  ancho_contaminado = 'rSDSS'
  lcentre = centers[estrecho]
  nn_estrecho = name_to_num[estrecho]
  nn_ancho_contaminado = name_to_num[ancho_contaminado]
  nn_ancho_sin_contaminar = name_to_num[ancho_sin_contaminar]
  nFilters = 12
  lPiv = lambdaPivot(nFilters)
  c = 2.99792458 * 10 **18 #Angst/s
  alpha660 = alpha(nn_estrecho)
  alphai = alpha(nn_ancho_sin_contaminar)
  alphaR=  alpha(nn_ancho_contaminado)
  beta660 = beta(nn_estrecho,lcentre)
  betaR = beta(nn_ancho_contaminado,lcentre)
  distance = ((ggL['pos'][:,0]**2 + ggL['pos'][:,1]**2 + ggL['pos'][:,2]**2)**0.5) / cosmo.h
  dl = (1+ggL['redshift'])*(distance)

  ggC['ObsMagDust'][:,nn_estrecho] = ggC['ObsMagDust'][:,nn_estrecho]\
    + 5*(np.log10(dl[:])+5) # apparent 
  ggL['ObsMagDust'][:,nn_estrecho] = ggL['ObsMagDust'][:,nn_estrecho]\
       + 5*(np.log10(dl[:])+5) # apparent
  ggL['ObsMagDust'][:,nn_ancho_contaminado] = ggL['ObsMagDust'][:,nn_ancho_contaminado]\
         + 5*(np.log10(dl[:])+5) # apparent 
  ggL['ObsMagDust'][:,nn_ancho_sin_contaminar] = ggL['ObsMagDust'][:,nn_ancho_sin_contaminar]\
         + 5*(np.log10(dl[:])+5) # apparent 

  Fr = (10**(-0.4*(ggL['ObsMagDust'][:,nn_ancho_contaminado] + 48.6)))\
        *c/(lPiv[nn_ancho_contaminado]**2) # erg/s cm2 A
  Fi = (10**(-0.4*(ggL['ObsMagDust'][:,nn_ancho_sin_contaminar] + 48.6)))\
        *c/(lPiv[nn_ancho_sin_contaminar]**2) # erg/s cm2 A
  F660 = (10**(-0.4*(ggL['ObsMagDust'][:,nn_estrecho] + 48.6)))\
        *c/(lPiv[nn_estrecho]**2) # erg/s cm2 A
  Fline = ((Fr-Fi)-((alphaR-alphai)/(alpha660-alphai)*(F660-Fi)))/((beta660*((alphai-alphaR)/(alpha660-alphai)))+betaR)
  M = (F660 - Fi - (beta660*Fline))/(alpha660-alphai)
  N = Fi - alphai*((F660-Fi-(beta660*Fline))/(alpha660-alphai))
  Fcont = (M * lcentre) + N # erg / s cm2 A
  Fcont = (lPiv[nn_estrecho]**2/c) * Fcont  # erg / s cm2 Hz
  mag_Line = ggL['ObsMagDust'][:,nn_estrecho]
  mag_Cont_3FM = -2.5*np.log10(Fcont) - 48.6
  mag_Cont_SAM = ggC['ObsMagDust'][:,nn_estrecho]
  
  dm3FM = mag_Cont_3FM - mag_Line
  dmSAM = mag_Cont_SAM - mag_Line

  return dm3FM, dmSAM
