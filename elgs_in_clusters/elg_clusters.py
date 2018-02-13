import os

jplus_dir         =  '/home/CEFCA/aaorsi/work/j-plus/'
elg_analysis_dir  = '/home/CEFCA/aaorsi/work/elg_jplus/'
elgdata           = '%s/out/elgs.dat' % elg_analysis_dir
redmapper_dir     = '/home/CEFCA/aaorsi/work/redmapper/'
redmapperdata     = redmapper_dir + 'redmapper_dr8_public_v6.3_catalog.fits'
tilesdata         = '%s/tiles/tiles_data.tsv' % elg_analysis_dir


cwd = os.getcwd()
import sys

sys.path.append(jplus_dir)
sys.path.append(elg_analysis_dir)

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits

import jplus
#import elg_analysis as elg
import elgtools as tools_elg
import learn_elgs as learn_elg
import pickle

from astropy import units as u
from astropy.coordinates import SkyCoord

tile_scale = 1.40

rmin = -2
rmax =np.log10(tile_scale/2.)
nbins = 20

rarr = np.linspace(rmin, rmax, nbins)
dr = rarr[1] - rarr[0]

t_info = np.loadtxt(tilesdata)
print 'tiles info read'

def get_distance(ra1, dec1, ra2, dec2):
  ra1 *= np.pi/180.
  ra2 *= np.pi/180.
  dec1 *= np.pi/180.
  dec2 *= np.pi/180.
  
# Angular distance for two sets of coordinates
  cosg = (np.cos(np.pi/2. - dec1)*np.cos(np.pi/2. - dec2) +
          np.sin(np.pi/2. - dec1)*np.sin(np.pi/2.-dec2)*np.cos(ra1 - ra2))
  
  return np.arccos(cosg)


def quick_dist(ra1, dec1, ra2, dec2, units='deg'):
  return np.sqrt( ((ra1 - ra2)*np.cos(dec1*np.pi/180))**2 + (dec1 - dec2)**2)

def haversine_dist(ra1, dec1, ra2, dec2):
  th1 = np.pi/2. - dec1 * np.pi/180.0
  th2 = np.pi/2. - dec2 * np.pi/180.0

  ph1 = ra1 * np.pi/180.0
  ph2 = ra2 * np.pi/180.0

  dph = np.abs(ph1 - ph2)
  dth = np.abs(th1 - th2)

  harg = np.sin(dph/2)**2 + np.cos(ph1)*np.cos(ph2) * np.sin(dth/2.)**2

  return 2 *np.arcsin(np.sqrt(harg)) * 180./np.pi  # Return distance in degrees


def skydist(ra1,dec1,ra2,dec2, units='deg'):
  
  n2 = len(ra2)

  if n2 != len(dec2):
    raise ValueError('ra, dec arrays have different number of elements')
  
  dist = np.zeros(n2)

  c1 = SkyCoord(ra1,dec1,unit=units)
  for i in range(n2):
    c2 = SkyCoord(ra2[i], dec2[i], unit=units)
    dist[i] = c1.separation(c2).radian

  return dist




def load_catalogues(elgfile = elgdata, centralobj = redmapperdata, centralobjtype = 'fits', 
                    zrange=[.3,.35]):


# Function to load elgs and redmapper catalogues

  elgs = pickle.load(open(elgfile))
  print 'ELG catalogue loaded'

  if centralobjtype == 'fits':
    raw_central = fits.open(centralobj)
    central     = raw_central[1].data
    zsel = np.where((central['z_lambda'] > zrange[0]) & (central['z_lambda'] < zrange[1]))[0]
    ncentral = len(zsel)
  else:
    central = centralobj
    ncentral = len(central['ra'])
    zsel     = np.arange(ncentral)

  
  return {'elgs':elgs, 'central':central,'zsel':zsel}

  
def get_density(datadic):  
  
  gal_arr = np.zeros(nbins)
  central = datadic['central']
  zsel    = datadic['zsel']
  elgs    = datadic['elgs']
  ncentral  = len(zsel)

  #dfunc = quick_dist  # Define which distance calculator to use
  dfunc = haversine_dist # get_distance  # Define which distance calculator to use
  
  print 'central objects loaded, %d clusters' % ncentral

  # Select only those central objects 
  # near a J-PLUS tile
  
  
  elgtiles = np.unique(elgs['tile_id'])
  ntiles = len(elgtiles)
  print 'ELG sample is contained in %d tiles' % ntiles
  print 'scanning tiles..'
 
  ktot = 0
  for i in range(ncentral):
    idc = zsel[i]

    dist = dfunc(central['ra'][idc], central['dec'][idc], t_info[:,1], 
                t_info[:,2])
    idist = np.argmin(dist)
    mindist = dist[idist]
    
    if mindist < tile_scale:
      idtt = np.where(elgtiles == t_info[idist,0])[0]
      
      if len(idtt) != 1:
        continue

      sel_tile = np.where(elgs['tile_id'] == t_info[idist,0])[0]
      ngals = len(sel_tile)
  
      ktot += 1

      for k in range(ngals):
        idg = sel_tile[k]
        dist = np.log10(dfunc(central['ra'][idc], central['dec'][idc], 
        elgs['coords'][idg,0], elgs['coords'][idg,1]))

        if dist < rmax:
          dbin = int((dist  - rmin)/dr)
          gal_arr[dbin:] += 1

#    if ktot == 79:
#      break

  print gal_arr/ktot

  print 'N centers: ',ktot
  return gal_arr/ktot


def plot_densities(d1,d2):
  import matplotlib.pyplot as plt
  vol_r = np.pi*(10**rarr)**2
  d1v = d1/vol_r
  d2v = d2/vol_r

  plt.semilogy(rarr,d1v,color='blue',label='ELG density around clusters',linewidth=3)
  plt.semilogy(rarr,d2v,color='red',label='ELG density around random')

  plt.legend(fontsize=20)

  plt.xlabel(r'$\log(\theta [{\rm deg}])$',fontsize=20)
  plt.ylabel(r'$\langle n \rangle (\leq r) [{\rm deg}^{-2}]$',fontsize=20)

  plt.xlim([rarr.min(),rarr.max()])
  plt.ylim([np.min([d1v,d2v]), d2v.max()])

  plt.show()

  return 1


#def elg_clusters():
d1 = load_catalogues()
tcenter = {'ra':t_info[:,1], 'dec':t_info[:,2]}
d2 = load_catalogues(centralobj=tcenter, centralobjtype='')

import pymangle

tiles = jplus.datasets.fetch_jplus_tile_list(upad=False,overwrite=True)
stars = jplus.datasets.fetch_jplus_objects_to_mask(overwrite=True)
mask  = jplus.mangle.Mask(tiles,stars, overwrite=True)
ran = mask.create_random(3e4)

elgs = d1['elgs']
good = mask.contains(elgs['coords'])

d3 = {}
d3['elgs'] = ran['coords']


import ipdb ; ipdb.set_trace()

g1 = get_density(d1)
g2 = get_density(d2)

plot_densities(g1,g2)

  
