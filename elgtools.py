# Many routines to analyse ELGs with J-PLUS datas 

import sys
sys.path.append('/home/CEFCA/aaorsi/work/j-plus/')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gsc

from scipy.interpolate import interp1d
import numpy as np
import jplus
import misc 


def make_sel_sigma(gal_jplus,jarr):
# construct a sigma curve to get rid of objects dominated by colour terms consistent with noise.

  nbins = len(jarr)
  jbin = jarr[1] - jarr[0]
  rj = gal_jplus['rJAVA'][:,0] - gal_jplus['J0660'][:,0]
  ij = gal_jplus['iJAVA'][:,0] - gal_jplus['J0660'][:,0]
  ij8 = gal_jplus['iJAVA'][:,0] - gal_jplus['J0861'][:,0]

  j = gal_jplus['J0660'][:,0]
  rj_err = np.sqrt(gal_jplus['rJAVA'][:,1]**2 + gal_jplus['J0660'][:,1]**2)
  ij_err = np.sqrt(gal_jplus['iJAVA'][:,1]**2 + gal_jplus['J0660'][:,1]**2)
  mean_rj = np.mean(rj)
  mean_ij = np.mean(ij)
  tol = 0.1
  sigma2_col_rj = np.zeros(nbins)
  sigma2_col_ij = np.zeros(nbins)
  for i in range(nbins):
      isel = np.where((j > jarr[i]-jbin/2.) & (j < jarr[i]+jbin/2.))[0]
      if len(isel) > 0:
          errcol_rj = rj_err[isel]
          sigma2_col_rj[i] = np.percentile(errcol_rj,99)

          errcol_ij = ij_err[isel]
          sigma2_col_ij[i] = np.percentile(errcol_ij,99)

          
  sigmarj = interp1d(jarr,sigma2_col_rj,fill_value = 99,bounds_error=False) # mags outside magarr are discarded
  sigmaij = interp1d(jarr,sigma2_col_rj,fill_value = 99,bounds_error=False)

  return sigmarj, sigmaij

def make_selection(gal_jplus, by_tile = True, ijlim = 0.6, rjlim = 0.6,snr_limit = 5.0,ij8lim = -10):

# Create a list of ELG candidates from gal_jplus.
# by_tile: means that this is performed for each tile independently, considering the SNR and errors of each tile.  Otherwise it is 
#          done for all objects at once.
# ijlim, rjlim: limits on i-j0660 and r-j0660 respectively. Default to fiducial values.
# ij8min      : limit on i - j0861, which could help too. Default to -10
# snr_min : minimum SNR of the NB filter.
  
  nbins = 20
  jarr = np.linspace(16,24,nbins)
  jbin = jarr[1] - jarr[0]
  
  snarr = np.zeros(nbins) # signal-to-noise array
  ntotgals = len(gal_jplus['tile_id'])
  kcount = np.zeros(ntotgals, dtype=np.int)
  kk = 0
  tiles_arr = np.unique(gal_jplus['tile_id'])
  ntiles = len(tiles_arr) if by_tile else 1
  print 'Total number of tiles: %d' % ntiles
  
  for it in range(ntiles):
    if by_tile:
      tile_i = tiles_arr[it]
      seltile = gal_jplus['tile_id'] == tile_i
      gal_tile = jplus.tools.select_object(gal_jplus, seltile)  # selecting objects in tile_i
    else:
      gal_tile = gal_jplus

    sigmarj, sigmaij = make_sel_sigma(gal_tile,jarr) # constructing sigma curves for that tile

    j = gal_tile['J0660'][:,0]
    rj = gal_tile['rJAVA'][:,0] - gal_tile['J0660'][:,0]
    ij = gal_tile['iJAVA'][:,0] - gal_tile['J0660'][:,0]
    ij8 = gal_tile['iJAVA'][:,0] - gal_tile['J0861'][:,0]
    jerr = gal_tile['J0660'][:,1]
    sn = 1./jerr
    
    for ib in range(nbins):
      isel = np.where((j > jarr[ib]-jbin/2.) & (j < jarr[ib]+jbin/2.))[0]
      snarr[ib] = np.median(sn[isel])
      mag_snr = interp1d(snarr,jarr,bounds_error=False, fill_value=99.0)
      jmin = mag_snr(snr_limit)

    
    icand = np.where((rj > rjlim) & (ij > ijlim)
            & (rj > sigmarj(j)) & (ij > sigmaij(j)) &
            (sn > snr_limit) & (ij8 > ij8lim))[0]
   
    ncand = len(icand)
    print 'mag for snr %f: %f. candidates in tile %d: %d' % (snr_limit,jmin,tile_i,ncand)

    if ncand > 0:
      kcount[kk:kk+ncand] = seltile[icand]  # From tile it, candidates icand
    
    kk += ncand

  gal_cand = jplus.tools.select_object(gal_jplus,kcount[0:kk])
  return gal_cand
 
def zline(lam_line, wf,tf):  # finds the redshift range of a line in a given filter 
  w10 = misc.quantile(wf, tf, 0.02)
  w90 = misc.quantile(wf, tf, 0.98)

  z10 = (w10 - lam_line) / lam_line
  z90 = (w90 - lam_line) / lam_line

  return [z10,z90]
  


def plot_colcol_sdss_jplus(gal_sdss, gal_jplus,zcoord='zspec',add_muse=True):

  d2,ind2 = jplus.tools.crossmatch_angular(gal_sdss['coords'],gal_jplus['coords'],max_distance=3e-4)
  m2 = ((d2 != np.inf))
  jplus_iz = jplus.tools.select_object(gal_jplus,m2)
  r_j = jplus_iz['rJAVA'][:,0] - jplus_iz['J0660'][:,0]
  i_j = jplus_iz['iJAVA'][:,0] - jplus_iz['J0660'][:,0]
  zax = gal_sdss['zspec']

  plt.figure('colcol')

  lab = 'SDSS Spec in J-PLUS'
  plt.plot(r_j,i_j,'.',color='gray',label = lab)
  
  plt.ylim([-0.5,2.5])
  plt.xlim([-0.5,2.5])

  plt.xlabel('rJAVA - J0660',fontsize=20)
  plt.ylabel('iJAVA - J0660',fontsize=20)

  # Now add composite spectra

  spec_elg = jplus.datasets.fetch_eboss_elg_composite()
  
  
  zarr = [0.005,0.8, 0.356, 0.30]
  w0   = [6560., 3727., 5007., 4860.0]
  wname = [r'H\alpha',r'[OII]',r'[OIII]',r'H\beta']
  nline = len(w0)
  new = 10
  ewarr = np.logspace(0.5,3.0,num=new)
  color = 'blue'

  j = 1

  ztrack = np.zeros([new,2])
  for iew in range(new):    
      i = ewarr[iew]
      specline = jplus.datasets.add_line_to_spec(spec_elg,EW=i,obsframe=True,wavelength0=w0[0])
      specc = {'w':specline['w'],'flux':specline['flux']}
      conv_mags = jplus.datasets.compute_jplus_photometry_singlespec(specc)
      r_j = conv_mags['rJAVA'] - conv_mags['J0660']
      i_j = conv_mags['iJAVA'] - conv_mags['J0660']
      g_r = conv_mags['gJAVA'] - conv_mags['rJAVA']
      ztrack[iew,0] = r_j
      ztrack[iew,1] = i_j
      
  plt.plot(ztrack[:,0],ztrack[:,1],color=color)
              

  ijlim = 0.6
  rjlim = 0.0

  plt.plot([ijlim,3],[ijlim,ijlim],linewidth=4,color='black')        
  plt.plot([rjlim,rjlim],[ijlim,3],linewidth=4,color='black')        
  plt.legend()       

  j0660 = jplus.datasets.fetch_jplus_filter('J0660')
  color = plt.cm.Paired(np.linspace(0.1,.9,nline))
  
  for il in range(nline):
    zr = zline(w0[il],j0660.wave,j0660.throughput)
    iz_sdss = np.where((gal_sdss['zspec'] > zr[0]) & 
              (gal_sdss['zspec'] < zr[1]))
              
    d2,ind2 = jplus.tools.crossmatch_angular(gal_sdss['coords'][iz_sdss],gal_jplus['coords'],max_distance=3e-4)
    m2 = ((d2 != np.inf))
    jplus_iz = jplus.tools.select_object(gal_jplus,m2)
    r_j = jplus_iz['rJAVA'][:,0] - jplus_iz['J0660'][:,0]
    i_j = jplus_iz['iJAVA'][:,0] - jplus_iz['J0660'][:,0]
    
    plt.plot(r_j,i_j,'o',color=color[il],markersize=12,label=r'$%s, %.2f<z<%.2f$' % (wname[il],zr[0],zr[1]) )


  if add_muse: # also add MUSE spectra in colour-colour plot
    
    for il in range(nline):
    
      zr = zline(w0[il],j0660.wave,j0660.throughput) # show OII emitters
      nz, speclist = get_musewide_spec(zr)
     
      for iz in range(nz):
        ww = speclist[iz]['WAVE_VAC']
        ff = speclist[iz]['FLUX']
        ff += np.abs(ff.min()) # avoid negative flux
        nw = len(ww)
        
        flux = np.zeros([nw,2]) # recreating flux structure for compatibility
        flux[:,0] = ff
        specc = {'w':ww,'flux':flux}
        conv_mags = jplus.datasets.compute_jplus_photometry_singlespec(specc)
        
        r_j = conv_mags['rJAVA'] - conv_mags['J0660']
        i_j = conv_mags['iJAVA'] - conv_mags['J0660']
        g_r = conv_mags['gJAVA'] - conv_mags['rJAVA']
    
        plt.plot([r_j,r_j],[i_j,i_j],'D',color=color[il],
        markersize=14,label='MUSE Wide-Survey %.2f <z<%.2f' % (zr[0],zr[1]) if iz == 0 else '')






  plt.legend()
  plt.show()
  return ijlim,rjlim


def find_xmatches(gal_orig, gal_match,zcond = None,zcoord = 'zspec'):
# Returns a cross-match between 2 galaxy catalogues
# zcond can receive a 2-element list specifying a redshift range to apply to gal_orig.
# zcoord is the name of the redshift coordinate in gal_orig.

  if zcond is not None:
    iz = ((gal_orig[zcoord] > zcond[0]) & (gal_orig[zcoord] < zcond[1]))
    galcat = jplus.tools.select_object(gal_orig,iz)
  else:
    galcat = gal_orig

  d2,ind2 = jplus.tools.crossmatch_angular(galcat['coords'],gal_match['coords'],max_distance=3e-4)
  m2 = ((d2 != np.inf))
  xmatch = jplus.tools.select_object(galcat,m2)

  return xmatch



def get_musewide_spec(zrange):
# downloads muse 1d spectra of objects at a given z-range and returns the spectra

  from astropy.io import fits
  import wget
  import subprocess
  import os.path

  musedir = '/home/CEFCA/aaorsi/work/MUSE-Wide_Survey/'
  url_root= 'http://data.muse-vlt.eu/MW_1-24/1d_spectra/'
  
  maintable = 'MW_1-24_main_table.fits'

  dd = fits.open(musedir + maintable)
  data = dd[1].data
  dd.close()

  iz = np.where((data['Z'] > zrange[0])
       & (data['Z'] < zrange[1]))[0]

  nz = len(iz)
  print 'Number of MUSE spectra in %.2f < z < %.2f: %d' % (zrange[0],zrange[1],nz)

  speclist = []

  for i in range(nz):
    finput = 'spectrum_%s.fits.gz' % data['UNIQUE_ID'][iz[i]]
    url = '%s%s' % (url_root,finput)
    
    if not os.path.isfile(finput):
      filename = wget.download(url)
      print subprocess.check_output(['gunzip','-f',filename])
    else:
      filename = finput

    fitsfile = filename[:-3]
    dd = fits.open(fitsfile)
    speclist.append(dd[1].data)
    dd.close()

    
  return nz, speclist






