# Many routines to analyse ELGs with J-PLUS datas 

#import sys
#sys.path.append('/home/CEFCA/aaorsi/work/j-plus/')
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gsc

from scipy.interpolate import interp1d
import numpy as np
import jplus
import misc 

import gc

def haversine_dist(ra1, dec1, ra2, dec2):
  th1 = np.pi/2. - dec1 * np.pi/180.0
  th2 = np.pi/2. - dec2 * np.pi/180.0

  ph1 = ra1 * np.pi/180.0
  ph2 = ra2 * np.pi/180.0

  dph = np.abs(ph1 - ph2)
  dth = np.abs(th1 - th2)

  harg = np.sin(dph/2)**2 + np.cos(ph1)*np.cos(ph2) * np.sin(dth/2.)**2

  return 2 *np.arcsin(np.sqrt(harg)) * 180./np.pi  # Return distance in degrees



def continuum_curve(data, tileid, sigma_threshold = 2, LineName = 'J0660', BroadLineName = 'rJAVA', BroadNoLineName='gJAVA',
Plot=False):
  
  import MockJPLUS as mtools
  from scipy import interpolate
  from matplotlib.colors import LogNorm

  igal = jplus.tools.select_object(data, data['tile_id'] == tileid)  # Galaxies in tile 

  dm, fline, fcont = mtools.gen_3fm(igal[LineName][:,0],  igal[BroadLineName][:,0], igal[BroadNoLineName][:,0], 
                           Broad_NoLineName=BroadNoLineName, Broad_LineName=BroadLineName)

  dm_err = mtools.gen_3fm_err(igal[LineName][:,0], igal[LineName][:,1], 
           igal[BroadLineName][:,0], igal[BroadLineName][:,1], 
           igal[BroadNoLineName][:,0], igal[BroadNoLineName][:,1],
           Broad_NoLineName=BroadNoLineName)

#  med_dm = np.median(dm)

  mr = [16,24]
  bsize= 0.1
  toldm = 0.1
  mag_arr = np.arange(mr[0],mr[1],bsize)
  
  nbins = len(mag_arr)
  
  scurve = np.zeros(nbins)

  for i in range(nbins):
    mask = (igal[BroadLineName][:,0] > mag_arr[i] - bsize/2.) & (igal[BroadLineName][:,0] < mag_arr[i] + bsize/2.)
    i_err = dm_err[mask]
    med_dm = np.median(dm[mask])
    nmask = len(i_err)
    scurve[i] = med_dm if nmask == 0 else np.percentile(dm_err[mask], 68.0)*sigma_threshold + med_dm

  
  sfunc = interpolate.interp1d(mag_arr, scurve, fill_value = 'extrapolate')
  #sfunc = interpolate.interp1d(mag_arr, scurve, fill_value = -1, bounds_error=False)
  sout = [mag_arr, scurve]

  idcand = []
  dmcand = []
  dmcand_err = []
  for j in range(len(dm)):
    if (dm[j] >= sfunc(igal[BroadLineName][j,0])) and (dm[j] > 0):
      idcand.append(igal['ids'][j])
      dmcand.append(dm[j])
      dmcand_err.append(dm_err[j])

  if Plot:
    import matplotlib.pyplot as plt

    #Plot.plot(igal[BroadLineName][:,0], dm, 'k,')
    counts, ybins, xbins, image = Plot.hist2d(igal[BroadLineName][:,0],dm,
                                    bins=50, cmap=plt.cm.Greys, normed=LogNorm(), cmin=0.01, 
                                    alpha=0.7,range=[[18,23],[-2,3]])

    
    m2 = np.arange(mr[0]-1, mr[1]+1, bsize)
    Plot.plot(m2, sfunc(m2), '-',color='yellow',linewidth=2)
    
    for j in range(len(dm)):
      if dm[j] >= sfunc(igal[BroadLineName][j,0]):
        plt.plot([igal[BroadLineName][j,0], igal[BroadLineName][j,0]],[dm[j],dm[j]],'.',color='magenta', markersize=2)

    Plot.set_xlim([18.1, 22.9])
    Plot.set_ylim([-0.49,2])
    Plot.text(0.1,0.75,'Tile:%d'%tileid,fontsize=12,transform=Plot.transAxes)
    


  return idcand, dmcand, dmcand_err, sout












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

def make_selection(gal_jplus, by_tile = True, ijlim = 0.6, rjlim = 0.6,snr_limit = 5.0,
ij8lim = -10,cstar_max = 1.0,makeplot = False):

# Create a list of ELG candidates from gal_jplus.
# by_tile: means that this is performed for each tile independently, considering the SNR and errors of each tile.  Otherwise it is 
#          done for all objects at once.
# ijlim, rjlim: limits on i-j0660 and r-j0660 respectively. Default to fiducial values.
# ij8min      : limit on i - j0861, which could help too. Default to -10
# snr_min : minimum SNR of the NB filter.
  
  nbins = 25
  jarr = np.linspace(16,24,nbins)
  jbin = jarr[1] - jarr[0]
  
  snarr = np.zeros(nbins) # signal-to-noise array
  ntotgals = len(gal_jplus['tile_id'])
  kcount = np.zeros(ntotgals, dtype=np.int)
  kk = 0
  tiles_arr = np.unique(gal_jplus['tile_id'])
  ntiles = len(tiles_arr) if by_tile else 1
  print 'Total number of tiles: %d' % ntiles

  npp = 0
  ix = 0
  iy = 0
  for it in range(ntiles):
    if by_tile:
      tile_i = tiles_arr[it]
      seltile = gal_jplus['tile_id'] == tile_i
      idseltile = np.where(seltile)[0]
      gal_tile = jplus.tools.select_object(gal_jplus, seltile)  # selecting objects in tile_i
    else:
      gal_tile = gal_jplus

    sigmarj, sigmaij = make_sel_sigma(gal_tile,jarr) # constructing sigma curves for that tile

    j = gal_tile['J0660'][:,0]
    rj = gal_tile['rJAVA'][:,0] - gal_tile['J0660'][:,0]
    ij = gal_tile['iJAVA'][:,0] - gal_tile['J0660'][:,0]
    ij8 = gal_tile['iJAVA'][:,0] - gal_tile['J0861'][:,0]
    jerr = gal_tile['J0660'][:,1]
    cstar = gal_tile['cstar']
    sn = 1./jerr
    
    for ib in range(nbins):
      isel = np.where((j > jarr[ib]-jbin/2.) & (j < jarr[ib]+jbin/2.))[0]
      snarr[ib] = np.median(sn[isel])
      mag_snr = interp1d(snarr,jarr,bounds_error=False, fill_value=99.0)
      jmin = mag_snr(snr_limit)

    
    icand = np.where(#(rj < rjlim) #& (ij > ijlim)
            (rj > sigmarj(j)) & (ij > sigmaij(j)) &
            (sn > snr_limit) & (ij8 > ij8lim)
            & (cstar < cstar_max))[0]
   
    ncand = len(icand)
    #print 'mag for snr %f: %f. candidates in tile %d: %d' % (snr_limit,jmin,tile_i,ncand)

    if makeplot and npp < 4:
      print ix, iy
      gs = gsc.GridSpec(2,2)
      gs.update(wspace=0.0,hspace=0.0)
      
      ax = plt.subplot(gs[ix,iy])
      ax.plot(gal_tile['J0660'][:,0],gal_tile['rJAVA'][:,0]-gal_tile['J0660'][:,0],',',color='gray')
      ax.plot(gal_tile['J0660'][icand,0],gal_tile['rJAVA'][icand,0]-gal_tile['J0660'][icand,0],
      '.',color='royalblue')
      ax.plot(jarr,sigmarj(jarr),color='red',linewidth=2)

      ax.set_ylim([-0.49,2.99])
      ax.set_xlim([16,21.5])
      ax.text(0.25,.75,'%d' % tile_i,transform=ax.transAxes)
      if iy > 0:
        ax.set_yticklabels([])
        
      iy = iy + 1 if ix == 1 else iy + 0 
      ix = ix + 1 if ix < 1 else 0
      npp += 1
    
    if npp == 4 and makeplot:
      plt.show()
    
    if ncand > 0:
      kcount[kk:kk+ncand] = idseltile[icand]  # From tile it, candidates icand
      kk += ncand
   

  print 'Total Number of candidates: %ld ' % kk 
  gal_cand = jplus.tools.select_object(gal_jplus,kcount[0:kk])
  return gal_cand
 
def zline(lam_line, wf,tf):  # finds the redshift range of a line in a given filter 
  w10 = misc.quantile(wf, tf, 0.001)
  w90 = misc.quantile(wf, tf, 0.999)

  z10 = (w10 - lam_line) / lam_line
  z90 = (w90 - lam_line) / lam_line

  return [z10,z90]
  


def plot_colcol_sdss_jplus(gal_sdss, gal_jplus,zcoord='zspec',xaxis = ['rJAVA','J0660'],yaxis=['iJAVA','J0660'],
add_muse = True, add_composite=True, add_sdsszline = True):

  d2,ind2 = jplus.tools.crossmatch_angular(gal_sdss['coords'],gal_jplus['coords'],max_distance=3e-4)
  m2 = ((d2 != np.inf))
  jplus_iz = jplus.tools.select_object(gal_jplus,m2)
  xax = jplus_iz[xaxis[0]][:,0] - jplus_iz[xaxis[1]][:,0]
  yax = jplus_iz[yaxis[0]][:,0] - jplus_iz[yaxis[1]][:,0]
  zax = gal_sdss['zspec']

  plt.figure
  plt.figure('colcol')

  lab = 'SDSS Spec in J-PLUS'
  plt.plot(xax,yax,'.',color='gray',label = lab)
  
  plt.ylim([-0.5,2.5])
  plt.xlim([-0.5,2.5])

  plt.xlabel('%s-%s' % (xaxis[0],xaxis[1]),fontsize=20)
  plt.ylabel('%s-%s' % (yaxis[0],yaxis[1]),fontsize=20)

  # Now add composite spectra

  spec_elg = jplus.datasets.fetch_eboss_elg_composite()
  
  
  zarr = [0.005,0.8, 0.356, 0.30]
  w0   = [6560., 3727., 5007., 4860.0,4959]
  wname = [r'H\alpha',r'[OII]',r'[OIII]_{\lambda 5007}',r'H\beta', r'[OIII]_{\lambda 4959}']
  nline = len(w0)
  new = 50
  ewarr = np.logspace(0.5,3.0,num=new)
  color = plt.cm.coolwarm(np.linspace(0.,1,new))

  j = 1

  icount = 0
  
  if add_composite:
    icount += 1
    idd = np.where( np.concatenate([xaxis,yaxis]) == 'J0660')[0]
    # Add composite spectrum only if plot includes the J0660 filter.
    if len(idd) > 0:
      ztrack = np.zeros([new,2])
      for iew in range(new):    
          i = ewarr[iew]
          specline = jplus.datasets.add_line_to_spec(spec_elg,EW=i,obsframe=True,wavelength0=w0[0])
          specc = {'w':specline['w'],'flux':specline['flux']}
          conv_mags = jplus.datasets.compute_jplus_photometry_singlespec(specc)
          
          xax = conv_mags[xaxis[0]] - conv_mags[xaxis[1]]
          yax = conv_mags[yaxis[0]] - conv_mags[yaxis[1]]


          ztrack[iew,0] = xax
          ztrack[iew,1] = yax
          
  #        plt.plot([xax,xax],[yax,yax],'.',color=color[iew])
      
      sc = plt.scatter(ztrack[:,0],ztrack[:,1],s = 25,c=ewarr,cmap=plt.cm.coolwarm)
      cb = plt.colorbar(sc)
      cb.ax.set_title('EW')

  ijlim = 0.6
  rjlim = 0.0

#  plt.plot([ijlim,3],[ijlim,ijlim],linewidth=4,color='black')        
#  plt.plot([rjlim,rjlim],[ijlim,3],linewidth=4,color='black')        
#  plt.legend()       

  j0660 = jplus.datasets.fetch_jplus_filter('J0660')
  color = plt.cm.Paired(np.linspace(0.1,.9,nline))
  
  if add_sdsszline:
    icount += 10
    for il in range(nline):
      zr = zline(w0[il],j0660.wave,j0660.throughput)
      iz_sdss = np.where((gal_sdss['zspec'] > zr[0]) & 
                (gal_sdss['zspec'] < zr[1]))
                
      d2,ind2 = jplus.tools.crossmatch_angular(gal_sdss['coords'][iz_sdss],gal_jplus['coords'],max_distance=3e-4)
      m2 = ((d2 != np.inf))
      jplus_iz = jplus.tools.select_object(gal_jplus,m2)
  #    r_j = jplus_iz['rJAVA'][:,0] - jplus_iz['J0660'][:,0]
  #    i_j = jplus_iz['iJAVA'][:,0] - jplus_iz['J0660'][:,0]
      
      xax = jplus_iz[xaxis[0]][:,0] - jplus_iz[xaxis[1]][:,0]
      yax = jplus_iz[yaxis[0]][:,0] - jplus_iz[yaxis[1]][:,0]

      plt.plot(xax,yax,'o',color=color[il],markersize=8,label=r'$%s, %.2f<z<%.2f$' % (wname[il],zr[0],zr[1]) )


  if add_muse: # also add MUSE spectra in colour-colour plot
    icount += 100 
    for il in range(nline):
    
      zr = zline(w0[il],j0660.wave,j0660.throughput) # show OII emitters
      speclist = get_musewide_spec(zr,strong='O2')
     
      nz = len(speclist)

      for iz in range(nz):
        ww = speclist[iz]['w']
        ff = speclist[iz]['flux'][:,0]
        nw = len(ww)
        for ii in range(nw):
          ff[ii] = 0.0 if ff[ii] < 0 else ff[ii]  # remove negative fluxes (???)

        
        flux = np.zeros([nw,2]) # recreating flux structure for compatibility
        flux[:,0] = ff
        specc = {'w':ww,'flux':flux}
        conv_mags = jplus.datasets.compute_jplus_photometry_singlespec(specc)
        
        #r_j = conv_mags['rJAVA'] - conv_mags['J0660']
        #i_j = conv_mags['iJAVA'] - conv_mags['J0660']
        #g_r = conv_mags['gJAVA'] - conv_mags['rJAVA']
    

        xax = conv_mags[xaxis[0]] - conv_mags[xaxis[1]]
        yax = conv_mags[yaxis[0]] - conv_mags[yaxis[1]]

        plt.plot([xax,xax],[yax,yax],'D',color=color[il],
        markersize=10,label='MUSE Wide-Survey $%.2f <z<%.2f$' % (zr[0],zr[1]) if iz == 0 else '')


  plt.legend(loc = 'upper left')
  plt.savefig('%s-%s.%s-%s_%d.pdf' % (xaxis[0],xaxis[1],yaxis[0],yaxis[1],icount),bbox_inches='tight')
  plt.close()
  return ijlim,rjlim


def find_xmatches(gal_orig, gal_match,zcond = None,zcoord = 'zspec',getz=True,max_distance=3e-4):
# Returns a cross-match between 2 galaxy catalogues
# zcond can receive a 2-element list specifying a redshift range to apply to gal_orig.
# zcoord is the name of the redshift coordinate in gal_orig.

  if zcond is not None:
    iz = ((gal_orig[zcoord] > zcond[0]) & (gal_orig[zcoord] < zcond[1]))
    galcat = jplus.tools.select_object(gal_orig,iz)
  else:
    galcat = gal_orig

  d2,ind2 = jplus.tools.crossmatch_angular(galcat['coords'],gal_match['coords'],max_distance=max_distance)
  m2 = ((d2 != np.inf))
  xmatch = jplus.tools.select_object(galcat,m2)
  
  if getz:
    match_clean = jplus.tools.select_object(gal_match, ind2[m2])
    xmatch[zcoord] = match_clean[zcoord]

  return xmatch



def get_musewide_spec(zrange,strong=None,name=''):
# downloads muse 1d spectra of objects at a given z-range and returns the spectra
# strong = 'O2', 'O3', etc... select objects with that primary line. None selects them all

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

  if strong is not None:
    lead_line = data['LEAD_LINE'] == strong
  else:
    lead_line = 1

  iz = np.where((data['Z'] > zrange[0])
       & (data['Z'] < zrange[1]) & lead_line)[0]

  nz = len(iz)
  print 'Number of MUSE spectra in %.2f < z < %.2f: %d' % (zrange[0],zrange[1],nz)

  speclist = []

  for i in range(nz):
    finput = 'spectrum_%s.fits.gz' % data['UNIQUE_ID'][iz[i]]
    url = '%s%s' % (url_root,finput)
   
    spdirout = '../elg_jplus/spec/muse/'
    sfile = '%s%s' % (spdirout,finput)
    sfile = sfile[:-3]
    if not os.path.isfile(sfile):
      print 'file %s not found, downloading...' % (sfile)
#      import pdb ; pdb.set_trace()
      filename = wget.download(url)
      print subprocess.check_output(['mv',filename,spdirout])
      filename = '../elg_jplus/spec/muse/%s' % filename
      print subprocess.check_output(['gunzip','-f',filename])
    else:
      filename = finput

    fitsfile = filename[:-3]
    dd = fits.open(fitsfile)
    ww = dd[1].data['WAVE_VAC']
    ff = dd[1].data['FLUX']

    nw = len(ww)
    flux = np.zeros([nw,2])
    flux[:,0] = ff * 1e-20


    speclist.append({'flux':flux,'w':ww,'z':data['Z'][iz[i]],'survey':'MUSE-Wide',
    'file':sfile,'name':name})
    dd.close()

    
  return speclist



def get_eboss_spec(zrange,typegal="'GALAXY'",name=''):
# Returns a list of spectra from eBOSS galaxies already downloaded
#
  from astropy.io import fits
  datadir = '/home/CEFCA/aaorsi/work/elg_jplus/spec/eboss/'
  
  listfile = datadir + 'optical_search_27334.csv'
  ldata = np.loadtxt(listfile,delimiter=',',dtype='S')
  nd = len(ldata)
  
  zlist = np.zeros(nd)
  idspec = []
  typespec = []
   


  for ii in range(nd):
    zlist[ii] = np.float(ldata[ii][7])
    idspec.append('%s-%s-%04d' % (ldata[ii][0],ldata[ii][1],np.int(ldata[ii][2])))
    typespec.append(ldata[ii][10])

  
  sel = np.where((np.array(typespec) == typegal) & (zlist >= zrange[0]) & (zlist <= zrange[1]))[0]
  
  ns = len(sel)

  if ns == 0:
    print 'WARNING: No eBOSS spectra found for %s between %f < z < %f' % (typegal, zrange[0],zrange[1])
    return []

  flux = np.zeros([ns,2])
  ww   = np.zeros(ns)

  specout = []

  print 'Retrieving eBOSS spectra...'
  for jj in range(ns):
    idj = sel[jj]
    ffile = '%sspec-%s.fits' % (datadir,idspec[jj])
    hdu = fits.open(ffile)

    ff    = 1e-17 * hdu[1].data['model'] # using model fluxes
    ww    = 10**(hdu[1].data['loglam'])
    nw    = len(ww)
    flux  = np.zeros([nw,2])
    flux[:,0] = ff

    hdu.close()
    specout.append({'flux':flux,'w':ww,'z':zlist[idj],'survey':'eBOSS',
    'file':ffile,'name':name})
    
  print '%d spectra recorded.' % ns
  return specout


def get_vvds_spec(zrange,survey = 'All',name=''):
# Returns a list of spectra from all VVDS surveys.

  from astropy.io import fits
  import glob

  datadir = '/home/CEFCA/aaorsi/work/elg_jplus/spec/vvds/'
  
  if survey == 'All':
    surveyarr = ['F10_WIDE','F14_WIDE','F22_WIDE','F02_DEEP','CDFS_DEEP','F02_UDEEP']
  else:
    surveyarr = survey

  nsurvey = len(surveyarr)
  specout = []

  for s_i in range(nsurvey):  # loops over each survey
    fdata = '%scesam_vvds_sp%s.txt' % (datadir, surveyarr[s_i])
    zcol = 5 if surveyarr[s_i] == 'F02_UDEEP' else 4
    specid, zarr = np.loadtxt(fdata,unpack=True,usecols=(0,zcol))
    #zarr = dd[4]
    #specid = dd[0]

    zsel = np.where((zarr >= zrange[0]) & (zarr <= zrange[1]))[0]
    nz = len(zsel)

    if nz == 0:
      print '0 galaxies in %f<z<%f found for %s' % (zrange[0],zrange[1],surveyarr[s_i])
      continue
    
    print '%d objects in %s' % (nz, surveyarr[s_i])

    for jj in range(nz):
      specfile = glob.glob('%s%s/1D/sc_%d*atm_clean*' % (datadir,surveyarr[s_i],specid[zsel[jj]]))
      if len(specfile) == 1:
        h = fits.open(specfile[0])
        fl = h[0].data
        ww   = h[0].header['CRVAL1'] + h[0].header['CDELT1'] * np.arange(2,h[0].header['NAXIS1']+2)
        nww = len(ww)
        flux = np.zeros([nww,2])
        flux[:,0] = fl

        specout.append({'flux': flux, 'w':ww,'z':zarr[zsel[jj]],'survey':'VVDS-'+surveyarr[s_i],
        'file':specfile,'name':name})
        h.close()



  return specout

  
def get_vipers(zrange,name='VIPERS'):
# Returns a list of spectra from all VVDS surveys.

  from astropy.io import fits
  import glob

  datadir = '/home/CEFCA/aaorsi/work/elg_jplus/spec/vipers/'
  
  surveyarr = ['W1', 'W4']

  nsurvey = len(surveyarr)
  specout = []

  for s_i in range(nsurvey):  # loops over each survey
    fdata = '%s/%s/VIPERS_%s_SPECTRO_PDR2.txt' % (datadir, surveyarr[s_i], surveyarr[s_i])
    zcol = 9
    specid, zarr = np.loadtxt(fdata,unpack=True,usecols=(2,zcol))
    zsel = np.where((zarr >= zrange[0]) & (zarr <= zrange[1]))[0]
    nz = len(zsel)

    if nz == 0:
      print '0 galaxies in %f<z<%f found for %s' % (zrange[0],zrange[1],surveyarr[s_i])
      continue
    
    print '%d objects in VIPERS %s' % (nz, surveyarr[s_i])

    for jj in range(nz):
      specfile = '%s%s/VIPERS_%d.fits' % (datadir,surveyarr[s_i],specid[zsel[jj]])
      with fits.open(specfile,memmap=False) as h:
        fl = h[1].data['fluxes']
        ww =  h[1].data['waves']
        nww = len(ww)
        flux = np.zeros([nww,2])
        flux[:,0] = fl
        specout.append({'flux': flux, 'w':ww,'z':zarr[zsel[jj]],'survey':surveyarr[s_i],
        'file':specfile,'name':name})
        del fl
        del ww
        del h


  return specout




def get_elg_photoz(gal_elgs,noJ0660 = False,overwrite = True, filterflag = None, emlines=True, suffix='whatever',
                  plot_photoz=False, gallibname='GISSEL'):
#2) Run LePhare on JPLUS Broad-band photometry, recalibrating tile by tile
  
  matplotlib.rcParams['figure.figsize'] = (8,6)
  #d,ind = jplus.tools.crossmatch_angular(gal_elgs['coords'],gal_eboss['coords'],max_distance=3e-4)
  #m2 = ((d != np.inf))
  #elg_eboss = jplus.tools.select_object(gal_elgs,m2)

  for ifilter in jplus.datasets.jplus_filter_names():
    if ifilter in gal_elgs:
      ind = np.isfinite(gal_elgs[ifilter][:,0])
      gal_elgs[ifilter][~ind,:] = [-99,-99]
      gal_elgs[ifilter][gal_elgs[ifilter][:,0]==99,:] = [-99,-99]
      gal_elgs[ifilter][gal_elgs[ifilter][:,0]==0,:] = [-99,-99]
      
      gal_elgs[ifilter][gal_elgs[ifilter][:,1]==0,1] = 0.1
      

  #for ifilter in jplus.datasets.jplus_filter_names():
  #    print ifilter
      #print elg_eboss[ifilter] 
  if filterflag == 'alljplus' or filterflag is None:
    fflag = [1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0]
  elif filterflag == 'bbjplus':
    fflag =  [1,1,1,1,1,1,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0]
  elif filterflag == 'noJ0660':
    fflag = [1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0]
  elif filterflag == 'SDSS_JPLUSBB':
    fflag = [0,1,1,1,1,0,1,0,1,0,1,0,1,1,1,1,1,0,0,0,0]
  else:
    print 'filterflag not recognised. Assuming alljplus'
    fflag = [1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0]

    
  #elg_eboss['redshift'] = np.zeros(len(elg_eboss['rJAVA'][:,0]))
  #gal_elgs['redshift'] = np.zeros(len(gal_elgs['rJAVA'][:,0]))

  if noJ0660:

    Lephare_noha = jplus.photoz.LePhare(gal_elgs, per_tile=False, outspec=False, recalibration=False,
                                 filterflag=bbfilters,
                                 suffix='_'+suffix,emlines=True,filename='noJ0660')

    Lephare_noha.prepare(overwrite=True)
    jp_photoz_noha = Lephare_noha.run(overwrite=True)

  Lephare_all = jplus.photoz.LePhare(gal_elgs, per_tile=False, outspec=False, recalibration=False,
                                 filterflag=fflag,
                                 suffix='_'+suffix,emlines=emlines,filename='allfilters', gallibname=gallibname)

  Lephare_all.prepare(overwrite=overwrite)
  jp_photoz_all = Lephare_all.run(overwrite=overwrite)

  if plot_photoz:
    if noJ0660:
      plt.hist(jp_photoz_noha['photoz'],range=[.01,.9],bins=50,alpha = 0.5,color='blue',label='No J0660')

    plt.hist(jp_photoz_all['photoz'],range=[.01,.9],bins=50,alpha = 0.5,color='red',label = 'all filters photo-z')

    plt.xlabel('redshift')

    plt.legend()
    plt.savefig('photoz.pdf',bbox_inches='tight')

  
  return jp_photoz_all

