
jplusdir = '/home/CEFCA/aaorsi/work/j-plus/'

import sys
sys.path.append(jplusdir)
import os
os.chdir(jplusdir)

import jplus
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gsc

import elgtools as elg

jplus.plotting.setup_text_plots(fontsize=10,usetex=True)
spec_elg = jplus.datasets.fetch_eboss_elg_composite()
matplotlib.rcParams['figure.figsize'] = (12,10)


#### Options #####

AddedPlots      = False  # Plot added stuff

LoadCatalogues  = True  # Should be true most of the time
PlotComposite   = False
ZpLephare       = True  # Recalibrate zero points using Lephare
FindXMatches    = True # Cross-match J-PLUS with SDSS, eBOSS targets, etc.
UseSDSSBB       = True  # Use SDSS Broad band filters instead of J-PLUS
PlotColCol      = True # Plot the colour-colour selection of ELG candidates
PlotColMags     = False # Plot color-magnitude diagrams
MakeELGsel      = True  # Create an ELG selection
GetPhotoz       = True # Get photometric redshifts with LePhare
ComputeTwoP     = True  # Compute the angular correlation function of the catalogue

BrowseObjImages = False  # Opens a browser with the object image of each candidate


if PlotComposite:
  z = 0.76
  spec_elg['w'] *= (1 + z)

  for fname in jplus.datasets.jplus_filter_names():
      j_filter = jplus.datasets.fetch_jplus_filter(fname)
      plt.fill(j_filter.wave, j_filter.throughput * 4, ec=jplus.datasets.filter_color(fname), fc=jplus.datasets.filter_color(fname), alpha=0.3)

  plt.ylim([0,3])
  plt.xlim([3000,11000])
  #plt.xlim([3000,9000])
  plt.plot(spec_elg['w'],spec_elg['flux'][:,0],color='blue')

  plt.show()


if LoadCatalogues:
  mag_excess = "AND (m.MAG_APER_3_0[jplus::rSDSS]- m.MAG_APER_3_0[jplus::J0660]) > 0"
  gal_jplus = jplus.datasets.fetch_jplus_objects(mag_type="aperMags", overwrite=False, 
                                                 object_name="allELGs", nchunks=10, mag_limit=[16,24],
                                                extra_where_conds=mag_excess,db='test3')

  gal_sdss  = jplus.datasets.fetch_sdss_qso(mag_type="aperMags", overwrite=False,mag_limit=[16,24])
  gal_phsdss  = jplus.datasets.fetch_sdss_photgalaxies(mag_type="aperMags", overwrite=False,mag_limit=[16,22])
  gal_eboss = jplus.datasets.fetch_eboss_galaxies(mag_type="modelMags", overwrite=False,mag_limit=[16,24])

  for dset in [gal_jplus, gal_eboss]:
      for ifilter in jplus.datasets.jplus_filter_names():
          if ifilter in gal_eboss:
              ind = np.isfinite(dset[ifilter][:,0])
              dset[ifilter][~ind,:] = [-99,-99]
              dset[ifilter][dset[ifilter][:,0]==99,:] = [-99,-99]
              dset[ifilter][dset[ifilter][:,0]==0,:] = [-99,-99]

  if ZpLephare:
    print 're-calibrating zero points using LePhare outputs'
    for ifilter in jplus.datasets.jplus_filter_names(only_nb=True):
        zerop = jplus.tools.compute_zeropoints(mag_type="aperMags", filter_name=ifilter)
        for tile, value in zerop[ifilter]:
            gal_jplus[ifilter][(gal_jplus['tile_id']==tile) & (gal_jplus[ifilter][:,0] != -99),0] -= value



if FindXMatches:
  
  sdssxeboss    = elg.find_xmatches(gal_sdss, gal_eboss, zcond=[.74,.76])
  gal_jplus_xx  = elg.find_xmatches(gal_jplus,sdssxeboss)
  gal_jplus_x   = elg.find_xmatches(gal_jplus, gal_eboss) 
  
  print 'Galaxies from eBOSS in JPLUS: ',len(gal_jplus_x['rJAVA'])
  print 'Galaxies in eBOSS, sdss Spec and J-PLUS:',len(gal_jplus_xx['rJAVA'])
  print 'Galaxies in SDSS Spec and eBOSS Targets', len(sdssxeboss['zspec'])


if UseSDSSBB:
  print 'Replacing J-PUS BBs with SDSS BBs...'
  d2,ind2 = jplus.tools.crossmatch_angular(gal_jplus['coords'],gal_phsdss['coords'],max_distance=3e-4)
  m2 = ((d2 != np.inf))
  
  for f_jplus, f_sdss in zip(jplus.datasets.jplus_filter_names(only_bb = True), jplus.datasets.sdss_filter_names()):
      if AddedPlots:
        print f_jplus, f_sdss
        plt.plot(gal_jplus[f_jplus][m2,0], gal_phsdss[f_sdss][ind2[m2],0]-gal_jplus[f_jplus][m2,0],'.')
        plt.show()
      gal_jplus[f_jplus][m2,:] = gal_phsdss[f_sdss][ind2[m2],:]


if PlotColCol:
  ijlim, rjlim = elg.plot_colcol_sdss_jplus(gal_sdss,gal_jplus)
# ijlim, rjlim = elg.plot_colcol_sdss_jplus(gal_sdss,gal_jplus,xaxis=['gJAVA','rJAVA'],yaxis=['iJAVA','zJAVA'])

  """
  ijlim, rjlim = elg.plot_colcol_sdss_jplus(gal_sdss,gal_jplus,add_muse=False,add_composite=False,add_sdsszline=False)
  ijlim, rjlim = elg.plot_colcol_sdss_jplus(gal_sdss,gal_jplus,xaxis=['gJAVA','rJAVA'],yaxis=['iJAVA','zJAVA'],
  add_muse=False,add_composite=False,add_sdsszline=False)
  
  ijlim, rjlim = elg.plot_colcol_sdss_jplus(gal_sdss,gal_jplus,add_muse=False,add_composite=False)
  ijlim, rjlim = elg.plot_colcol_sdss_jplus(gal_sdss,gal_jplus,xaxis=['gJAVA','rJAVA'],yaxis=['iJAVA','zJAVA'],
  add_muse=False,add_composite=False)

  ijlim, rjlim = elg.plot_colcol_sdss_jplus(gal_sdss,gal_jplus,add_muse=False)
  ijlim, rjlim = elg.plot_colcol_sdss_jplus(gal_sdss,gal_jplus,xaxis=['gJAVA','rJAVA'],yaxis=['iJAVA','zJAVA'],
  add_muse=False)
  """

if MakeELGsel:
  ijlim = 0.5
  rjlim = 0.5
  gal_elgs = elg.make_selection(gal_jplus,ijlim = ijlim, rjlim = rjlim,makeplot=False)  

  nelgs= len(gal_elgs['tile_id'])
  print 'Total number of OII emitter candidates: %d' % nelgs


"""
if PlotColMags:       
  gs = gsc.GridSpec(2,2)
  gs.update(wspace=0.3,hspace=0.3,right=0.85)
  axarr = []
  axarr.append(plt.subplot(gs[0,0]))
  axarr.append(plt.subplot(gs[0,1]))
  axarr.append(plt.subplot(gs[1,0]))

  axarr[0].plot(gal_jplus['J0660'][:,0],gal_jplus['rJAVA'][:,0]-gal_jplus['J0660'][:,0],',',color='gray')
  axarr[1].plot(gal_jplus['rJAVA'][:,0]-gal_jplus['J0660'][:,0],gal_jplus['iJAVA'][:,0] - gal_jplus['J0660'][:,0],',',color='gray')
  axarr[2].plot(gal_jplus['J0660'][:,0],gal_jplus['gJAVA'][:,0]-gal_jplus['rJAVA'][:,0],',',color='gray')

  axarr[0].plot(gal_jplus_x['J0660'][:,0],gal_jplus_x['rJAVA'][:,0]-gal_jplus_x['J0660'][:,0],',',color='blue')
  axarr[1].plot(gal_jplus_x['rJAVA'][:,0]-gal_jplus_x['J0660'][:,0],gal_jplus_x['iJAVA'][:,0] - gal_jplus_x['J0660'][:,0],',',color='blue')
  axarr[2].plot(gal_jplus_x['J0660'][:,0],gal_jplus_x['gJAVA'][:,0]-gal_jplus_x['rJAVA'][:,0],',',color='blue')

  axarr[0].plot(gal_jplus_xx['J0660'][:,0],gal_jplus_xx['rJAVA'][:,0]-gal_jplus_xx['J0660'][:,0],'*',color='pink',markersize=10)
  axarr[1].plot(gal_jplus_xx['rJAVA'][:,0]-gal_jplus_xx['J0660'][:,0],gal_jplus_xx['iJAVA'][:,0] - gal_jplus_xx['J0660'][:,0],'*',
                color='pink',markersize=10)
  axarr[2].plot(gal_jplus_xx['J0660'][:,0],gal_jplus_xx['gJAVA'][:,0]-gal_jplus_xx['rJAVA'][:,0],'*',color='pink',markersize=10)

  axarr[0].plot(jarr,sigma2_col_rj,color='red',linewidth=2)
  axarr[0].plot(gal_jplus['J0660'][icand,0],gal_jplus['rJAVA'][icand,0]-gal_jplus['J0660'][icand,0],',',color='cyan')
  axarr[1].plot(gal_jplus['rJAVA'][icand,0]-gal_jplus['J0660'][icand,0],gal_jplus['iJAVA'][icand,0] - gal_jplus['J0660'][icand,0],
                ',',color='cyan')
  axarr[2].plot(gal_jplus['J0660'][icand,0],gal_jplus['gJAVA'][icand,0]-gal_jplus['rJAVA'][icand,0],',',color='cyan')

  spec_elg = jplus.datasets.fetch_eboss_elg_composite()
  new = 5
  ewarr = np.logspace(0.5,3,num=new)
  amparr = np.logspace(-1,2.,num=10)
  color = plt.cm.coolwarm(np.linspace(0.1,0.9,new))
  z = 0.76
  for iew in range(new):
      i = ewarr[iew]
      print i
      ij = 0
      for j in amparr:
          #print 'EW = ',i
          #print 'amp = ',j
          specline = jplus.datasets.add_line_to_spec(spec_elg,EW=i,obsframe=True, z=z)
          specline['flux'][:,0] *= j
          specline['w'] *= (1+z)
          #print specline['w'].min(),specline['w'].max()
          conv_mags = jplus.datasets.compute_jplus_photometry_singlespec(specline)
          
          axarr[0].plot([conv_mags['J0660'],conv_mags['J0660']],         [(-conv_mags['J0660'] + conv_mags['rJAVA']),(-conv_mags['J0660'] + conv_mags['rJAVA'])],'.',
                       color=color[iew])
          
          r_j = conv_mags['rJAVA'] - conv_mags['J0660']
          i_j = conv_mags['iJAVA'] - conv_mags['J0660']
          g_r = conv_mags['gJAVA'] - conv_mags['rJAVA']
          axarr[1].plot([r_j,r_j],[i_j,i_j],'.',color=color[iew])
          label = "$EW = %.2f \AA{}$" % (i) if ij == 0 else ''
          axarr[2].plot([conv_mags['J0660'],conv_mags['J0660']],[g_r,g_r],'.',color=color[iew],label=label)
          ij += 1
          
          
          

  axarr[0].set_ylim([-0.5,3])
  axarr[0].set_xlim([16,21.5])
  axarr[0].set_xlabel('J0660',fontsize=10)
  axarr[0].set_ylabel('r - J0660',fontsize=10)

  axarr[1].set_ylim([-1,3])
  axarr[1].set_xlim([-1,3])
  axarr[1].set_xlabel('r - J0660',fontsize=10)
  axarr[1].set_ylabel('i - J0660',fontsize=10)

  axarr[2].set_ylim([-1,4])
  axarr[2].set_xlim([16,24])
  axarr[2].set_xlabel('J0660',fontsize=10)
  axarr[2].set_ylabel('g - r',fontsize=10)
  axarr[2].legend(loc='center left', bbox_to_anchor=(1.25, 0.5),fontsize=15)

  plt.show()
"""

if GetPhotoz:
  zphot = get_elg_photoz(gal_elgs)
  import pdb ; pdb.set_trace()

if BrowseObjImages:
  import webbrowser as wb
  fracbrowse = .01  # fraction of objects to browse. Set to 1 to look at them all
  ncand = len(gal_elgs['tile_id'])
  nfrac = int(fracbrowse * ncand)
  flag = np.zeros(nfrac)
  ids = np.random.permutation(np.arange(ncand))[0:nfrac]
  
  print '0: bad; 1: ok; 2: excellent!'
  for iw in range(nfrac):
    objid = gal_elgs['object_id']
    tileid = gal_elgs['tile_id']
    url = 'http://upad.cefca.es/catalogues/jplus-test-2/object_query.html?image=%d&number=%d' % (tileid[ids[iw]],objid[ids[iw]])
    wb.open(url,new=0)
    flag[iw] = input('Flag object?: ')


if ComputeTwoP:

  import CosmoBolognaLib as cbl
  import pymangle
  #%matplotlib inline
  #%matplotlib nbagg
  #

  tiles = jplus.datasets.fetch_jplus_tile_list(overwrite=False)
  stars = jplus.datasets.fetch_jplus_objects_to_mask(overwrite=False)
  #tiles = jplus.select_object(tiles, tiles['depth']<20.5)
  mask = jplus.mangle.Mask(tiles, stars,overwrite=False)
  ran = mask.create_random(10000)

  #hp = jplus.healpix.HealpixMap(ran['coords'], nside=256)
  #hp.plot(title="Random Distribution")

  good = mask.contains(gal_elgs['coords'])## & (gal_elgs['rJAVA'][:,0] < 20.5)
  #gal_elgs= jplus.tools.select_object(gal_elgs, good)

  #good = (mask.contains(gal_sdss['coords']) & (gal_sdss['rSDSS'][:,0] < 21.5) 
  #        & (gal_sdss['zspec'] > 0.7) & (gal_sdss['zspec'] < 0.8 ) ) 
  #gal_elgs= jplus.tools.select_object(gal_sdss, good)

  plt.figure(3)

  plt.plot(ran['coords'][:,0],ran['coords'][:,1],',',color='gray')
  plt.plot(gal_elgs['coords'][:,0],gal_elgs['coords'][:,1],'o',markersize=3,color='royalblue',label='ELG candidates')
  plt.legend()
  plt.xlabel('RA')
  plt.ylabel('DEC')
  plt.savefig('footprint.pdf',bbox_inches='tight')


  # In[27]:


  # In[1]:

  cosmology = cbl.Cosmology()
  ran['redshift'] = np.ones(len(ran['coords'][:,0]))
  gal_elgs['redshift'] = np.ones(len(gal_elgs['coords'][:,0]))

  cat_objs = cbl.Catalogue(cbl.EnumTypes._Galaxy_, cbl.EnumTypes._observedCoordinates_, 
                            gal_elgs['coords'][:,0],gal_elgs['coords'][:,1],gal_elgs['redshift'], cosmology)
  ran_objs = cbl.Catalogue(cbl.EnumTypes._RandomObject_, cbl.EnumTypes._observedCoordinates_, 
                            ran['coords'][:,0],ran['coords'][:,1],ran['redshift'], cosmology)

  angMin = 0.01                #// minimum angular separation 
  angMax = 3.                  #// maximum angular separation 
  nbins = 10                      #// number of bins
  shift = 0.5                  #// shift used to set the bin centre 
  angularUnits = cbl.EnumTypes._degrees_
  twopt = cbl.TwoPointCorrelation1D_angular(cat_objs, ran_objs,cbl.EnumTypes._linear_, angMin, angMax, nbins, shift, 
                                            angularUnits)
  cbl.set_ObjectRegion_SubBoxes(cat_objs,ran_objs,3,3,3)
  twopt.measure(cbl.EnumTypes._Jackknife_,'./')
  twopt.write('./', 'test');

  plt.figure(4)

  plt.errorbar(twopt.xx(), twopt.xi1D(), twopt.error1D(), fmt='o',color='royalblue', label="2pt monopole")
  plt.legend()
  plt.xlabel(r'$\theta [deg]$',fontsize=20)
  plt.ylabel(r'$w(\theta)$',fontsize=20)
  plt.savefig('w_elgcand.pdf',bbox_inches='tight')


