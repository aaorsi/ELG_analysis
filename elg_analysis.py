
import os
jplusdir = '/home/CEFCA/aaorsi/work/j-plus/'
elgdir   = os.getcwd()

import sys
sys.path.append(jplusdir)
os.chdir(jplusdir)

import jplus
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gsc

import elgtools as elg

import pickle

jplus.plotting.setup_text_plots(fontsize=10,usetex=True)
spec_elg = jplus.datasets.fetch_eboss_elg_composite()
matplotlib.rcParams['figure.figsize'] = (12,10)


#### Options #####

LoadELGs        = True   # Overrides ELG generation, and reads it from a python file.
fout = '%s/out/elgs.dat' % elgdir


AddedPlots      = False  # Plot added stuff
LoadCatalogues  = False  # Should be true most of the time
PlotComposite   = False
ZpLephare       = False  # Recalibrate zero points using Lephare
FindXMatches    = False # Cross-match J-PLUS with SDSS, eBOSS targets, etc.
UseSDSSBB       = False  # Use SDSS Broad band filters instead of J-PLUS
PlotColCol      = False # Plot the colour-colour selection of ELG candidates
PlotColMags     = False # Plot color-magnitude diagrams
MakeELGsel      = False  # Create an ELG selection
GetPhotoz       = False # Get photometric redshifts with LePhare
ComputeTwoP     = False  # Compute the angular correlation function of the catalogue

BrowseObjImages = False  # Opens a browser with the object image of each candidate

GetTrainSet     = True
LoadTrainSet    = True
tfout = '%s/out/trainspec.dat' % elgdir


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

  gal_sdss  = jplus.datasets.fetch_sdss_galaxies(mag_type="aperMags", overwrite=False,mag_limit=[16,24])
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
  
  sdssxeboss    = elg.find_xmatches(gal_sdss, gal_eboss, zcond=[.74,.76],zcoord='redshift')
  gal_jplus_xx  = elg.find_xmatches(gal_jplus,sdssxeboss)
  gal_jplus_x   = elg.find_xmatches(gal_jplus, gal_eboss) 
  jplus_sdss    = elg.find_xmatches(gal_jplus,gal_sdss,zcoord='redshift')

  print 'Galaxies from eBOSS in JPLUS: ',len(gal_jplus_x['rJAVA'])
  print 'Galaxies in eBOSS, sdss Spec and J-PLUS:',len(gal_jplus_xx['rJAVA'])
  print 'Galaxies in SDSS Spec and eBOSS Targets', len(sdssxeboss['redshift'])
  #print 'Galaxies in SDSS Spec and JPLUS', len(jplus_sdss['redshift'])


if UseSDSSBB and LoadELGs == False:
  print 'Replacing J-PUS BBs with SDSS BBs...'
  d2,ind2 = jplus.tools.crossmatch_angular(gal_jplus['coords'],gal_phsdss['coords'],max_distance=3e-4)
  m2 = ((d2 != np.inf))
  
  for f_jplus, f_sdss in zip(jplus.datasets.jplus_filter_names(only_bb = True), jplus.datasets.sdss_filter_names()):
      if AddedPlots:
        print f_jplus, f_sdss
        plt.plot(gal_jplus[f_jplus][m2,0], gal_phsdss[f_sdss][ind2[m2],0]-gal_jplus[f_jplus][m2,0],'.')
        plt.show()
      gal_jplus[f_jplus][m2,:] = gal_phsdss[f_sdss][ind2[m2],:]


if PlotColCol and LoadELGs == False:
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

ijlim = 0.25
rjlim = 0.25
if MakeELGsel and LoadELGs is False:
  gal_elgs = elg.make_selection(gal_jplus,ijlim = ijlim, rjlim = rjlim,makeplot=False)  

  nelgs= len(gal_elgs['tile_id'])
  print 'Total number of OII emitter candidates: %d' % nelgs
  
  with open(fout,'wb') as outfile:
    pickle.dump(gal_elgs,outfile,protocol=pickle.HIGHEST_PROTOCOL)
  
  print 'file %s created' % fout


if LoadELGs:
  gal_elgs = pickle.load(open(fout))

  

nelgs = len(gal_elgs['rJAVA'][:,0])


if GetTrainSet:
  
  if LoadTrainSet:
    alls       = pickle.load(open(tfout))
    allspec    = alls['allspec']
    photo_spec = alls['photo_spec']
    nall = len(allspec)

  else:
    # Construct a training set of J-PLUS photometry for ELGs at different redshifts
    # using SDSS, eBOSS, MUSE and VVDS.

    #             Ha    OII     OIII   OII    Hb
    
    j0660 = jplus.datasets.fetch_jplus_filter('J0660')
    
    linelist = np.array([6563.0, 3727.0,5007.,4861.0])
    nline = len(linelist)
    
    allspec = []

    sdss_jplusphot = []

    for il in range(nline):
    
      zr = elg.zline(linelist[il],j0660.wave,j0660.throughput)
      
      muse_spec   = elg.get_musewide_spec(zr)
      eboss_spec  = elg.get_eboss_spec(zr)
      vvds_spec   = elg.get_vvds_spec(zr)

      nmuse = len(muse_spec)
      neboss = len(eboss_spec)
      nvvds = len(vvds_spec)

      for im in range(nmuse):
        allspec.append(muse_spec[im])

      for im in range(neboss):
        allspec.append(eboss_spec[im])
        
      for im in range(nvvds):
        allspec.append(vvds_spec[im])

    
    nall = len(allspec)
    print '%d spectra from all surveys' % nall
    print 'Convolving spectra with J-PLUS filters...'
    photo_spec = []
  #  plt.figure(6)

    for i in range(nall):
      conv = jplus.datasets.compute_jplus_photometry_singlespec(allspec[i])
      photo_spec.append(conv)

    print 'writing original spec file'
    dicspec = {'allspec':allspec,'photo_spec':photo_spec}
    with open(tfout,'wb') as outfile:
      pickle.dump(dicspec,outfile,protocol=pickle.HIGHEST_PROTOCOL)

#    plt.plot(allspec[i]['w'],allspec[i]['flux'])

#    for band in photo_spec[i].iterkeys():
#      filt = jplus.datasets.fetch_jplus_filter(band)
#      w = filt.avgwave()
#      val = photo_spec[i][band]
#      plt.plot([w,w],[val,val],'o')

#    plt.ylim([30,15])
#    plt.title(r'%s' % allspec[i]['survey'])
#    plt.text(6000,20,'z= %f' % allspec[i]['z'])
#    plt.show()


  subtrain = [] # This contains the convolved spectra of suitable objects
  zzlist   = []
  for i in range(nall):

    nozero        = ((photo_spec[i]['J0660'] != 99) &
#                    (photo_spec[i]['gJAVA'] != 99) &
                    (photo_spec[i]['rJAVA'] != 99) & 
                    (photo_spec[i]['iJAVA'] != 99)) 
#                    (photo_spec[i]['zJAVA'] != 99) &
#                    (photo_spec[i]['J0861'] != 99)) 

    maglimits     = ((photo_spec[i]['rJAVA'] - photo_spec[i]['J0660'] > rjlim) or
                    (photo_spec[i]['iJAVA'] - photo_spec[i]['J0660'] > ijlim))

   
    if nozero and maglimits:
      subtrain.append(photo_spec[i])
      zzlist.append(allspec[i]['z'])

  ntrain = len(subtrain)
  print 'number of objects in training set: %d' % ntrain
 

  from sklearn.neural_network import MLPRegressor
  from sklearn.preprocessing import StandardScaler 
  from sklearn import tree
  from sklearn import svm

  Colors_arr = []
  z_arr      = []

  # preparing training data
  for i in range(ntrain):
    Colors_arr.append([subtrain[i]['rJAVA'] - subtrain[i]['J0660'], 
                       subtrain[i]['iJAVA'] - subtrain[i]['J0660'], 
                       subtrain[i]['rJAVA'] - subtrain[i]['iJAVA']])

    z_arr.append(zzlist[i])


  scaler = StandardScaler()
  scaler.fit(Colors_arr)
  Colors_arr = scaler.transform(Colors_arr)

  
  # preparing jplus elg data
  Colors_elgs = []
  for i in range(nelgs):
    Colors_elgs.append([gal_elgs['rJAVA'][i,0] - gal_elgs['J0660'][i,0], 
                       gal_elgs['iJAVA'][i,0] - gal_elgs['J0660'][i,0], 
                       gal_elgs['rJAVA'][i,0] - gal_elgs['iJAVA'][i,0]])




  Colors_elgs = scaler.transform(Colors_elgs) 

  clf         = MLPRegressor(solver='adam')
  clf         = clf.fit(Colors_arr,z_arr)
  zelgs_neural= clf.predict(Colors_elgs)

  t_clf       = tree.DecisionTreeRegressor()
  t_clf       = t_clf.fit(Colors_arr,z_arr)
  zelgs_tree  = t_clf.predict(Colors_elgs)
  
  s_clf       = svm.SVR()
  s_clf       = s_clf.fit(Colors_arr,z_arr)
  zelgs_svr   = s_clf.predict(Colors_elgs)


  
  import matplotlib.pyplot as plt


  plt.hist(zelgs_neural,200,range=[0,1],label='Neural z ELGs',normed=True,color='blue')
  plt.hist(zelgs_tree,200,range=[0,1],label='Decision tree z ELGs',normed=True,color='red')
  plt.hist(zelgs_svr,200,range=[0,1],label='SVR z ELGs',normed=True,color='magenta')
  plt.hist(z_arr,200,range=[0,1],label='Training set',normed=True,color='green')
  plt.legend()
  
  plt.show()
  
  
  
  import ipdb ; ipdb.set_trace()


if GetPhotoz:
  zphot = elg.get_elg_photoz(gal_elgs)
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


