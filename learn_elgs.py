# Routines to find ELGs using Machine learning algorithms 

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


def LoadSample(tfout, overwrite=False, filtername = 'J0660', linelist = 'x', linename = 'x',sdssxjplus = False,zrange=False):
  filt = jplus.datasets.fetch_jplus_filter(filtername)
  
  if linelist == 'x': 
    linelist = np.array([6563.0, 5007., 4861.0,3727.0])
    linename = ['Halpha', 'OIII', 'Hbeta', 'OII']

  nline = len(linelist)
  meanzline = np.zeros(nline)
  
  allspec = []

  for il in range(nline):
    
    if zrange:
      zr = zrange
    else:
      zr = elg.zline(linelist[il],filt.wave,filt.throughput)
    
    print zr
    
    meanzline[il] = np.mean(zr)
 

  if overwrite == False:
    alls       = pickle.load(open(tfout))
    allspec    = alls['allspec']
    mag_spec = alls['photo_spec']
    nall = len(allspec)

  else:
    # Construct a training set of J-PLUS photometry for ELGs at different redshifts
    # using SDSS, eBOSS, MUSE and VVDS.

    #             Ha    OII     OIII   OII    Hb
    
    for il in range(nline):
   
      if zrange:
        zr = zrange
      else:
        zr = elg.zline(linelist[il],filt.wave,filt.throughput)

      muse_spec   = elg.get_musewide_spec(zr,name=linename[il])
      eboss_spec  = elg.get_eboss_spec(zr,name=linename[il])
      vvds_spec   = elg.get_vvds_spec(zr,name=linename[il])

      nmuse   = len(muse_spec)
      neboss  = len(eboss_spec)
      nvvds   = len(vvds_spec)

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

    if sdssxjplus:
    # Adding jplus magnitudes of xmatches with SDSS spectra
      print 'Loading SDSS-SpecObj and JPLUS catalogues ...'
      gal_sdss  = jplus.datasets.fetch_sdss_galaxies(mag_type="aperMags", 
      overwrite=False,mag_limit=[16,24])
      
      mag_excess = "AND (m.MAG_APER_3_0[jplus::rSDSS]- m.MAG_APER_3_0[jplus::J0660]) > 0"
      gal_jplus = jplus.datasets.fetch_jplus_objects(mag_type="aperMags", overwrite=False, 
                                                 object_name="allELGs", nchunks=10, mag_limit=[16,24],
                                                extra_where_conds=mag_excess,db='test3')

      print 'done with that.'
      allxmatch = elg.find_xmatches(gal_jplus, gal_sdss,zcoord='redshift')
      
      for il in range(nline):
        zr = elg.zline(linelist[il],filt.wave,filt.throughput)
        zsel = np.where((allxmatch['redshift'] > zr[0] ) & (allxmatch['redshift'] < zr[1]))[0]
        nzsel = len(zsel)
        print '%d galaxies from SDSS-Spec at %f<z<%f' % (nzsel, zr[0],zr[1])
        if nzsel > 0:
          for jj in range(nzsel):
            idd = zsel[jj]
            allspec.append({'z': allxmatch['redshift'][idd], 'name':linename[il]})
            
            magjj = {}
            for kk in allxmatch.keys():
              nll = len(allxmatch[kk])
              if nll < idd:
                magjj[kk] = allxmatch[kk] # Stores things like the date
              else:
                magjj[kk] = allxmatch[kk][idd]

            photo_spec.append(magjj)

           


    mag_spec = {}
    
    for i in photo_spec[0].keys():
      mag_ = np.zeros([nall,2])
      for j in range(nall):
        if 'date' in photo_spec[j]: # SDSS-Spec 
          mag_[j,0] = photo_spec[j][i,0]
        else:
          mag_[j,0] = photo_spec[j][i]

      mag_spec[i] = mag_


    print 'writing original spec file'
    dicspec = {'allspec':allspec,'photo_spec':mag_spec}

    with open(tfout,'wb') as outfile:
      pickle.dump(dicspec,outfile,protocol=pickle.HIGHEST_PROTOCOL)

  return allspec, mag_spec


def apply_condition(photo_spec, allspec,rjlim = 0.25, ijlim = 0.25):
  subtrain = {} # This contains the convolved spectra of suitable objects
  zzlist   = []
  namelist = []

  nall = len(photo_spec[photo_spec.keys()[0]])

  for i in photo_spec.keys():
    subtrain[i] = []


  for i in range(nall):

    nozero        = ((photo_spec['J0660'][i,0] != 99) &
                    (photo_spec['gJAVA'][i,0] != 99) &
                    (photo_spec['rJAVA'][i,0] != 99) & 
                    (photo_spec['iJAVA'][i,0] != 99) & 
                    (photo_spec['zJAVA'][i,0] != 99) &
                    (photo_spec['J0861'][i,0] != 99)) 

    maglimits     = ((photo_spec['rJAVA'][i,0] - photo_spec['J0660'][i,0] > rjlim) or
                    (photo_spec['iJAVA'][i,0] - photo_spec['J0660'][i,0] > ijlim))

     
    if nozero and maglimits:
      for ii in photo_spec.keys():
        subtrain[ii].append([photo_spec[ii][i,0],0])

      zzlist.append(allspec[i]['z'])
      namelist.append(allspec[i]['name'])
  
  for i in subtrain.keys():
    subtrain[i] = np.asarray(subtrain[i])

  ntrain = len(subtrain[subtrain.keys()[0]])
  print 'number of objects in training set: %d' % ntrain

  return subtrain, zzlist, namelist


def prepare_sample(samp):

  Colors_samp = []
  nsamp = len(samp['rJAVA'][:,0])
  for i in range(nsamp):
    Colors_samp.append([samp['rJAVA'][i,0] - samp['J0660'][i,0], 
                       samp['iJAVA'][i,0] - samp['J0660'][i,0], 
                       samp['rJAVA'][i,0] - samp['iJAVA'][i,0],
                       samp['zJAVA'][i,0] - samp['iJAVA'][i,0],
                       samp['zJAVA'][i,0] - samp['J0861'][i,0],
                       samp['gJAVA'][i,0] - samp['rJAVA'][i,0]])



  return Colors_samp
  

def learning_elgs(Traindata, Trainfeature, Testdata, EstimatorType = 'Classifier', Plot = True, Scale = True):

  import matplotlib.pyplot as plt
  
  ntrain = len(Trainfeature)
  nelgs  = len(Testdata)
  linename = np.unique(Trainfeature)
  nline = len(linename)

  mlcolors = plt.cm.Set1(np.linspace(0.,1,10))

  from sklearn.neural_network import MLPRegressor
  from sklearn.preprocessing import StandardScaler 
  from sklearn import tree
  from sklearn import svm

  from sklearn.neural_network import MLPClassifier
  from sklearn.ensemble import RandomForestClassifier
  from sklearn.gaussian_process import GaussianProcessClassifier
  from sklearn.gaussian_process.kernels import RBF

  if Scale:
    scaler = StandardScaler()
    scaler.fit(Traindata)
    Traindata = scaler.transform(Traindata)
    Testdata = scaler.transform(Testdata) 

  if EstimatorType == 'Regression':

    clf         = MLPRegressor(solver='adam')
    clf         = clf.fit(Traindata, Trainfeature)
    zelgs_neural= clf.predict(Testdata)

    t_clf       = tree.DecisionTreeRegressor()
    t_clf       = t_clf.fit(Traindata,Trainfeature)
    zelgs_tree  = t_clf.predict(Testdata)
    
    s_clf       = svm.SVR()
    s_clf       = s_clf.fit(Traindata,Trainfeature)
    zelgs_svr   = s_clf.predict(Testdata)

    plt.hist(zelgs_neural,200,range=[0,1],label='Neural z ELGs',normed=True,color='blue',alpha=.35)
    plt.hist(zelgs_tree,200,range=[0,1],label='Decision tree z ELGs',normed=True,color='red',alpha=.35)
    plt.hist(zelgs_svr,200,range=[0,1],label='SVR z ELGs',normed=True,color='magenta',alpha=.35)
    plt.hist(Trainfeature,200,range=[0,1],label='Training set',normed=True,color='green')
    plt.legend()
  
    plt.show()

  if EstimatorType == 'Classifier':

    C = 1.0
    kernel = 1.0 * RBF([1.0,1.0]) # for GPC

    classifiers = { 'MLP'                       : MLPClassifier(solver='lbfgs'), 
                    'Random Forest'             : RandomForestClassifier(n_estimators=10),
                    'SVC'                       : SVC(),
                    'L1 logistic'               : LogisticRegression(C=C, penalty='l1'),
                    'L2 logistic (OvR)'         : LogisticRegression(C=C, penalty='l2'),
                    'L2 logistic (Multinomial)' : LogisticRegression(C=C, solver='lbfgs',multi_class='multinomial'),
                    'GPC'                       : GaussianProcessClassifier(kernel)
                    }

    for index, (name, classifier) in enumerate(classifiers.items()):
    
      classifier.fit(Traindata, Trainfeature)
      y_pred = classifier.predict(Testdata)
      classif_rate = np.mean(y_pred.ravel() == y.ravel()) * 100.0
      print "classif_rate for %s: %f" % (name, classif_rate)
      
      probas = classifier.predict_proba(

      
    totlines_train    = np.zeros(nline) 
    totlines_neural   = np.zeros(nline)
    totlines_rforest  = np.zeros(nline)
    totlines_svm      = np.zeros(nline)

    for il in range(nline):
      totlines_train[il]    = np.float(len(np.where(np.array(Trainfeature) == linename[il])[0]))/ntrain
      totlines_neural[il]   = np.float(len(np.where(lelgs_neural == linename[il])[0]))/nelgs
      totlines_rforest[il]  = np.float(len(np.where(lelgs_rforest == linename[il])[0]))/nelgs
      totlines_svm[il]      = np.float(len(np.where(lelgs_svm == linename[il])[0]))/nelgs

 
    if Plot:
      width=0.35
      fig, ax = plt.subplots()
    
      ax.bar(np.arange(nline),totlines_train,width/3.,label='Training set',color='gray')
      ax.bar(np.arange(nline) + width/3.,totlines_neural,width/3.,label='NN ELGs',
      color=mlcolors[0],alpha=.85)
      
      ax.bar(np.arange(nline) + 2*width/3.,totlines_rforest,width/3.,label='Random forest ELGs',
      color=mlcolors[1],alpha=.85)
      
      ax.bar(np.arange(nline) + width,totlines_svm,width/3.,label='SVM ELGs',
      color=mlcolors[2],alpha=.85)
      
      ax.set_xticks(np.arange(nline) + width)
      ax.set_xticklabels((linename))
      ax.legend(loc='upper left')
      plt.savefig('learn_elgs.pdf',bbox_inches='tight')
  
  
  return [totlines_neural, totlines_rforest, totlines_svm]

   
  
  



