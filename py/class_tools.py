# Several routines to deal with the ELG classification

import itertools
from sklearn.ensemble import ExtraTreesClassifier
import matplotlib.pyplot as plt
import numpy as np

def class_to_int(istr):
    if istr == 'Halpha':
        return 0
    elif istr == 'OIII+Hbeta':
        return 1
    if istr == 'OII':
        return 2
    elif istr == 'contaminant':
        return 3
    else:
        print '%s not recognised'%istr
        return -99


def class_to_int_binary(istr, bclass):
    if istr == bclass:
        return 1
    else:
        return 0


def CV_split(data, kfold, ith, binary_class = None):
# data: the dataset
# kfold: k-fold CV. (e.g. kfold = 5 splits dataset into 5 folds)
# ith  : return ith k-fold as test, the rest as training.


    new_train = {}
    test_set = {}
    for k in data.keys():
        test_set[k] = [] # initialize every key as an empty list
        new_train[k] = []

    if binary_class == 'OII':
        classes = ['OII', 'All']
    else:
        classes = np.unique(data['class'])

    print classes
    for ic in classes:

        if binary_class == 'OII' and ic != 'OII':
            sel = np.where(data['class'] != 'OII')[0]
        else:
            sel = np.where(data['class'] == 'OII')[0]

        nsel = len(sel)
        ntest = int(nsel*test_frac)
        ntrain = nsel - ntest
        fold_size = int(sel/kfold)
        id0 = ith*fold_size
        id1 = fold_size*(ith+1) if ith < kfold-1 else nsel-1 #last fold reaches to the end of the array to prevent numerical inaccuracies 
        id_test = sel[id0:id1]
        id_train = np.delete(sel, np.arange(id0,id1))

        for k in test_set.keys():
            for jd in id_test:
                test_set[k].append(data[k][jd])
            for jd in id_train:
                new_train[k].append(data[k][jd])

    for k in test_set.keys():
        test_set[k] = np.array(test_set[k]) # Turn lists into arrays
        new_train[k] = np.array(new_train[k])
    return test_set, new_train



# keeping a fraction away for later validation

def split_trainset(data, test_frac=0.2, binary_class = None):

    new_train = {}
    test_set = {}
    for k in data.keys():
        test_set[k] = [] # initialize every key as an empty list
        new_train[k] = []

    if binary_class == 'OII':
        classes = ['OII', 'All']
    else:
        classes = np.unique(data['class'])

    print classes
    for ic in classes:

        if binary_class == 'OII' and ic != 'OII':
            sel = np.where(data['class'] != 'OII')[0]
        else:
            sel = np.where(data['class'] == 'OII')[0]

        nsel = len(sel)
        ntest = int(nsel*test_frac)
        ntrain = nsel - ntest
        perm = np.random.permutation(sel)
        id_test = perm[0:ntest]
        id_train = perm[ntest:]
        for k in test_set.keys():
            for jd in id_test:
                test_set[k].append(data[k][jd])
            for jd in id_train:
                new_train[k].append(data[k][jd])

    for k in test_set.keys():
        test_set[k] = np.array(test_set[k]) # Turn lists into arrays
        new_train[k] = np.array(new_train[k])
    return test_set, new_train



def get_features(data):
# Function to generate features from galaxy properties

    featnames = ['dm_j0660','J0378','J0395','J0410','J0430','J0515',
             'J0660','J0861','uSDSS','gSDSS','rSDSS','iSDSS','zSDSS']

    feat_arr = []
    err_arr = []

    terms = featnames[1:] # create all colour combinations
    nterms = len(terms)
    ncomb = int(nterms*(nterms-1)/2.)
    print 'All colours:', ncomb
    comb = list(itertools.combinations(terms,2))
    lcomb = list(comb)
    colournames = ['%s - %s'%(x[0], x[1]) for x in list(comb)]

    ndata = len(data['obj'])
    for x in range(ndata):
        fx = []
        ex = []
        for y in featnames:
            fx.append(data[y][x,0])
            ex.append(data[y][x,1])
        for z in range(ncomb):
            fx.append(data[lcomb[z][0]][x,0] - data[lcomb[z][1]][x,0])
            ex.append(np.sqrt((data[lcomb[z][0]][x,1])**2 + (data[lcomb[z][1]][x,1])**2))

        feat_arr.append(fx)
        err_arr.append(ex)
    class_arr = data['class']
    featnames += colournames
#    print featnames, len(featnames)
    nfeat = len(featnames)

    return feat_arr, err_arr, featnames


def feature_importance(x_train, y_train, featnames, figsize=(20,10), n_estimators= 2000):

    forest = ExtraTreesClassifier(n_estimators=n_estimators,
                                  random_state=1323,n_jobs=4)

    x_train = np.array(x_train)
    forest.fit(x_train, y_train)
    importances = forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_],
                 axis=0)
    indices = np.argsort(importances)[::-1]
    inames = [r'$%s$'%featnames[x] for x in indices]

    # Print the feature ranking
    print("Feature ranking:")

    for f in range(x_train.shape[1]):
        print("%d. feature %d (%f): %s" % (f + 1, indices[f], importances[indices[f]], inames[f]))


    # Plot the feature importances of the forest
    plt.figure(1, figsize=figsize)
    plt.title("Feature importances")
    plt.bar(range(x_train.shape[1]), importances[indices],
           color="r", yerr=std[indices], align="center")
    plt.xticks(range(x_train.shape[1]),inames,rotation='vertical')
    plt.xlim([-1, x_train.shape[1]])
    plt.show()

    return indices, importances[indices], inames



def resample_errors(dset, nr=10, scale_up_down = True, exclude_class = None, ref_mag = 'rSDSS', verbose=False,
                          balance_set = False):
# nr: Number of times each object is re-sampled
# scale_up_down: Each new sample will have magnitudes scaled to fall within the range occupied by its class.
# exclude_class: To skip a particular class. 
# ref_mag: The scaling will be done within the range of this magnitude


  ngal = len(dset['obj'])

  if scale_up_down:
  # Before the resampling, find the range covered by ref_mag in each class
    mag_range = {}
    classes = np.unique(dset['class'])
    
    for ic in classes:
      sel = np.where(dset['class'] == ic)[0]
      min_val = dset[ref_mag][sel,0].min()
      max_val = dset[ref_mag][sel,0].max()
      if verbose:
        print 'Class %s range of %s'% (ic, ref_mag)
        print '%.3f - %.3f'% (min_val, max_val)

      mag_range[ic] = [min_val, max_val, len(sel)]

 
  ids_shuffled = np.random.permutation(np.arange(ngal))
  
  rs_dset = {} #dictionary storing new dataset
  for key in dset.keys():
      rs_dset[key] = [] # all keys copied


  for i in range(ngal):
    imag = dset[ref_mag][i,0]
    iclass = dset['class'][i]
    if balance_set:
#      cnr = int(nr*mn/mag_range[iclass][2])
      cnr = int(np.ceil((nr - mag_range[iclass][2])/mag_range[iclass][2]))
    else:
      cnr = nr

    for key in dset.keys():
        rs_dset[key].append(dset[key][i]) #First add the original object
          
    if exclude_class == iclass:
        continue

    for j in range(cnr):
        
        if scale_up_down:   
            dmag = np.random.uniform(mag_range[iclass][0], mag_range[iclass][1]) - imag
        else:
            dmag = 0
        
        for key in dset.keys():
            if len(dset[key][i]) ==2: # if it's a magnitude (this probably also affects the coordinates, but who cares)
                rs_dset[key].append([dmag + np.random.normal(dset[key][i,0],
                                      np.max([dset[key][i,1],0])),dset[key][i,1]])

            else:
                rs_dset[key].append(dset[key][i])

  for key in dset.keys():
      rs_dset[key] = np.array(rs_dset[key])


  
  return rs_dset
  



