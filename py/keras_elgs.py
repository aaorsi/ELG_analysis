# Script based on local_learn_elgs.ipynb

import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation
from keras.optimizers import SGD
import numpy as np
import jplus
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import pickle

#load databases

import pickle
dset = pickle.load(open('dataset.data'))
print dset.keys()
print

training_features = np.asarray(dset['tfeatures'])
training_class    = np.asarray(dset['tclass'])
test_features     = np.asarray(dset['vfeatures'])
test_class        = np.asarray(dset['vclass'])

x_train, x_test, y_train, y_test = train_test_split(
training_features, training_class, test_size=0.2, random_state=None)  # leave out 10% out for test/validation, whatever that's called

print 'Unique classes in training set: ',np.unique(y_train)
print 'Number of objects in training set:', len(y_train)

Scaledata = True
if Scaledata:
    print 'scaling data...',
    scaler = StandardScaler()
    scaler.fit(x_train)
    x_train = scaler.transform(x_train)
    x_test  = scaler.transform(x_test)
    print 'done'

nfeat = len(x_train[0])
print nfeat
#blippi

#Example sets

def class_to_int(istr):
    if istr == 'Halpha':
        return 0
    if istr == 'OII':
        return 1
    elif istr == 'OIII+Hbeta':
        return 2
    elif istr == 'contaminant':
        return 3
    else:
        print '%s not recognised'%istr
        return -99


y_train_int = np.asarray([class_to_int(x) for x in y_train])
y_test_int  = [class_to_int(x) for x in y_test]


y2_train = keras.utils.to_categorical(y_train_int, num_classes=len(np.unique(y_train_int)))
y2_test = keras.utils.to_categorical(y_test_int, num_classes=len(np.unique(y_test_int)))

# Do a grid to search for best parameters

from keras.wrappers.scikit_learn import KerasClassifier
from sklearn.model_selection import GridSearchCV
from keras import optimizers

# All parameter gradients will be clipped to
# a maximum norm of 1.
sgd = optimizers.SGD(lr=0.01, clipnorm=1.)

import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.regularizers import l2

def build_classifier(optim='sgd',drop=0.0,Neurons=5,init='uniform',lambd=0.0):
   #   from keras.optimizers import SGD
    model = Sequential()
    model.add(Dense(Neurons, input_dim=nfeat, kernel_initializer=init, activation='sigmoid', W_regularizer=l2(lambd)))
    model.add(Dropout(drop))
    model.add(Dense(Neurons, kernel_initializer=init, activation='sigmoid', W_regularizer=l2(lambd)))
    model.add(Dropout(drop))
    model.add(Dense(4, kernel_initializer=init, activation='sigmoid',W_regularizer=l2(lambd)))
    # Compile model
    model.compile(loss='sparse_categorical_crossentropy', optimizer=optim, metrics=['accuracy'])
    return model

seed = 7
np.random.seed(seed)

model = KerasClassifier(build_fn=build_classifier, verbose=0)

optim = ['sgd']#,'adam']#,'rmsprop','sgd'] #['rmsprop', 'adam']
init = ['uniform']#['glorot_uniform', 'normal', 'uniform']
epochs = [10,20]
batches = [200]#,200,500]#,500]#, 500]
lambd   = [0]#,1e-1]#,5e-3,1e-3]
dropout = [0,.1]#1e-3,1e-2,1e-1]#,0.5,0.75]
Neurons = [20]
param_grid = dict(optim=optim, epochs=epochs, batch_size=batches, init=init, lambd=lambd,drop=dropout,
                 Neurons=Neurons)
grid = GridSearchCV(estimator=model, param_grid=param_grid,n_jobs=-1,scoring='f1_macro')
grid_result = grid.fit(np.asarray(x_train), y_train_int)

best= grid_result.best_estimator_
print grid_result.best_params_
print best.get_params()
print grid_result.scorer_
print grid_result.best_score_
y_pred = best.predict(x_train, verbose=1)
if len(np.unique(y_pred)) <  4:
  print 'something is wrong with the classification in y_pred, check it out'
  print 'unique labels predicted:',np.unique(y_pred)

from sklearn.metrics import f1_score

print  'F1_Scores for best estimator'
print f1_score(np.asarray(y_train_int),np.asarray(y_pred),average=None)

import pandas as pd
df = grid_result.cv_results_
# summarize results
print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_))
means = grid_result.cv_results_['mean_test_score']
stds = grid_result.cv_results_['std_test_score']
params = grid_result.cv_results_['params']

for mean, stdev, param in zip(means, stds, params):
    print("%f (%f) with: %r" % (mean, stdev, param))



