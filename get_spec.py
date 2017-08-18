# A routine to produce a python file with spectra and convolved magnitudes

from learn_elgs import *
from elgtools import *
import pickle

linelist = [1216.0]
linename = ['Lyalpha']

zrange = [2.0, 2.5]
zrange = False

spec, mags = LoadSample('lyalpha_spec',overwrite= True, filtername = 'J0395', linelist = linelist, linename = linename,
zrange = zrange)




