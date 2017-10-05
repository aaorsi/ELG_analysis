# A routine to produce a python file with spectra and convolved magnitudes

from learn_elgs import *
from elgtools import *
import pickle

linelist = [3727.0, 4860.0, 5007.0, 4959.,6563.0]
linename = ['OII3727', 'Hb', 'OIII5007', 'OIII4959','Ha']

zrange = [2.0, 2.5]
zrange = False

spec, mags = LoadSample('J0660_spec',overwrite= True, filtername = 'J0660', linelist = linelist, linename = linename,
zrange = zrange)




