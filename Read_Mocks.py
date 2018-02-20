import numpy as np
import pylab as pl
from astropy.cosmology import *
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0 = 73.0, Om0 = 0.25)
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os.path
from astropy.coordinates import *
from struct import *

def readmock_chunk_SAM(filename, zspace = None):

        if(os.path.isfile(filename)==True):
                fd = open(filename,'rb')
                ngals = np.fromfile(fd,dtype='i4',count=1)
                gal = np.dtype([('Type',np.int32),('Mvir',np.float32),('pos',np.float32,3),('vel',np.float32,3),\
                ('sfr',np.float32),('sfr_inst',np.float32),('BulgeMass',np.float32),\
                ('DiskMass',np.float32),('Time',np.float32),('redshift',np.float32),\
                ('BlackholeMass',np.float32),('MetalColdGas',np.float32),('ColdGas',np.float32),('MassWeightAge',np.float32),\
                ('ObsMagDust',np.float32,12)])
                distance = np.zeros(ngals)
                GalArr = np.fromfile(fd,dtype=gal,count=ngals)
                distance = np.sqrt(GalArr['pos'][:,0]**2 + GalArr['pos'][:,1]**2+GalArr['pos'][:,2]**2)
                Nlyc_ = np.log10(1.35) + np.log10(GalArr['sfr_inst'][:]) + 53.0
                zzarr = np.linspace(0,20,num=1000)
                darr  = cosmo.comoving_distance(zzarr).value * cosmo.h
                dz = interp1d(darr,zzarr,kind='linear')
                if zspace:
                        vr = ((GalArr['pos'][:,0]*GalArr['vel'][:,0])/distance + (GalArr['pos'][:,1]*GalArr['vel'][:,1])/distance \
                                + (GalArr['pos'][:,2]*GalArr['vel'][:,2])/distance)
                        a_exp = 1 / (1 + GalArr['redshift'][:])
                        H_z = cosmo.H(GalArr['redshift'][:]).value/cosmo.h
                        s = distance + vr / (a_exp*H_z)
                        zsel = np.where(s >= 0)
                        GalArr = GalArr[zsel]
                        ngg = len(GalArr)
                        distance = distance[zsel]
                        vr = vr[zsel]
                        a_exp = a_exp[zsel]
                        H_z = H_z[zsel]
                        s = distance + vr / (a_exp*H_z)
                        redshiftPluspeculiar = dz(s)
                        Nlyc_ = Nlyc_[zsel]

                if zspace:
                        return GalArr, Nlyc_, ngg, s, redshiftPluspeculiar
                else:
                        return GalArr, Nlyc_, ngals
        else:
                print 'Error', filename


def readmock_chunk_PythonCut(dirname, zspace = None):

        if(os.path.isfile(dirname)==True):
                n =  np.dtype([('Type',np.int32)])
                ngals = np.zeros(1,dtype=n)
                fd = open(dirname,'rb')
                ngals = np.load(fd)
                gal = np.dtype([('Type',np.int32),('Mvir',np.float32),('pos',np.float32,3),('vel',np.float32,3),\
                ('sfr',np.float32),('sfr_inst',np.float32),('BulgeMass',np.float32),\
                ('DiskMass',np.float32),('Time',np.float32),('redshift',np.float32),\
                ('BlackholeMass',np.float32),('MetalColdGas',np.float32),('ColdGas',np.float32),('MassWeightAge',np.float32),\
                ('ObsMagDust',np.float32,12)])
                #gal = np.dtype([('Type',np.int32),('Mvir',np.float32),('pos',np.float32,3),('vel',np.float32,3),\
                #('sfr',np.float32),('BulgeMass',np.float32),('DiskMass',np.float32),('Time',np.float32),('redshift',np.float32),\
                #('mag',np.float32),('MetalColdGas',np.float32),('ColdGas',np.float32),('MassWeightAge',np.float32),\
                #('ObsMagDust',np.float32,12)])
                distance = np.zeros(ngals)
                GalArr = np.zeros(ngals,dtype=gal)
                GalArr = np.load(fd)
                Nlyc_ = np.log10(1.35) + np.log10(GalArr['sfr_inst'][:]) + 53.0
                zzarr = np.linspace(0,20,num=1000)
                darr  = cosmo.comoving_distance(zzarr).value * cosmo.h
                dz = interp1d(darr,zzarr,kind='linear')
                if zspace:
                        distance = np.sqrt(GalArr['pos'][:,0]**2 + GalArr['pos'][:,1]**2+GalArr['pos'][:,2]**2)
                        vr = ((GalArr['pos'][:,0]*GalArr['vel'][:,0])/distance + \
                                        (GalArr['pos'][:,1]*GalArr['vel'][:,1])/distance + \
                                        (GalArr['pos'][:,2]*GalArr['vel'][:,2])/distance)
                        a_exp = 1 / (1 + GalArr['redshift'][:])
                        H_z = cosmo.H(GalArr['redshift'][:]).value/cosmo.h
                        s = distance + vr / (a_exp*H_z)
                        zsel = np.where((s >= 0))
                        GalArr = GalArr[zsel]
                        ngg = len(GalArr)
                        distance = distance[zsel]
                        vr = vr[zsel]
                        a_exp = a_exp[zsel]
                        H_z = H_z[zsel]
                        s = distance + vr / (a_exp*H_z)
                        redshiftPluspeculiar = dz(s)
                        Nlyc_ = Nlyc_[zsel]
                if zspace:
                        return GalArr, Nlyc_, ngg, s , redshiftPluspeculiar
                else:
                        return GalArr, Nlyc_, ngals
        else:
                print 'Error with the path file', dirname

def readmock_chunk_SAM_provisional(filename, zspace = None):

	if(os.path.isfile(filename)==True):
                fd = open(filename,'rb')
                ngals = np.fromfile(fd,dtype='i4',count=1)
                gal = np.dtype([('Type',np.int32),('Mvir',np.float32),('pos',np.float32,3),('vel',np.float32,3),\
                ('sfr',np.float32),('BulgeMass',np.float32),\
                ('DiskMass',np.float32),('Time',np.float32),('redshift',np.float32),\
                ('mag',np.float32),('MetalColdGas',np.float32),('ColdGas',np.float32),('MassWeightAge',np.float32),\
                ('ObsMagDust',np.float32,12)])
                distance = np.zeros(ngals)
                GalArr = np.fromfile(fd,dtype=gal,count=ngals)
                distance = np.sqrt(GalArr['pos'][:,0]**2 + GalArr['pos'][:,1]**2+GalArr['pos'][:,2]**2)
                Nlyc_ = np.log10(1.35) + np.log10(GalArr['sfr'][:]) + 53.0
                zzarr = np.linspace(0,20,num=1000)
                darr  = cosmo.comoving_distance(zzarr).value * cosmo.h
                dz = interp1d(darr,zzarr,kind='linear')
                if zspace:
                        vr = ((GalArr['pos'][:,0]*GalArr['vel'][:,0])/distance + (GalArr['pos'][:,1]*GalArr['vel'][:,1])/distance \
                                + (GalArr['pos'][:,2]*GalArr['vel'][:,2])/distance)
                        a_exp = 1 / (1 + GalArr['redshift'][:])
                        H_z = cosmo.H(GalArr['redshift'][:]).value/cosmo.h
                        s = distance + vr / (a_exp*H_z)
                        zsel = np.where(s >= 0)
                        GalArr = GalArr[zsel]
                        ngg = len(GalArr)
                        distance = distance[zsel]
                        vr = vr[zsel]
                        a_exp = a_exp[zsel]
                        H_z = H_z[zsel]
                        s = distance + vr / (a_exp*H_z)
                        redshiftPluspeculiar = dz(s)
                        Nlyc_ = Nlyc_[zsel]

                if zspace:
                        return GalArr, Nlyc_, ngg, s, redshiftPluspeculiar
                else:
                        return GalArr, Nlyc_, ngals
        else:
                print 'Error', filename

