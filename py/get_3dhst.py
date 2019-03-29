# A set of routines to analyse and cross-match 3D-HST data with J-PLUS

from astropy.io import fits

def load_3DHST_cats(RootDir='/home/CEFCA/aaorsi/work/3D-HST/'):
# Load 3D-HST catalogues containing positions, redshifts and line fluxes
  
  GlobalCatDir = '%s/aegis_3dhst.v4.1.cats/Catalog/' % RootDir
  LineCatDir   = '%s/aegis_3dhst_v4.1.5_catalogs/' % RootDir

  Gfile = '%s/aegis_3dhst.v4.1.cat.FITS' % GlobalCatDir
  Lfile = '%s/aegis_3dhst.v4.1.5.linefit.concat.fits' % LineCatDir


  GlobalCat = fits.open(Gfile)
  dcat = GlobalCat[1].data

  LineCat = fits.open(Lfile)
  lcat = LineCat[1].data
  
  return dcat, lcat




def select_3DHST_z(zmin, zmax, LineName=None, LineMin=-1, RootDir=False):
# Routine to select objects from 3D-HST based on redshift and, alternatively, a particular line above a line flux

    
  dcat, lcat = load_3DHST_cats(RootDir) if RootDir else load_3DHST_cats() 
  
# dcat IDs start at 1, dcat[id][0] = 1

  if LineName:
    llim = lcat[LineName] > LineMin
    obj = (lcat['z_max_grism'] >= zmin) & (lcat['z_max_grism'] <= zmax) & llim
  else:
    obj = (lcat['z_max_grism'] >= zmin) & (lcat['z_max_grism'] <= zmax)

  
  objindex = lcat['number'][obj] - 1

  outcat = {'ra': dcat['ra'][objindex], 
           'dec':dcat['dec'][objindex],
           'z': lcat['z_max_grism'][obj]
           }

  if LineName:
    outcat[LineName] = lcat[LineName][obj]


  return outcat











