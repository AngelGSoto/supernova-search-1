from astropy.io import fits
from astropy.table import Table

f = fits.open('./data/STRIPE82-0001.fits')[1].data

mask = f['CLASS_STAR_R'] < 0.5
mask &= f['R_auto'] > 13.
mask &= f['e_R_auto'] < 0.15
mask &= (f['PhotoFlag_R'] == 0.) | (f['PhotoFlag_R'] == 2.)

cols = [f['ID'][mask], f['RA'][mask], f['DEC'][mask], f['FWHM_R'][mask]]
names = ['ID', 'RA', 'DEC', 'FWHM_R']
tab = Table(cols, names=names)

tab.write('./data/selected-gals.csv', overwrite=True)
