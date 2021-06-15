from astropy.io import fits
from astropy.table import Table

f = fits.open("STRIPE82-0001.fits")[1].data

mask = f['CLASS_STAR_G'] < 0.1
mask &= f['G_auto'] > 16
mask &= f['e_G_auto'] < 0.15

cols = [f['ID'][mask], f['RA'][mask], f['DEC'][mask]]
names = ['ID', 'RA', 'DEC']
tab = Table(cols, names=names)

tab.write('selected-gals.csv')