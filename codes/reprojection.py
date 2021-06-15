from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from reproject import reproject_interp
#from astropy.utils.data import get_pkg_data_filename

hdu1 = fits.open('frame-r-006372-5-0073.fits')[0]
hdu2 = fits.open('iDR3.STRIPE82-0001.000260-crop.fits')[0]

ax1 = plt.subplot(1,2,1, projection=WCS(hdu2.header))
ax1.imshow(hdu1.data, origin='lower', vmin=-2.e-4, vmax=5.e-4)
ax1.coords['ra'].set_axislabel('RA')
ax1.coords['dec'].set_axislabel('DEC')

ax2 = plt.subplot(1,2,2, projection=WCS(hdu2.header))
ax2.imshow(hdu2.data, origin='lower', vmin=0, vmax=3)
ax2.coords['ra'].set_axislabel('RA')
ax2.coords['dec'].set_axislabel('DEC')
#ax2.coords['dec'].set_axislabel_position('r')
#ax2.coords['dec'].set_ticklabel_position('r')

plt.show()
 
array, footprint = reproject_interp(hdu1, hdu2.header)
#print(array, footprint)

ax1 = plt.subplot(1,2,1, projection=WCS(hdu1.header))
ax1.imshow(array, origin='lower', vmin=-2.e-4, vmax=5.e-4)
ax1.coords['ra'].set_axislabel('Right Ascension')
ax1.coords['dec'].set_axislabel('Declination')

ax2 = plt.subplot(1,2,2, projection=WCS(hdu2.header))
ax2.imshow(hdu2.data, origin='lower', vmin=0, vmax=3)
ax2.coords['ra'].set_axislabel('Right Ascension')
ax2.coords['dec'].set_axislabel('Declination')
# ax2.coords['dec'].set_axislabel_position('r')
# ax2.coords['dec'].set_ticklabel_position('r')

plt.show()

#fits.writeto('teste.fits', array, hdu1.header, overwrite=True)