import numpy as np
from scipy.ndimage import gaussian_filter
from scipy import misc
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.colors import LogNorm

fig = plt.figure()
#hdu = fits.open('frame-r-006372-5-0073.fits')
hdu = fits.open('iDR3.STRIPE82-0001.000260-crop.fits')
data = hdu[0].data
#plt.gray()  # show the filtered result in grayscale
ax1 = fig.add_subplot(121)  # left side
ax2 = fig.add_subplot(122)  # right side
ascent = data
result = gaussian_filter(data, sigma=1)
#print('antes: ', ascent)
#print( )
#print('depois: ', result)
ax1.imshow(ascent, vmin=0, vmax=3)
ax2.imshow(result, vmin=0, vmax=3)
plt.show()
