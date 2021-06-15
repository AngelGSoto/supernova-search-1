import numpy as np
from scipy.ndimage import gaussian_filter
from scipy import misc
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.colors import LogNorm

fig = plt.figure()
hdu = fits.open('frame-r-006372-5-0073.fits')
data = hdu[0].data
#plt.gray()  # show the filtered result in grayscale
ax1 = fig.add_subplot(121)  # left side
ax2 = fig.add_subplot(122)  # right side
ascent = data
result = gaussian_filter(data, sigma=5)
print('antes: ', ascent)
print( )
print('depois: ', result)
ax1.imshow(ascent, origin='lower', norm=LogNorm())
ax2.imshow(result, vmin=2.e3, vmax=3.e3, interpolation='nearest', origin='lower', norm=LogNorm())
plt.show()
