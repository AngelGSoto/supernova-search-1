import numpy as np
from scipy.ndimage import gaussian_filter
from scipy import misc
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS

a = np.arange(50, step=2).reshape((5,5))
print(a)

gaussian_filter(a, sigma=1)
print(a)

#fig = plt.figure()
fig = fits.open('frame-r-006372-5-0073.fits')
plt.gray()  # show the filtered result in grayscale
ax1 = fig.add_subplot(121)  # left side
ax2 = fig.add_subplot(122)  # right side
#ax2 = plt.subplot(1,2,2, projection=WCS(hdu2.header))
#ax2.imshow(hdu2.data, origin='lower', vmin=0, vmax=3)
ascent = fig
result = gaussian_filter(fig, sigma=5)
ax1.imshow(ascent)
ax2.imshow(result)
plt.show()