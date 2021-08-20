from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import convolve, Gaussian2DKernel, Box2DKernel
from astropy.visualization import MinMaxInterval
from scipy.ndimage import gaussian_filter

id98 = 'iDR3.STRIPE82-0001.000098_0.152506_-1.409033'
id112 = 'iDR3.STRIPE82-0001.000112_359.389440_-1.408599'
id128 = 'iDR3.STRIPE82-0001.000128_359.817735_-1.408191'
id148 = 'iDR3.STRIPE82-0001.000148_359.515791_-1.407501'
id165 = 'iDR3.STRIPE82-0001.000165_359.898317_-1.406894'

id = id128

hdu_splus = fits.open(id + '-crop.fits')
hdu_sdss = fits.open(id + '-rep.fits')
#hdu_sdss = fits.open('frame-r-007778-6-0316.fits')


### SPLUS ###

fig = plt.figure()
fig.suptitle('SPLUS')
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

# Resampling
ax1.imshow(hdu_splus[0].data, origin='lower')
ax1.set_title('Resampled', fontsize=10)

# Normalization
max_splus = hdu_splus[0].data.max()
norm_splus = hdu_splus[0].data / max_splus
ax2.imshow(norm_splus, origin='lower')
ax2.set_title('Normalized', fontsize=10)

# Gaussian
gauss_kernel = Gaussian2DKernel(1)
result1 = convolve(hdu_splus[0].data, gauss_kernel)
max_splus = result1.max()
result1 = result1 / max_splus
ax3.imshow(result1, origin='lower')
ax3.set_title('Gaussian', fontsize=10)
plt.show()

'''
# Gaussian
interval1 = MinMaxInterval()
vmin1, vmax1 = interval1.get_limits(hdu_splus[0].data)
result1 = gaussian_filter(hdu_splus[0].data, sigma=1)
max_gauss1 = result1.max()
result1 = result1 / max_gauss1
ax3.imshow(result1, origin='lower')
ax3.set_title('Gaussian', fontsize=10)
plt.show()
'''

print('\nSPLUS')
print('Original (min and max): %s and %s' % (hdu_splus[0].data.min(), hdu_splus[0].data.max()))
print('Final (min and max): %s and %s' % (result1.min(), result1.max()))
print('----------------------------------------------------')


### SDSS ###

fig = plt.figure()
fig.suptitle('SDSS')
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

# Resampling
ax1.imshow(hdu_sdss[0].data, origin='lower')
ax1.set_title('Resampled', fontsize=10)

# Normalization
max_norm = hdu_sdss[0].data.max()
norm_sdss = hdu_sdss[0].data / max_norm
ax2.imshow(norm_sdss, origin='lower')
ax2.set_title('Normalized', fontsize=10)

# Gaussian
gauss_kernel = Gaussian2DKernel(1)
result2 = convolve(hdu_sdss[0].data, gauss_kernel)
max_sdss = result2.max()
result2 = result2 / max_sdss
ax3.imshow(result2, origin='lower')
ax3.set_title('Gaussian', fontsize=10)
plt.show()

'''
# Gaussian
interval2 = MinMaxInterval()
vmin2, vmax2 = interval2.get_limits(hdu_sdss[0].data)
result2 = gaussian_filter(hdu_sdss[0].data, sigma=1)
max_gauss2 = result2.max()
result2 = result2 / max_gauss2
ax3.imshow(result2, origin='lower')
ax3.set_title('Gaussian', fontsize=10)
plt.show()
'''

print('\nSDSS')
print('Original (min and max): %s and %s' % (hdu_sdss[0].data.min(), hdu_sdss[0].data.max()))
print('Final (min and max): %s and %s' % (result2.min(), result2.max()))
print('----------------------------------------------------')


### RESIDUE ###

res = hdu_splus[0].data - hdu_sdss[0].data
res_norm = norm_splus - norm_sdss
res_gauss = result1 - result2

fig = plt.figure()
fig.suptitle('Residue (SPLUS - SDSS)')
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
ax1.imshow(res, origin='lower')
ax1.set_title('Residue', fontsize=8)
ax2.imshow(res_norm, origin='lower')
ax2.set_title('Normalized Residue', fontsize=8)
ax3.imshow(res_gauss, origin='lower')
ax3.set_title('Gaussian Residue', fontsize=8)
#fig.colorbar(im)
plt.show()


print('\nResidue (min and max): %s and %s' % (res.min(), res.max()))
print('Standard Deviation: %s \nMedian: %s' % (np.std(res), np.median(res)))
print('----------------------------------------------------')
print('\nResidue (Normalized) (min and max): %s and %s' % (res_norm.min(), res_norm.max()))
print('Standard Deviation: %s \nMedian: %s' % (np.std(res_norm), np.median(res_norm)))
print('----------------------------------------------------')
print('\nResidue (Gauss) (min and max): %s and %s' % (res_gauss.min(), res_gauss.max()))
print('Standard Deviation: %s \nMedian: %s' % (np.std(res_gauss), np.median(res_gauss)))

fig = plt.figure()
fig.suptitle(id)
ax1 = fig.add_subplot(331)
ax2 = fig.add_subplot(332)
ax3 = fig.add_subplot(333)
ax4 = fig.add_subplot(334)
ax5 = fig.add_subplot(335)
ax6 = fig.add_subplot(336)
ax7 = fig.add_subplot(337)
ax8 = fig.add_subplot(338)
ax9 = fig.add_subplot(339)

ax1.imshow(hdu_splus[0].data, origin='lower')
ax1.set_title('Resampled SPLUS', fontsize=8)
ax2.imshow(norm_splus, origin='lower')
ax2.set_title('Normalized SPLUS', fontsize=8)
ax3.imshow(result1, origin='lower')
ax3.set_title('Gaussian SPLUS', fontsize=8)

ax4.imshow(hdu_sdss[0].data, origin='lower')
ax4.set_title('Resampled SDSS', fontsize=8)
ax5.imshow(norm_sdss, origin='lower')
ax5.set_title('Normalized SDSS', fontsize=8)
ax6.imshow(result2, origin='lower')
ax6.set_title('Gaussian SDSS', fontsize=8)

im7 = ax7.imshow(res, origin='lower')
ax7.set_title('Resampled Residue', fontsize=8)
plt.colorbar(im7, ax=ax7)
im8 = ax8.imshow(res_norm, origin='lower')
ax8.set_title('Normalized Residue', fontsize=8)
plt.setp(ax8.get_yticklabels(), visible=False)
plt.colorbar(im8, ax=ax8)
im9 = ax9.imshow(res_gauss, origin='lower')
ax9.set_title('Gaussian Residue', fontsize=8)
plt.setp(ax9.get_yticklabels(), visible=False)
plt.colorbar(im9, ax=ax9)

plt.tight_layout()
plt.show()