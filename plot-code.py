from astropy.io import fits
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.convolution import convolve, Gaussian2DKernel

tablefile = './data/selected-gals.csv'
table = pd.read_csv(f'{tablefile}', nrows=50)

for i in range(len(table)):
    ra = table['RA'][i]
    dec = table['DEC'][i]
    id_1 = table['ID'][i]
    id_ = id_1[2:len(id_1) - 1]
    filename = '%s_%.6f_%.6f' % (id_, ra, dec)

    hdu_splus = fits.open('./data/splus/' + filename + '-crop.fits')
    hdu_sdss = fits.open('./data/sdss/' + filename + '/' + filename + '-rep.fits')

    print('\n*****************************************************')
    print('Objeto: ', filename)
    print('RA: ', ra)
    print('DEC: ', dec)
    print('----------------------------------------------------')

    ### SPLUS ###
    '''
    fig = plt.figure()
    fig.suptitle('SPLUS')
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    '''

    # Resampling
    '''
    ax1.imshow(hdu_splus[0].data, origin='lower')
    ax1.set_title('Resampled', fontsize=10)
    '''

    # Normalization
    max_splus = hdu_splus[0].data.max()
    norm_splus = hdu_splus[0].data / max_splus
    '''
    ax2.imshow(norm_splus, origin='lower')
    ax2.set_title('Normalized', fontsize=10)
    '''

    # Gaussian
    gauss_kernel = Gaussian2DKernel(1)
    result1 = convolve(hdu_splus[0].data, gauss_kernel)
    max_splus = result1.max()
    result1 = result1 / max_splus
    '''
    ax3.imshow(result1, origin='lower')
    ax3.set_title('Gaussian', fontsize=10)
    plt.show()
    '''

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

    '''
    fig = plt.figure()
    fig.suptitle('SDSS')
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    '''

    # Resampling
    '''
    ax1.imshow(hdu_sdss[0].data, origin='lower')
    ax1.set_title('Resampled', fontsize=10)
    '''

    # Normalization
    max_norm = hdu_sdss[0].data.max()
    norm_sdss = hdu_sdss[0].data / max_norm
    '''
    ax2.imshow(norm_sdss, origin='lower')
    ax2.set_title('Normalized', fontsize=10)
    '''

    # Gaussian
    gauss_kernel = Gaussian2DKernel(1)
    result2 = convolve(hdu_sdss[0].data, gauss_kernel)
    max_sdss = result2.max()
    result2 = result2 / max_sdss
    '''
    ax3.imshow(result2, origin='lower')
    ax3.set_title('Gaussian', fontsize=10)
    plt.show()
    '''

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

    '''
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
    '''


    print('\nResidue (min and max): %s and %s' % (res.min(), res.max()))
    print('Standard Deviation: %s \nMedian: %s' % (np.std(res), np.median(res)))
    print('----------------------------------------------------')
    print('\nResidue (Normalized) (min and max): %s and %s' % (res_norm.min(), res_norm.max()))
    print('Standard Deviation: %s \nMedian: %s' % (np.std(res_norm), np.median(res_norm)))
    print('----------------------------------------------------')
    print('\nResidue (Gauss) (min and max): %s and %s' % (res_gauss.min(), res_gauss.max()))
    print('Standard Deviation: %s \nMedian: %s' % (np.std(res_gauss), np.median(res_gauss)))

    '''
    #PLOT 9 IMAGENS

    fig = plt.figure(figsize=(6.5, 6))
    fig.suptitle(filename, y=0.965)
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
    ax1.set_title('Resampled SPLUS', fontsize=9)
    plt.setp(ax1.get_xticklabels(), visible=False)

    ax2.imshow(hdu_sdss[0].data, origin='lower')
    ax2.set_title('Resampled SDSS', fontsize=9)
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)

    im3 = ax3.imshow(res, origin='lower')
    ax3.set_title('Resampled Residue', fontsize=9)
    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    plt.colorbar(im3, ax=ax3)

    ax4.imshow(norm_splus, origin='lower')
    ax4.set_title('Normalized SPLUS', fontsize=9)
    plt.setp(ax4.get_xticklabels(), visible=False)

    ax5.imshow(norm_sdss, origin='lower')
    ax5.set_title('Normalized SDSS', fontsize=9)
    plt.setp(ax5.get_yticklabels(), visible=False)
    plt.setp(ax5.get_xticklabels(), visible=False)

    im6 = ax6.imshow(res_norm, origin='lower')
    ax6.set_title('Normalized Residue', fontsize=9)
    plt.setp(ax6.get_yticklabels(), visible=False)
    plt.setp(ax6.get_xticklabels(), visible=False)
    plt.colorbar(im6, ax=ax6)

    im7 = ax7.imshow(result1, origin='lower')
    ax7.set_title('Gaussian SPLUS', fontsize=9)

    im8 = ax8.imshow(result2, origin='lower')
    ax8.set_title('Gaussian SDSS', fontsize=9)
    plt.setp(ax8.get_yticklabels(), visible=False)

    im9 = ax9.imshow(res_gauss, origin='lower')
    ax9.set_title('Gaussian Residue', fontsize=9)
    plt.setp(ax9.get_yticklabels(), visible=False)
    plt.colorbar(im9, ax=ax9)

    plt.tight_layout()
    plt.savefig('./results/' + filename + '.png')
    #plt.show()
    '''

    #'''
    #PLOT 3 IMAGENS

    fig = plt.figure(figsize=(6.5, 2.5))
    fig.suptitle(filename, y=0.93)
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    im1 = ax1.imshow(result1, origin='lower')
    ax1.set_title('SPLUS', fontsize=9)

    im2 = ax2.imshow(result2, origin='lower')
    ax2.set_title('SDSS', fontsize=9)
    plt.setp(ax2.get_yticklabels(), visible=False)

    im = ax3.imshow(np.arange(100).reshape((10, 10)))
    im3 = ax3.imshow(res_gauss, origin='lower')
    ax3.set_title('Residue (SPLUS - SDSS)', fontsize=9)
    plt.setp(ax3.get_yticklabels(), visible=False)
    divider = make_axes_locatable(ax3)
    cax = ax3.inset_axes([1.05, 0, 0.05, 1], transform=ax3.transAxes)
    plt.colorbar(im3, cax=cax)

    plt.tight_layout()
    plt.savefig('./results/gaussian/' + filename + '.png')
    #plt.show()
    #'''