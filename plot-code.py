from astropy.io import fits
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.convolution import convolve, Gaussian2DKernel
import os

'''
From the candidates table (analysis.py code), plots the process of the method and the final results
'''

tablefile = './results/candidates.csv'
table = pd.read_csv(f'{tablefile}')

for i in range(len(table)):
    ra = table['RA'][i]
    dec = table['DEC'][i]
    id_1 = table['ID'][i]
    id_ = id_1[2:len(id_1) - 1]
    filename = '%s_%.6f_%.6f' % (id_, ra, dec)

    hdu_splus = fits.open('./data/splus/' + filename + '-crop.fits')

    pasta = './data/sdss/' + filename
    try:
        caminhos = [os.path.join(pasta, nome) for nome in os.listdir(pasta)]
    except FileNotFoundError:
        continue
    arquivos = [arq for arq in caminhos if os.path.isfile(arq)]
    sdss_files = [arq for arq in arquivos if arq.lower().endswith(".fits")]

    for i in range(len(sdss_files)):
        filename_sdss = sdss_files[i]
        hdu_sdss = fits.open(filename_sdss)

        print('\n*****************************************************')
        print('Objeto: ', filename)
        print('RA: ', ra)
        print('DEC: ', dec)
        print('----------------------------------------------------')

        ### SPLUS ###

        # Normalization
        max_splus = hdu_splus[0].data.max()
        norm_splus = hdu_splus[0].data / max_splus

        # Gaussian
        gauss_kernel = Gaussian2DKernel(1)
        result1 = convolve(hdu_splus[0].data, gauss_kernel)
        max_splus = result1.max()
        result1 = result1 / max_splus

        ### SDSS ###

        # Normalization
        max_norm = hdu_sdss[0].data.max()
        norm_sdss = hdu_sdss[0].data / max_norm

        # Gaussian
        gauss_kernel = Gaussian2DKernel(1)
        result2 = convolve(hdu_sdss[0].data, gauss_kernel)
        max_sdss = result2.max()
        result2 = result2 / max_sdss

        ### RESIDUE ###

        res = hdu_splus[0].data - hdu_sdss[0].data
        res_norm = norm_splus - norm_sdss
        res_gauss = result1 - result2

        ### PLOT 9 IMAGES ###

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
        plt.savefig('./results/images/' + filename + '.png')
        #plt.show()


        ### PLOT 3 IMAGES ###

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
        plt.savefig('./results/images/gaussian/' + filename_sdss.split('/')[-1] + '.png')
        #plt.show()
