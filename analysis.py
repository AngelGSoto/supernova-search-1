from astropy.io import fits
from astropy.table import Table
import numpy as np
import pandas as pd
from astropy.convolution import convolve, Gaussian2DKernel
import splusdata
import getpass
import os
import gzip
import getpass

'''
From the object selection table, analyses the residue and selects candidates
'''

username = str(input("Login: "))
password = getpass.getpass("Password: ")
conn = splusdata.connect(username, password)

newtable = []

tablefile = './data/selected-gals-vac.csv'
table = pd.read_csv(f'{tablefile}', nrows=400)

for i in range(len(table)):
    ra = table['RA'][i]
    dec = table['DEC'][i]
    id_1 = table['ID'][i]
    id_ = id_1[2:len(id_1) - 1]
    fwhm = table['FWHM_R'][i]
    filename = '%s_%.6f_%.6f' % (id_, ra, dec)

    margin = int((6/0.55) * fwhm)

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

        # SPLUS #

        max_splus = hdu_splus[0].data.max()
        norm_splus = hdu_splus[0].data / max_splus

        gauss_kernel = Gaussian2DKernel(1)
        gauss_splus = convolve(hdu_splus[0].data, gauss_kernel)
        max_splus = gauss_splus.max()
        gauss_splus = gauss_splus / max_splus

        # SDSS #

        max_norm = hdu_sdss[0].data.max()
        norm_sdss = hdu_sdss[0].data / max_norm

        gauss_kernel = Gaussian2DKernel(1)
        gauss_sdss = convolve(hdu_sdss[0].data, gauss_kernel)
        max_sdss = gauss_sdss.max()
        gauss_sdss = gauss_sdss / max_sdss

        # RESIDUE #

        res = hdu_splus[0].data - hdu_sdss[0].data
        res_norm = norm_splus - norm_sdss
        res_gauss = gauss_splus - gauss_sdss

        if (np.abs(res_gauss.max() + res_gauss.min()) > 0.1):
            candidate = True

            # Twelve band images
            img = conn.twelve_band_img(ra, dec, margin, noise=0.15, saturation=0.15)
            img.save('./results/colored-stamp/' + filename + '.png')

            data = [id_, ra, dec, res_gauss.max(), res_gauss.min(), np.abs(res_gauss.max() + res_gauss.min())]
            newtable.append(data)

# CANDIDATE TABLE #

cols = ['ID', 'RA', 'DEC', 'MAX', 'MIN', 'RES']
df = pd.DataFrame(newtable, columns=cols)
df.to_csv('./results/candidates.csv', index=False)