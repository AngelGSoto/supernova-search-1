from astropy.io import fits
import numpy as np
import pandas as pd
from astropy.convolution import convolve, Gaussian2DKernel
import splusdata
import getpass

username = str(input("Login: "))
#username = 'juliamoliveira'
password = getpass.getpass("Password: ")
conn = splusdata.connect(username, password)

newtable = []

tablefile = './data/selected-gals-vac.csv'
table = pd.read_csv(f'{tablefile}', nrows=10)

for i in range(len(table)):
    ra = table['RA'][i]
    dec = table['DEC'][i]
    id_1 = table['ID'][i]
    id_ = id_1[2:len(id_1) - 1]
    fwhm = table['FWHM_R'][i]
    filename = '%s_%.6f_%.6f' % (id_, ra, dec)

    margin = int((6/0.55) * fwhm)

    hdu_splus = fits.open('./data/splus/' + filename + '-crop.fits')
    hdu_sdss = fits.open('./data/sdss/' + filename + '/' + filename + '-rep.fits')

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

        img = conn.twelve_band_img(ra, dec, margin, noise=0.15, saturation=0.15)
        img.save('./results/colored-stamp/' + filename + '.png')

        data = [id_, ra, dec, res_gauss.max(), res_gauss.min(), np.abs(res_gauss.max() + res_gauss.min())]
        newtable.append(data)

# CANDIDATE TABLE #

cols = ['ID', 'RA', 'DEC', 'MAX', 'MIN', 'RES']
df = pd.DataFrame(newtable, columns=cols)
df.to_csv('./results/candidates.csv', index=False)
