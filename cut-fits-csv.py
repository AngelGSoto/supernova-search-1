'''
Cuting images FITS
Based on extract-image.py from Henney program and pyFIST.py
'''
from __future__ import print_function
import os
from astropy.io import fits
from astropy import wcs
from astropy.wcs import WCS
from astropy import coordinates as coord
from astropy import units as u
import sys
import pandas as pd

#inputs
table = './data/selected-gals.csv'
outfolder_SPLUS = './data/splus'
margin = 10.

csv = pd.read_csv(f'./{table}')

# for i in range(len(csv)):
for i in range(50):
    ra = csv['RA'][i]
    dec = csv['DEC'][i]
    id_ = csv['ID'][i]

name = '%s_%.6f_%.6f' % (id_, ra, dec)

file = name + ".fz"

path1 = outfolder_SPLUS
try:
    hdu = fits.open(os.path.join(path1, file))
except FileNotFoundError:
    hdu = fits.open(file)

crop_coords_unit = u.degree
crop_c = coord.SkyCoord(ra, dec, unit=(u.deg, u.deg))
w = wcs.WCS(hdu[1].header)

##########################################################
## Find minimum and maximum RA, DEC ######################
##########################################################

margin = margin * u.arcsec

# I had ignore the cos(delta) factor I mean, considering cos(delta)~1 (I should fix that)
ra1 = coord.Angle(crop_c.ra.min() - margin)
ra2 = coord.Angle(crop_c.ra.max() + margin)
dec1 = coord.Angle(crop_c.dec.min() - margin)
dec2 = coord.Angle(crop_c.dec.max() + margin)

###########################################################
## Rectangle in RA, Dec that encloses object with margin ##
###########################################################
coords = [
    [ra1.deg, dec1.deg],
    [ra1.deg, dec2.deg],
    [ra2.deg, dec1.deg],
    [ra2.deg, dec2.deg],
]

##########################################################
## Convert to pixel coords and find enclosing rectangle ##
##########################################################
pix_coords = w.wcs_world2pix(coords, 0)
x = pix_coords[:, 0]
y = pix_coords[:, 1]
i1, i2 = int(x.min()), int(x.max()) + 1
j1, j2 = int(y.min()), int(y.max()) + 1

ny, nx = hdu[1].data.shape
i1 = max(0, i1)
i2 = min(i2, nx - 1)
j1 = max(0, j1)
j2 = min(j2, ny - 1)
print("Extracted image window: [{}:{}, {}:{}]".format(i1, i2, j1, j2))

#########################################################
## Extract window from image and adjust WCS info ########
#########################################################
outhdu = fits.PrimaryHDU(
    data=hdu[1].data[j1:j2, i1:i2],
    header=hdu[1].header
)
outhdu.header["CRPIX1"] -= i1
outhdu.header["CRPIX2"] -= j1

####################
# Save the new file##
####################
outfile = file.replace(".fits", "-crop.fits")
# new_hdu = fits.PrimaryHDU(hdu[0].data, header=hdu[0].header)
outhdu.writeto(outfile, output_verify="fix", overwrite=True)
