'''
Cuting images FITS
Based on extract-image.py from Henney program and pyFIST.py
'''
from __future__ import print_function
from astropy.io import fits
from astropy import wcs
from astropy import coordinates as coord
from astropy import units as u 
import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description="""Cut images from fits files""")

parser.add_argument("source", type=str,
                    default="1000001-JPLUS-02363-v2_J0660",
                    help="Name of source (prefix for files) ")

parser.add_argument("-t", "--table", type=str,
                    default="./data/selected-gals.csv",
                    help="Coordinates table [csv]")

parser.add_argument("-m", "--margin", type=float,
                    default=10.0,
                    help="Margin around object (arcseconds)")

parser.add_argument("-d", "--debug", action="store_true",
                    help="Print out verbose debugging info about each line in region file")

args = parser.parse_args()
regionfile = './data/splus/' + args.source + ".fits"

hdu = fits.open(regionfile)

crop_coords_unit=u.degree

table = args.table
csv = pd.read_csv(f'./{table}')

for i in range(1):
    ra = csv['RA'][i]
    dec = csv['DEC'][i]
    id_ = csv['ID'][i]

    name = '%s_%.6f_%.6f' % (id_, ra, dec)
    filename = name + ".fits"

    crop_c = coord.SkyCoord(ra, dec, unit="deg")

    w = wcs.WCS(hdu[0].header)

    ##########################################################
    ## Find minimum and maximum RA, DEC ######################
    ##########################################################
    margin = args.margin * u.arc

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
    x = pix_coords[:,0]
    y = pix_coords[:,1]
    i1, i2 = int(x.min()), int(x.max()) + 1
    j1, j2 = int(y.min()), int(y.max()) + 1

    ny, nx = hdu[0].data.shape
    i1 = max(0, i1)
    i2 = min(i2, nx-1)
    j1 = max(0, j1)
    j2 = min(j2, ny-1)
    print("Extracted image window: [{}:{}, {}:{}]".format(i1, i2, j1, j2))

    #########################################################
    ## Extract window from image and adjust WCS info ########
    #########################################################
    outhdu = fits.PrimaryHDU(
       data=hdu[0].data[j1:j2, i1:i2],
       header=hdu[0].header
      )
    outhdu.header["CRPIX1"] -= i1
    outhdu.header["CRPIX2"] -= j1
    
    #################### 
    #Save the new file##
    ####################
    outfile = regionfile.replace(".fits", "-crop.fits")
    #new_hdu = fits.PrimaryHDU(hdu[0].data, header=hdu[0].header)
    outhdu.writeto(outfile, output_verify="fix", overwrite=True)
