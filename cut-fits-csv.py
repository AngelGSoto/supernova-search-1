'''
Cuting images FITS
Based on extract-image.py from Henney program and pyFIST.py
'''
from __future__ import print_function
import numpy as np
import json
import os
from astropy.io import fits
from astropy import wcs
from astropy.wcs import WCS
from astropy import coordinates as coord
from astropy import units as u 
import argparse
import sys
import pandas as pd

parser = argparse.ArgumentParser(
    description="""Cut images from fits files""")

parser.add_argument("source", type=str,
                    default="1000001-JPLUS-02363-v2_J0660",
                    help="Name of source (prefix for files) ")

parser.add_argument("-t", "--table", type=str,
                    default="selected-gals.csv",
                    help="Coordinates table [csv]")

parser.add_argument("-m", "--margin", type=float,
                    default=10.0,
                    help="Margin around object (arcseconds)")

parser.add_argument("-d", "--debug", action="store_true",
                    help="Print out verbose debugging info about each line in region file")

args = parser.parse_args()
regionfile = args.source + ".fits"

path1 = "../"
try:
    hdu = fits.open(os.path.join(path1, regionfile))
except FileNotFoundError:
    hdu = fits.open(regionfile)
    
crop_coords_unit=u.degree

# Definition 
def HMS(angle): 
    """
    Convert angle (which has astropy.units units) to an HMS string
    """
    return coord.Angle(angle).to_string(u.hour, sep=":")

def DMS(angle): 
    """
    Convert angle (which has astropy.units units) to a DMS string
    """
    return coord.Angle(angle).to_string(u.degree, sep=":")

table = args.table
csv = pd.read_csv(f'./{table}')
    
#for i in range(len(csv)):
for i in range(1):
    ra = csv['RA'][0]
    dec = csv['DEC'][0]
    print(ra, dec)

    crop_c = coord.SkyCoord(ra, dec, unit="deg")
    print(crop_c)

    w = wcs.WCS(hdu[0].header)
    #print(w)

    ##########################################################
    ## Find minimum and maximum RA, DEC ######################
    ##########################################################
    margin = args.margin * u.arcsecond
    # I had ignore the cos(delta) factor I mean, considering cos(delta)~1 (I should fix that)
    ra1 =  coord.Angle(crop_c.ra.min() - margin ) 
    ra2 =  coord.Angle(crop_c.ra.max() + margin ) 
    dec1 = coord.Angle(crop_c.dec.min() - margin) 
    dec2 = coord.Angle(crop_c.dec.max() + margin) 

    print("RA range:", HMS(ra1), HMS(ra2))
    print("Dec range:", DMS(dec1), DMS(dec2))

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