'''
Based on the progam of Gabriel.
Original version: covert.py
'''
from astropy.io import fits, ascii
import os
import argparse

def fz2fits(image):
    """
    It converts SPLUS images
    from .fz to .fits
    """
    datos = fits.open(image)[1].data
    heada = fits.open(image)[1].header
    imageout = image[:-2] + 'fits'
    print ('Creating file: ')
    print (imageout)
    fits.writeto(imageout, datos, heada, overwrite=True)

parser = argparse.ArgumentParser(
    description="""Convert file.fz to file.fits""")

parser.add_argument("fzfile", type=str,
                    default="MC0095_F378_swp",
                    help="Name of file, taken the prefix ")

cmd_args = parser.parse_args()
fzfile_ = cmd_args.fzfile + ".fz"

# Using the definition
fz2fits(fzfile_)

