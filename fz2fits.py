# Based on the program of Gabriel.
# Original version: covert.py
from astropy.io import fits

def fz2fits(image):
    data = fits.open(image)[1].data
    header = fits.open(image)[1].header
    imageout = image[:-2] + 'fits'
    fits.writeto(imageout, data, header, overwrite=True)

