from __future__ import print_function
import splusdata
import pandas as pd
import numpy as np
import sys
import os
import sqlcl
from astropy.io import fits
from astropy import wcs, coordinates as coord, units as u
from astropy.wcs import WCS


def main():
    # Change for different SDSS cuts
    width = 0.075
    tablefile = './data/selected-gals.csv'
    table = pd.read_csv(f'{tablefile}', nrows=50)
    bands = ['r']
    outfolder_SDSS = './data/sdss'
    outfolder_SPLUS = './data/splus'

    #splusCuts(table)
    #sdssCuts(width, table, bands, outfolder_SDSS)

    #importar programa de converter pra fits

    cutFits(table, outfolder_SPLUS)


# --------------------------------------------------------------------
# SPLUS CUTS
# Code by Gustavo Schwarz (www.github.com/Schwarzam/splusdata)

def splusCuts(table):
    conn = splusdata.connect('juliamoliveira', '10203040')

    for key, value in table.iterrows():
        hdu = conn.get_cut(value.RA, value.DEC, 128, 'R')
        hdu.writeto('./data/splus/%s_%.6f_%.6f.fz' % (value.ID, value.RA, value.DEC))

    print()
    print('SPLUS stamps have been downloaded.')
    print()


# --------------------------------------------------------------------
# SDSS CUTS
# Adapted from a code by Ehsan Kourkchi (www.github.com/ekourkchi/SDSS_get)

def sdssCuts(width, table, bands, outfolder):
    def xcmd(cmd, verbose):
        if verbose: print('\n' + cmd)

        tmp = os.popen(cmd)
        output = ''
        for x in tmp: output += x
        if 'abort' in output:
            failure = True
        else:
            failure = tmp.close()

    class cd:
        """Context manager for changing the current working directory"""

        def __init__(self, newPath):
            self.newPath = os.path.expanduser(newPath)

        def __enter__(self):
            self.savedPath = os.getcwd()
            os.chdir(self.newPath)

        def __exit__(self, etype, value, traceback):
            os.chdir(self.savedPath)

    # ...................................................................

    def getSDSSfields(ra, dec, size):  # all in degree
        delta = size
        ra_max = ra + (delta / np.cos(abs(np.radians(dec))))
        ra_min = ra - (delta / np.cos(abs(np.radians(dec))))
        dec_max = dec + delta
        dec_min = dec - delta

        querry = """

         SELECT
         fieldID,
         run, 
         camCol, 
         field,
         ra, 
         dec,
         run,
         rerun 
         FROM Field   
         """

        querry += "WHERE ra BETWEEN " + str(ra_min) + " and " + str(ra_max) + " and dec BETWEEN " + str(
            dec_min) + " and " + str(dec_max)

        lines = sqlcl.query(querry).readlines()
        N = len(lines)

        field_lst = []
        for i in np.arange(2, N):
            line = str(lines[i])
            line = line.split(',')
            run = line[1]
            camcol = line[2]
            field = line[3]
            ra_ = line[4]
            dec_ = line[5]
            field_lst.append([run, camcol, field])

        return field_lst

    # ...................................................................

    def getSDSSfiles(fieldInfo, band, folder):
        run = fieldInfo[0]
        camcol = fieldInfo[1]
        field = fieldInfo[2]

        fileName = 'frame-' + band + '-''{0:06d}'.format(int(run)) + '-' + camcol + '-' + '{0:04d}'.format(
            int(field)) + '.fits.bz2'
        filename = 'frame-' + band + '-''{0:06d}'.format(int(run)) + '-' + camcol + '-' + '{0:04d}'.format(
            int(field)) + '.fits'
        http = 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/'

        http += run + '/'
        http += camcol + '/'
        http += fileName

        with cd(folder):
            xcmd("wget " + http, True)
            xcmd("bzip2 -d " + fileName, True)
            fits_prep(filename, edit_zp=True)

    # ...................................................................

    def removefix(filename):
        name_list = filename.split('.')
        N = len(name_list)
        if name_list[N - 1] == 'fits':
            name = ''
            for i in range(N - 2):
                name = name + name_list[i] + '.'
                name += name_list[N - 2]

        return name

    # ...................................................................
    # preparing fits files before stacking
    # Editing the zero points
    def fits_prep(filename, edit_zp=True):
        name = removefix(filename)
        xcmd("cp " + filename + ' ' + name + '.tmp.fits', True)

        hdulist = fits.open(filename)
        prihdr = hdulist[1].header

        if edit_zp:
            hdulist = fits.open(name + '.tmp.fits', mode='update')
            prihdr = hdulist[0].header
            prihdr['BSCALE'] = 1.
            prihdr['BZERO'] = 0.

            expTime = float(prihdr['EXPTIME'])
            nMgy = float(prihdr['NMGY'])

            new_zp = 16.40006562  # mJy / pix

            zeroPoint = 22.5
            alpha = 10 ** (-0.4 * (zeroPoint - new_zp))

            # hdulist[0].data = data = hdulist[0].data * (3631E-6 )
            hdulist[0].data = data = hdulist[0].data * (alpha)

            prihdr['EXPTIME'] = 1.0
            prihdr['ZP'] = new_zp

            hdulist.flush()

            cmd = 'mv ' + name + '.tmp.fits ' + filename
            xcmd(cmd, True)

    # ...................................................................

    # Downloading SDSS images and then making mosaics
    def sdssget(objName, ra, dec, size, bands, repository):
        fields = getSDSSfields(ra, dec, size)
        folder = repository + 'tmp'
        if os.path.isdir(folder):
            xcmd('rm -rf ' + folder, True)
        xcmd('mkdir ' + folder, True)

        if not os.path.isdir(repository + objName):
            xcmd('mkdir ' + repository + objName, True)

        for band in bands:
            for i in range(len(fields)):
                print(' ')
                print(' # ' + str(i + 1) + '/' + str(len(fields)) + '  ... obj: ' + objName + '  band: ' + band)

                getSDSSfiles(fields[i], band, folder)

            xcmd('mv ' + folder + '/*fits ' + repository + objName, True)
            xcmd('gzip ' + repository + objName + '/*fits', True)

        xcmd('rm -rf ' + folder, True)

    if __name__ == '__main__':
        warning = True

        for i in range(len(table)):
            ra = table['RA'][i]
            dec = table['DEC'][i]
            id = table['ID'][i]

            id_st = '{:s}'.format(id)
            ra_st = '{:.6f}'.format(ra)
            dec_st = '{:.6f}'.format(dec)
            if dec < 0:
                objName = id_st + '_' + ra_st + '_' + '-' + dec_st + '_sdss'
            else:
                objName = id_st + '_' + ra_st + '_' + '+' + dec_st + '_sdss'
            print("Using this name for the object: " + objName)

            if outfolder == None:
                print("\n[Warning] No output folder were provided ...")
                print("[Warning] Using the current directory to store results ...")
                warning = True
                repository = '.'
            else:
                repository = outfolder
                if not os.path.isdir(repository):
                    print("\n[Error] no such directory: " + repository)
                    print("Use -h option for help ... \n", file=sys.stderr)
                    exit(1)

            sdssget(objName, ra, dec, width, bands, repository + '/')

        print()
        print('SDSS stamps have been downloaded.')
        print()


# --------------------------------------------------------------------
# CUTTING FITS IMAGES
# Based on extract-image.py from Henney program and pyFIST.py
# Adapted from Luis Angel Soto's code (www.github.com/AngelGSoto)

def cutFits(table, outfolder, margin):

    #margin =

    for i in range(len(table)):
        ra = table['RA'][0]
        dec = table['DEC'][0]
        id_ = table['ID'][0]

        name = '%s_%.6f_%.6f' % (id_, ra, dec)

        filename = name + ".fits"

        path1 = outfolder
        try:
            hdu = fits.open(os.path.join(path1, filename))
        except FileNotFoundError:
            hdu = fits.open(filename)

        crop_coords_unit = u.degree
        crop_c = coord.SkyCoord(ra, dec, unit=(u.deg, u.deg))
        w = wcs.WCS(hdu[0].header)

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

        ny, nx = hdu[0].data.shape
        i1 = max(0, i1)
        i2 = min(i2, nx - 1)
        j1 = max(0, j1)
        j2 = min(j2, ny - 1)
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
        # Save the new file##
        ####################
        #outfile = filename.replace(".fits", "-crop.fits") ---- para não reescrever

        outhdu.writeto(filename, output_verify="fix", overwrite=True)

# --------------------------------------------------------------------

main()