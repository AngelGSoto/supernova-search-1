import splusdata
import pandas as pd
import sys
import os
import numpy as np
from astropy.io import fits as pyfits
from optparse import OptionParser
import sqlcl

def main():
    # Change for different SDSS cuts
    width = 0.075
    table = './data/selected-gals.csv'
    bands = ['r']
    outfolder = './data/sdss'

    splusCuts(table)
    sdssCuts(width, table, bands, outfolder)


# --------------------------------------------------------------------
# SPLUS CUTS

def splusCuts(table):
    conn = splusdata.connect('juliamoliveira', '10203040')

    df = pd.read_csv(table, nrows=50)

    for key, value in df.iterrows():
        hdu = conn.get_cut(value.RA, value.DEC, 128, 'R')
        hdu.writeto('./data/splus/%s_%.6f_%.6f.fz' % (value.ID, value.RA, value.DEC))

    print()
    print('SPLUS stamps have been downloaded.')
    print()


# --------------------------------------------------------------------
# SDSS CUTS

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

        hdulist = pyfits.open(filename)
        prihdr = hdulist[1].header

        if edit_zp:
            hdulist = pyfits.open(name + '.tmp.fits', mode='update')
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

        csv = pd.read_csv(f'{table}')

        # for i in range(len(csv)):
        for i in range(50):
            ra = csv['RA'][i]
            dec = csv['DEC'][i]
            id = csv['ID'][i]

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

main()