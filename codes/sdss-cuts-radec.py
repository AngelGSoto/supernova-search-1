import sys
import os
import numpy as np
import pyfits
from optparse import OptionParser
import sqlcl    # tools to querry SDSS DR12

#################################################

def arg_parser():
    parser = OptionParser()
    
    parser.add_option('-r', '--ra',
                      type='float', action='store',
                      help="""Right Ascension of the object center [degree]""")
    parser.add_option('-d', '--dec',
                      type='float', action='store',
                      help="""Declination of the object center [degree]""")    
    parser.add_option("-w", "--width",
                      type='float', action='store',
                      help="""The width of fied of view [degree]""") 
    parser.add_option("-n", "--name",
                      type='string', action='store',
                      help="""object name (optional)""")   
    parser.add_option("-o", "--outfolder",
                      type='string', action='store',
                      help="""output folder (optional)""")    
    parser.add_option("-b", "--bands",
                      type='string', action='store',
                      help="""SDSS bands (optional)""")   
   
    (opts, args) = parser.parse_args()
    return opts, args

#################################################

def xcmd(cmd,verbose):

  if verbose: print('\n'+cmd)

  tmp=os.popen(cmd)
  output=''
  for x in tmp: output+=x
  if 'abort' in output:
    failure=True
  else:
    failure=tmp.close()
  if False:
    print('execution of %s failed' % cmd)
    print('error is as follows',output)
    sys.exit()
  else:
    return output

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
        
#################################################

def getSDSSfields(ra, dec, size):  # all in degree
  
  
  delta = 0.5*size+0.13
  ra_max  =  ra+1.5*delta
  ra_min  =  ra-1.5*delta
  dec_max =  dec+delta
  dec_min =  dec-delta
  
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
  
  querry += "WHERE ra BETWEEN "+str(ra_min)+" and "+str(ra_max)+" and dec BETWEEN "+str(dec_min)+" and "+str(dec_max)
  
  lines = sqlcl.query(querry).readlines()
  N = len(lines)

  field_lst = []
  for i in np.arange(2,N):
      line = str(lines[i])
      line = line.split(',')
      run    = line[1]
      camcol = line[2]
      field  = line[3]
      ra_    = line[4]
      dec_   = line[5]
      field_lst.append([run, camcol, field])
  
  return field_lst

#################################################

def getSDSSfiles(fieldInfo, band, folder):   
  
  run    = fieldInfo[0]
  camcol = fieldInfo[1]
  field  = fieldInfo[2]
  
  fileName = 'frame-'+band + '-''{0:06d}'.format(int(run))+'-'+camcol+'-'+'{0:04d}'.format(int(field))+'.fits.bz2'
  filename = 'frame-'+band + '-''{0:06d}'.format(int(run))+'-'+camcol+'-'+'{0:04d}'.format(int(field))+'.fits'
  http = 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/'
  
  http += run + '/'
  http += camcol + '/'
  http += fileName
  
  with cd(folder):
    xcmd("wget " + http, True)
    xcmd("bzip2 -d " + fileName, True)
    fits_prep(filename, edit_zp=True)

#################################################

def removefix(filename):
  
  
  name_list = filename.split('.')
  N = len(name_list)
  if name_list[N-1] == 'fits':
    
    name =''
    for i in range(N-2):
      name = name + name_list[i] + '.'
    name += name_list[N-2]

  return name

#################################################
# preparing fits files before stacking
# Editing the zero points 
def fits_prep(filename, edit_zp=True):
  
  name = removefix(filename)
  xcmd("cp " + filename + ' ' + name+'.tmp.fits', True) 
  
  hdulist = pyfits.open(filename)
  prihdr = hdulist[1].header    
      
  if edit_zp:
     
     hdulist = pyfits.open(name+'.tmp.fits', mode='update')
     prihdr = hdulist[0].header
     prihdr['BSCALE'] = 1.
     prihdr['BZERO'] = 0.

     expTime = float(prihdr['EXPTIME'])
     nMgy = float(prihdr['NMGY'])

     new_zp = 16.40006562  # mJy / pix
     
     zeroPoint = 22.5
     alpha = 10**(-0.4*(zeroPoint-new_zp))
     
     #hdulist[0].data = data = hdulist[0].data * (3631E-6 )
     hdulist[0].data = data = hdulist[0].data * (alpha )

     prihdr['EXPTIME'] = 1.0
     prihdr['ZP'] = new_zp

     hdulist.flush() 

     cmd = 'mv '+name+'.tmp.fits ' + filename
     xcmd(cmd, True)    
     
#################################################

# Downloading SDSS images and then making mosaics
def sdssget(objName, ra, dec, size, repository): 
        
  fields = getSDSSfields(ra, dec, size)
  
  folder = repository+'tmp'
  if os.path.isdir(folder):
      xcmd('rm -rf '+folder, True)
  xcmd('mkdir '+folder, True)
  
  if not os.path.isdir(repository+objName):
      xcmd('mkdir '+repository+objName, True)
  
  for band in ['u','g','r','i','z']:
      
      for i in range(len(fields)):
        print(' ')  
        print(' # '+str(i+1)+'/'+str(len(fields))+'  ... obj: '+objName+'  band: '+band)
        
        getSDSSfiles(fields[i], band, folder)
      
      xcmd('mv '+folder+'/*fits '+repository+objName, True)
      xcmd('gzip '+repository+objName+'/*fits', True)

  xcmd('rm -rf '+folder, True)

if __name__ == '__main__':
  warning = False
  
  if (len(sys.argv) < 2): 
    print("\nNot enough input arguments ...")
    print("Use -h option for help ... \n", file=sys.stderr)
    exit(1)
  
  
  opts, args =  arg_parser()
  print("\n------------------------------------")
  print("Input Arguments (provided by User)")
  print("------------------------------------")
  print("R.A.       [deg]:", opts.ra)
  print("Dec.       [deg]:", opts.dec)
  print("FOV-width  [deg]:", opts.width)
  print("Object Name  (optional):", opts.name)
  print("Output Foldr (optional):", opts.outfolder)
  print("SDSS Bands   (optional):", opts.bands)
  print("------------------------------------")
  
  ra  = opts.ra 
  dec = opts.dec
  width = opts.width
  
  if ra == None or dec == None or width == None:
    print("\nWrong input parameters ...")
    print("Use -h option for more information ... \n", file=sys.stderr)
    exit(1)      


  bands = []

  if opts.bands != None:
      for b in opts.bands:
        if b in ['u', 'g', 'r', 'i', 'z'] and not b in bands:
            bands.append(b)
  if len(bands) == 0:
     print("\n[Warning] No output SDSS band is given ...")
     print("[Warning] Downloading all bands (ugriz) ...")
     bands = ['u', 'g', 'r', 'i', 'z']
     warning = True
  
  
  if  opts.outfolder == None:
      print("\n[Warning] No output folder were provided ...")
      print("[Warning] Using the current directory to store results ...")
      warning = True
      repository = '.'
  else:
      repository = opts.outfolder
      if not os.path.isdir(repository):
        print("\n[Error] no such directory: " + repository)
        print("Use -h option for help ... \n", file=sys.stderr)
        exit(1)
         
      
      
  if  opts.name == None:
      print("\n[Warning] No object name were given ...")
      warning = True
      ra_st = '{:.4f}'.format(ra)
      dec_st = '{:.4f}'.format(dec)
      if dec < 0:
         objName = 'Obj'+ra_st+'-'+dec_st
      else: 
         objName = 'Obj'+ra_st+'+'+dec_st
      print("Using this name for the object: "+objName)
  else:
      objName = opts.name  
   
  
  if warning:
      
    print("\n------------------------------------")
    print("[Warning] Use the following Parameters ?")
    print("------------------------------------")
    print("R.A.       [deg]:", ra)
    print("Dec.       [deg]:", dec)
    print("FOV-width  [deg]:", width)
    print("Object Name     :", objName)
    if repository == '.':
       print("Output Foldr    : <current folder>")
    else: print("Output Foldr    :", repository)
    print("SDSS Bands    :", bands)
    print("------------------------------------")  
      
    input_var = input("Continue (Y, n)?: ")
    if input_var.upper() == 'Y' or input_var.upper() == '': 
       warning = False
  
  if not warning:  
      sdssget(objName, ra, dec, width, repository+'/')
  else: 
      print("\nTry again ...")
      print("Use the '-h' option for more information ... \n", file=sys.stderr)