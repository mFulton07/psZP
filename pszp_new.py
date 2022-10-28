# Minimal demonstration of using AstroPy specutils to read a plot reduced
# spectra from Liverpool Telescope SPRAT level 2 reduced FITS files
Version='1.0.0'

import sys
import os
import argparse
import statistics
import numpy as np
import pandas as pd
import requests
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord

def PS1catalog(ra,dec,magmin,magmax):

    queryurl = 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/stack?'
    queryurl += 'ra='+str(ra)
    queryurl += '&dec='+str(dec)
    queryurl += '&radius=0.208'
    queryurl += '&columns=[raStack,decStack,gPSFMag,rPSFMag,iPSFMag,zPSFMag,yPSFMag,rKronMag]'
    queryurl += '&nDetections.gte=20&pagesize=25000'

    print('\nQuerying PS1 DR2 Stack for reference stars...\n')

    query = requests.get(queryurl)
    results = query.json()

    if len(results['data']) > 1:
    
        data = np.array(results['data'])

        # Star-galaxy separation: star if rPSFmag - rKronMag < 0.1
        star_data = np.array([x[:-1] for x in data if (x[3] - x[7] < 0.1) and x[7] != -999.0])


        # Remove stars with magnitudes outside the mag range
        for star in star_data:
          for idx, mag in enumerate(star[2:6]):
            if magmax <= mag <= magmin:
              continue
            else:
              star[idx+2] = None
        
        
        # Below is a bit of a hack to remove duplicates
        catalog = SkyCoord(ra=star_data[:,0]*u.degree, dec=star_data[:,1]*u.degree)
        unique_star_data = []
        indices = np.arange(len(star_data))
        used = []
        for star in star_data:
            source = SkyCoord(ra=star[0]*u.degree, dec=star[1]*u.deg)
            d2d = source.separation(catalog)
            catalogmsk = d2d < 2.5*u.arcsec
            indexmatch = indices[catalogmsk]
            for j in indexmatch:
                if j not in used:
                    unique_star_data.append(star_data[j])
                    for k in indexmatch:
                        used.append(k)

        np.savetxt('ref_stars.dat',unique_star_data,fmt='%.8f\t%.8f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f', header='ra\tdec\tg\tr\ti\tz\ty', comments='')
        print('Success! Reference star file created: ref_stars.dat\n')

    else:
        sys.exit('Field not in PS1 DR2! Exiting.\n')



def binMjds(epochs):

    # Create a list of the dection mjds
    binning_mjds = []
    for epoch in epochs:
        binning_mjds.append(epoch[0])
    sorted(binning_mjds)

    # Reduce binning mjds list to unique first occuring mjds, with 0.125 day bins
    binsize = 0.125
    unique_binning_mjds = sorted(list(dict.fromkeys(binning_mjds)))
    for mjd_1 in unique_binning_mjds:
        for mjd_2 in reversed(unique_binning_mjds):
            if 0.0 < abs(mjd_1 - mjd_2) < binsize:
                unique_binning_mjds.remove(mjd_2)
 
    # Bin the epochs keyed on the binning mjds
    binned_epochs = {}
    for mjd in unique_binning_mjds:
        binned_mjd_dets = []
        for epoch in epochs:
            if abs(epoch[0] - mjd) < binsize:
                binned_mjd_dets.append(epoch)
        avr_bin_mjd = statistics.mean([x[0] for x in binned_mjd_dets])
        binned_epochs[avr_bin_mjd] = binned_mjd_dets
            
    return binned_epochs



def binFilters(epochs_binned):

  binned_filtered_epochs = {}
  filtered_epochs = {}
  binning_filters = ['g','r','i','z','y','w']
  for mjdbin, epochs in epochs_binned.items():
    for filterbin in binning_filters:
      skycell_list = []
      for epoch in epochs:
        if epoch[1] == filterbin:
          skycell_list.append(epoch[2])
      if len(skycell_list) > 0:
        filtered_epochs[filterbin] = skycell_list
    binned_filtered_epochs[mjdbin] = filtered_epochs

  return binned_filtered_epochs



def loadSkycellObjects(epochs_skycell_objects):

  for mjdbin, filterlist in epochs_skycell_objects.items():
    for filterbin, skycelllist in filterlist.items():
      skycells_objects_list = []
      for skycell in skycelllist:
        
        # Open up the skycell file and low in the fits data
        hdul = fits.open(skycell)
        data = hdul[1].data
        for row in range(len(data)):
            object = {
                'PSF.RA':round(float(data['RA_PSF'][row]),8),
                'PSF.DEC':round(float(data['DEC_PSF'][row]),8),
                'PSF.CAL_PSF_MAG':round(float(data['CAL_PSF_MAG'][row]),4),
                'PSF.INST_PSF_MAG_SIG':round(float(data['PSF_INST_MAG_SIG'][row]),4),
                'FPA.ZP':hdr['HIERARCH FPA.ZP']
            }
            # Rejecting object PSFs...
            # Where the calibration failed
            # Where the calibration mag error is less than 3-sigma (mag_err <= 0.3)
            if not np.isnan(object['PSF.CAL_PSF_MAG']):
              if object['PSF.INST_PSF_MAG_SIG'] <= 0.3:
                skycells_objects_list.append(object)
      epochs_skycell_objects[mjdbin][filterbin] = skycells_objects_list

  return epochs_skycell_objects



def calSkycellOffsets(epochs_skycell_objects, refstars):

  for mjdbin, filterlist in epochs_skycell_objects.items():
    for filterbin, skycellobjects in filterlist.items():
      for object in skycellobjects:
        c1 = SkyCoord(ra=object['PSF.RA']*u.degree, dec=object['PSF.DEC']*u.degree)
        if filterbin == 'g':
          f=0
        elif filterbin == 'r':
          f=1
        elif filterbin == 'i':
          f=2
        elif filterbin == 'z':
          f=3
        elif filterbin == 'y':
          f=4
        elif filterbin == 'w':
          f=5
        for star in refstars[f]:
          c2 = SkyCoord(ra=star[0]*u.degree, dec=star[1]*u.degree)
          c1c2_sep = c2.separation(c1) * 3600
          a = c1c2_sep.value()
          if c1c2_sep.value() <= 0.1:
            object['CAL.MAG_OFFSET'] = object['PSF.CAL_PSF_MAG'] - star[2]
            object['CAL.MAG_TRUE'] = star[2]
            break
  return



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Pan-STARRS Zeropoint Correction Tool')

    parser.add_argument('--directory', '-d', dest='directory', help='Object directory inside Bundles containing the .cmf files. Mandatory.', default='2022cmc_z_035_bundle', nargs=1, type=str)
    parser.add_argument('--coords', '-c', dest='coords', help='RA and DEC coordinates of object in degrees. Mandatory.', default=[288.26459,19.77341], nargs=2, type=float)
    parser.add_argument('--verbose', '-v', dest='verbose', action='store_true', help='Turn on verbose mode. (default: Off).')

    args = parser.parse_args()

    print(' ')
    print('=======================================================')
    print(' Welcome to PSZP: Pan-STARRS Zeropoint Correction Tool ')  
    print(f'                        V{Version}                        ')
    print('         Written by Michael Fulton (2021-2022)         ')        
    print('=======================================================')
    print(' ')

    # For Debugging
    args.verbose = True


    # Set up thw working directory
    try:
      os.chdir(os.getcwd() + '/Bundles/' + args.directory)
    except:
      sys.exit('\nObject directory invalid! Please make sure your directory exists inside the Bundles directory.\n')


    # Set up the coordinates
    try:
        coords = np.array([float(i) for i in args.coords])
    except:
        sys.exit('\nObject coordinates invalid! Please supply via -c RA DEC in degrees.\n')
    if args.verbose and coords[1] > -50:
        print(f'\nCoordinates found! RA={coords[0]}, DEC={coords[1]}.\n')
    else:
        sys.exit('\nObject coordinates invalid! Please supply a DEC greater than -50 degrees.\n')
      

    # Create catalog of Pan-STARRS reference stars using object coordinates
    if not os.path.exists(os.getcwd() + '/ref_stars.dat'):
      PS1catalog(coords[0],coords[1],21,16)
      refcat = pd.read_csv("ref_stars.dat", sep="\t") 
    else:
      print('\nReference star catalog already exists. If you wish to make a new one, delete the current "ref_stars.dat" file and rerun the command.\n')
      refcat = pd.read_csv("ref_stars.dat", sep="\t")


    # Split the catalog into the separate filters for ease of use later
    refcat_g = []
    refcat_r = []
    refcat_i = []
    refcat_z = []
    refcat_y = []
    refcat_w = []
    refcat = refcat.reset_index()  # make sure indexes pair with number of rows
    for index, row in refcat.iterrows():
      if not np.isnan(row['g']):
        refcat_g.append([row['ra'], row['dec'], row['g']])
      if not np.isnan(row['r']):
        refcat_r.append([row['ra'], row['dec'], row['r']])
      if not np.isnan(row['i']):
        refcat_i.append([row['ra'], row['dec'], row['i']])
      if not np.isnan(row['z']):
        refcat_z.append([row['ra'], row['dec'], row['z']])
      if not np.isnan(row['y']):
        refcat_y.append([row['ra'], row['dec'], row['y']])
      if not np.isnan(row['g']) and not np.isnan(row['r']) and not np.isnan(row['i']):
        psw = round((row['g'] + row['r'] + row['i']) / 3.0, 3)
        refcat_w.append([row['ra'], row['dec'], psw])

    
    # Create a list of all the epochs from the .cmf skycell files in the working directory
    curdir = os.getcwd()
    all_files = [files for files in os.listdir(f'{os.getcwd()}') if os.path.isfile(os.path.join(f'{os.getcwd()}', files) ) ]
    skycells_list = [name for name in all_files if '.cmf' in name]
    print(f'\n{len(skycells_list)} skycell files found in the {args.directory} directory.\n')
    epochs = []
    for skycell in skycells_list:
      # Extract the header from the skycell file
      hdul = fits.open(skycell)
      hdr = hdul[0].header
      epoch_mjd = round(hdr['MJD-OBS'],5)
      epoch_filter = hdr['HIERARCH FPA.FILTER'].strip('.012')
      epochs.append([epoch_mjd, epoch_filter, skycell])

    
    # separate the epochs into filter and MJD bins
    epochs_binned = binMjds(epochs)
    epochs_filtered_binned = binFilters(epochs_binned)


    # Load in the skycell objects for each of the epochs
    epochs_skycell_objects = loadSkycellObjects(epochs_filtered_binned)

    offsets = calSkycellOffsets(epochs_skycell_objects, [refcat_g, refcat_r, refcat_i, refcat_z, refcat_y, refcat_w])



    # Parse the string names of the extensions into extension numbers
    if args.extn_str in ['2','SPEC_NONSS']: 
      extension = 2
      unitName = "adu"
    elif args.extn_str in ['3','SPEC_SS']: 
      extension = 3
      unitName = "adu"
    elif args.extn_str in ['4','NORMFLUX']: 
      # Relative flux normalized to 5500A is dimensionless
      extension = 4
      unitName = "Normalised"
    elif args.extn_str in ['5','FLUX']: 
      extension = 5
      unitName = "erg/s/cm2/A"
    
if args.verbose:
  print (args)

# Convert the host redshift into a number if one is provided
try:
  host_z = float(args.redshift_str)
except:
  host_z = 'none'

# Convert the bin factor into a number if one is provided
try:
  bf = float(args.binfactor_str)
except:
  bf = 'none'

# Convert the bin factor into a number if one is provided
try:
  hpx = float(args.hotpixel_str)
except:
  hpx = 'none'