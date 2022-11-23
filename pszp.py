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
import matplotlib.pyplot as plt
import csv
from astropy.io import fits
from astropy.coordinates import SkyCoord
from alive_progress import alive_bar

def PS1catalog(ra,dec):

    queryurl = 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/stack?'
    queryurl += 'ra='+str(ra)
    queryurl += '&dec='+str(dec)
    queryurl += '&radius=0.164'
    queryurl += '&columns=[raStack,decStack,gPSFMag,rPSFMag,iPSFMag,zPSFMag,yPSFMag,rKronMag]'
    queryurl += '&nDetections.gte=25&pagesize=10000'

    print('\nQuerying PS1 DR2 Stack for reference stars...\n')

    query = requests.get(queryurl)
    results = query.json()

    if len(results['data']) > 1:
    
        data = np.array(results['data'])

        # Star-galaxy separation: star if rPSFmag - rKronMag < 0.1
        star_data = np.array([x[:-1] for x in data if (x[3] - x[7] < 0.1) and x[7] != -999.0])

        # Remove stars with magnitudes outside the desired mag range:
        # Too bright (m>16) and we will see saturation effects in brighter stars
        # Too faint (m>22) and we willl see noise amongst the fainter stars + distant galaxies
        for star in star_data:
          for idx, mag in enumerate(star[2:]):
            if 16 <= mag <= 22:
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
        avr_bin_mjd = round(statistics.mean([x[0] for x in binned_mjd_dets]),3)
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
            # Where the calibration mag value is outside of the 16.5 - 20.5 range (Too bright and we will see saturation effects, too faint and we will see noise)
            if not np.isnan(object['PSF.CAL_PSF_MAG']):
              if object['PSF.INST_PSF_MAG_SIG'] <= 0.3:
                if 15.5 <= object['PSF.CAL_PSF_MAG'] <= 22.5:
                  skycells_objects_list.append(object)

      epochs_skycell_objects[mjdbin][filterbin] = skycells_objects_list

  return epochs_skycell_objects



def calSkycellOffsets(epochs_skycell_objects, refstars):

  for mjdbin, filterlist in epochs_skycell_objects.items():
    for filterbin, skycellobjects in filterlist.items():

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
      
      calibrated_objects_list = []
      extreme_offsets_outliers_list = []
      with alive_bar(len(skycellobjects), title=f'Calibrating {mjdbin}_{filterbin}: ') as bar:   
        for object in skycellobjects:
          for star in refstars[f]:
            match_offset_deg = ((star[0] - object['PSF.RA'])**2 + (star[1] - object['PSF.DEC'])**2)**0.5
            if (match_offset_deg*3600 <= 0.1):
              object['CAL.MAG_TRUE'] = star[2]
              object['CAL.MAG_OFFSET'] = object['PSF.CAL_PSF_MAG'] - star[2]
              calibrated_objects_list.append(object)
              if abs(object['CAL.MAG_OFFSET']) > 0.5:
                extreme_offsets_outliers_list.append(object)
              break
          bar()
      mag_offsets_list = [x['CAL.MAG_OFFSET'] for x in calibrated_objects_list]
      mean_clip = np.nanmean(mag_offsets_list)
      sigma_clip = 4 * np.nanstd(mag_offsets_list)
      calibrated_clipped_objects_list = [object for object in calibrated_objects_list if (mean_clip - sigma_clip <= object['CAL.MAG_OFFSET'] <= mean_clip + sigma_clip)]
      epochs_skycell_objects[mjdbin][filterbin] = calibrated_clipped_objects_list

  return epochs_skycell_objects, extreme_offsets_outliers_list



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Pan-STARRS Zeropoint Correction Tool')

    parser.add_argument('--directory', '-d', dest='directory', help='Object directory inside Bundles containing the .cmf files. Mandatory.', type=str)
    parser.add_argument('--coords', '-c', dest='coords', help='RA and DEC coordinates of object in degrees. Mandatory.', nargs=2, type=float)
    parser.add_argument('--verbose', '-v', dest='verbose', action='store_true', help='Turn on verbose mode. (default: Off).')

    args = parser.parse_args()

    print(' ')
    print('=======================================================')
    print(' Welcome to PSZP: Pan-STARRS Zeropoint Correction Tool ')  
    print(f'                        V{Version}                        ')
    print('         Written by Michael Fulton (2021-2022)         ')        
    print('=======================================================')
    print(' ')


    # Set up thw working directory
    try:
      os.chdir(os.getcwd() + '/Bundles/' + args.directory)
      print(f'\nWorking path: {os.getcwd()}/.\n')
    except FileNotFoundError as not_found:
      sys.exit(f'\nPath "{not_found.filename}" invalid! Please make sure your sub-directory exists inside the Bundles directory.\n')


    # Set up the coordinates
    try:
        coords = np.array([float(i) for i in args.coords])
    except:
        sys.exit('Object coordinates invalid! Please supply via -c RA DEC in degrees.\n')
    if args.verbose and coords[1] > -30:
        print(f'Coordinates parsed: RA={coords[0]}, DEC={coords[1]}.\n')
    else:
        sys.exit('Object coordinates invalid! Please supply a DEC greater than -30 degrees.\n')
      

    # Create catalog of Pan-STARRS reference stars using object coordinates
    if not os.path.exists(os.getcwd() + '/ref_stars.dat'):
      PS1catalog(coords[0],coords[1])
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
    print(f'\n{len(skycells_list)} skycell files found in the {args.directory} sub-directory. Loading in files...\n')
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
    print(f'Skycell objects loaded.\n')


    # Calculate magnitude offsets between skycell objects and reference stars
    offset_objects, outlier_objects = calSkycellOffsets(epochs_skycell_objects, [refcat_g, refcat_r, refcat_i, refcat_z, refcat_y, refcat_w])
    print(f'Skycell offsets calculated.\n')


    # Plot offsets
    print(f'\nPlotting offsets...\n')
    for mjdbin, filterlist in offset_objects.items():
      for filterbin, skycellobjects in filterlist.items():
        
        # Plot lightcurve comparisons
        fig, ax = plt.subplots()
        fig.set_size_inches(11.7, 8.3, forward=True)
        plt.tight_layout()
        
        ax.set_xlim([15.0, 22.0])
        ax.set_ylim([-0.6, +0.6])
        ax.set_xlabel('Reference Star Mag', fontsize=24)
        ax.set_ylabel('Mag Offset [Skycell - Reference]', fontsize=24)
        ax.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=16)
        ax.set_title(f'Difference Magnitude Offsets for: MJD {int(mjdbin)} {filterbin}-band', fontsize=24)

        x = [object['CAL.MAG_TRUE'] for object in skycellobjects]
        y = [object['CAL.MAG_OFFSET'] for object in skycellobjects]
        ax.scatter(x, y, s=100, marker='x', color='black', zorder=1)

        offest_mean = round(np.nanmean(y),4)
        offest_median = round(np.nanmedian(y),4)
        offest_stdev = round(np.nanstd(y),4)
        
        ax.text(15.05, 0.03, 'Observed\nFainter', color='r', ha='left', va='bottom', fontsize=16)
        ax.axhline(linewidth=3, color='r')
        ax.text(15.05, -0.03, 'Observed\nBrighter', color='r', ha='left', va='top', fontsize=16)
        ax.text(18.50, 0.55, f'Mean-offset: {offest_mean}    Median-offset: {offest_median}    Stdev-offset: {offest_stdev}', color='black', ha='center', va='center', fontsize=16, bbox=dict(facecolor='grey', alpha=0.67))

        plt.savefig(f'{curdir}/{mjdbin}_{filterbin}.png', bbox_inches='tight', dpi=600)

        print(f'Saved figure: {curdir}/{mjdbin}_{filterbin}.png.')

    
    
        # Also write the more extreme offsets to a text file
        if len(outlier_objects) > 0:
          with open(f'{curdir}/{mjdbin}_{filterbin}_extreme_offsets.txt', 'w', newline='') as output_file:
            keys = outlier_objects[0].keys()
            dict_writer = csv.DictWriter(output_file, keys)
            dict_writer.writeheader()
            dict_writer.writerows(outlier_objects)
          print(f'Recorded outlier offsets: {curdir}/{mjdbin}_{filterbin}_extreme_offsets.txt.')

    print(f'\nFinished!\n')