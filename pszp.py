# Minimal demonstration of using AstroPy specutils to read a plot reduced
# spectra from Liverpool Telescope SPRAT level 2 reduced FITS files
Version='1.2'

import sys
import os
import argparse
import numpy as np
import pandas as pd
import requests
import astropy.units as u
import matplotlib.pyplot as plt
import csv
from astropy.io import fits
from astropy.coordinates import SkyCoord
from alive_progress import alive_bar

plt.rcParams['font.family']='Times New Roman'
plt.rcParams['figure.figsize']=8,6
plt.rcParams['figure.autolayout']=True
plt.rcParams['mathtext.fontset']='dejavuserif' 

def PS1catalog(ra,dec,path):
    
    # radius is in degrees
    # limiting to no more than 10,000 catalog sources with >25 detections
    queryurl = 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/stack.json?'
    queryurl += 'ra='+str(ra)
    queryurl += '&dec='+str(dec)
    queryurl += '&radius=0.164'
    queryurl += '&columns=[raStack,decStack,gPSFMag,rPSFMag,iPSFMag,zPSFMag,yPSFMag,rKronMag]'
    queryurl += '&nDetections.gte=26&pagesize=10000'

    print('\nQuerying PS1 DR2 Stack for reference stars...\n')
    #print(queryurl)

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

        # Below is a bit of a hack to remove duplicate star records
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

        np.savetxt(path+'/ref_stars.dat',unique_star_data,fmt='%.8f\t%.8f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f', header='ra\tdec\tg\tr\ti\tz\ty', comments='')
        print('Success! Reference star file created: ref_stars.dat\n')

    else:
        sys.exit('Field not in PS1 DR2! Exiting.\n')



def binFilters(epochs):

  binned_filtered_epochs = {}
  binning_filters = ['g','r','i','z','y','w']

  for filterband in binning_filters:

    filterband_mjds_list = []
    filterband_file_list = []
    for epoch in epochs:
       if epoch[1] == filterband:
          filterband_mjds_list.append(epoch[0])
          filterband_file_list.append(epoch[2])

    if len(filterband_file_list) > 0:
      mjdbin = round(np.nanmean(filterband_mjds_list),3)
      filebin = {filterband:filterband_file_list}
      binned_filtered_epochs[mjdbin] = filebin

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
            # Where the calibration mag value is outside of the 15.5 - 22.5 range (Too bright and we will see saturation effects, too faint and we will see noise)
            if not np.isnan(object['PSF.CAL_PSF_MAG']):
              if 15.5 <= object['PSF.CAL_PSF_MAG'] <= 22.5:
                skycells_objects_list.append(object)

      epochs_skycell_objects[mjdbin][filterbin] = skycells_objects_list

  return epochs_skycell_objects



def calSkycellOffsets(epochs_skycell_objects, refstars):

  epochs_extreme_offsets = {}
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
      with alive_bar(len(skycellobjects), title=f'Calibrating {mjdbin} {filterbin}-band: ') as bar:   
        for object in skycellobjects:
          for star in refstars[f]:
            # Perform offset check if skycell object is within 0.1 arcsec of catalog reference star
            match_offset_deg = ((star[0] - object['PSF.RA'])**2 + (star[1] - object['PSF.DEC'])**2)**0.5
            if (match_offset_deg*3600 <= 0.1):
              object['CAL.MAG_TRUE'] = star[2]
              object['CAL.MAG_OFFSET'] = object['PSF.CAL_PSF_MAG'] - star[2]
              calibrated_objects_list.append(object)
              if abs(object['CAL.MAG_OFFSET']) > 0.5:
                extreme_offsets_outliers_list.append(object)
              # Once a calibration has been performed move onto next skycell object.
              break
          bar()
      mag_offsets_list = [x['CAL.MAG_OFFSET'] for x in calibrated_objects_list]
      # remove calbirated objects from the main list if they lie outside 3 stdev from the mean as these are typically galaxies that didn't get removed in the Star-galaxy separation check 
      # record these extreme outliers in a text file for manual checking later
      mean_clip = np.nanmean(mag_offsets_list)
      sigma_clip = 3.03 * np.nanstd(mag_offsets_list)
      calibrated_clipped_objects_list = [object for object in calibrated_objects_list if (mean_clip - sigma_clip <= object['CAL.MAG_OFFSET'] <= mean_clip + sigma_clip)]
      epochs_skycell_objects[mjdbin][filterbin] = calibrated_clipped_objects_list
      extreme_offsets_key = str(mjdbin) + '_' + str(filterbin) +'_extreme_offsets'
      epochs_extreme_offsets[extreme_offsets_key] = extreme_offsets_outliers_list

  return epochs_skycell_objects, epochs_extreme_offsets

def writeFileExtremeOffsets(outlier_objects):
  print(' ')
  if len(outlier_objects) > 0:
    for dictkey, dictrow in outlier_objects.items():
      if len(dictrow) > 0:
        with open(f'{curdir}/{dictkey}.txt', 'w', newline='') as output_file:
          keys = dictrow[0].keys()
          dict_writer = csv.DictWriter(output_file, keys)
          dict_writer.writeheader()
          dict_writer.writerows(dictrow)
        print(f'Recorded outlier offsets: {curdir}/{dictkey}.txt.')
  
  return

def writeFileAllOffsets(offset_collection, path):
  if len(offset_collection) > 0:
    for dictkey, dictrow in offset_collection.items():
      if len(dictrow) > 0:
        if not os.path.exists(path + '/offsets.txt'):
          with open(f'{path}/offsets.txt', 'w', newline='') as output_file:
            keys = dictrow[0].keys()
            dict_writer = csv.DictWriter(output_file, keys)
            dict_writer.writeheader()
            dict_writer.writerows(dictrow)
            print(' ')
            print(f'Recorded final offsets: {path}/offsets.txt.')
        else:
          with open(f'{path}/offsets.txt', 'a', newline='') as output_file:
            keys = dictrow[0].keys()
            dict_writer = csv.DictWriter(output_file, keys)
            dict_writer.writerows(dictrow)
            print(' ')
            print(f'Appended final offsets: {path}/offsets.txt.')
  return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Pan-STARRS Zeropoint Correction Tool')

    parser.add_argument('--directory', '-d', dest='directory', help='Object directory inside Bundles containing the .cmf files. Mandatory.', type=str)
    parser.add_argument('--coords', '-c', dest='coords', help='RA and DEC coordinates of object in degrees. Mandatory.', nargs=2, type=float)

    args = parser.parse_args()

    print(' ')
    print('=======================================================')
    print(' Welcome to PSZP: Pan-STARRS Zeropoint Correction Tool ')  
    print(f'                        V{Version}                        ')
    print('         Written by Michael Fulton (2021-2022)         ')        
    print('=======================================================')
    print(' ')


    # Set up the working directory
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
    if coords[1] > -30:
        print(f'Coordinates parsed: RA={coords[0]}, DEC={coords[1]}.\n')
    else:
        sys.exit('Object coordinates invalid! Please supply a DEC greater than -30 degrees.\n')
      

    # Create catalog of Pan-STARRS reference stars using object coordinates
    # Save catalog in parent (backpath) directory, or the child (currpath) directory if the parent directory is "Bundles"
    currpath = os.getcwd()
    os.chdir('..')
    backpath = os.getcwd()
    os.chdir(currpath)
    backpath_dirname = os.path.basename(backpath)
    currpath_dirname = backpath_dirname + '/' + os.path.basename(currpath)
    if backpath_dirname != 'Bundles':
      if not os.path.exists(backpath + '/ref_stars.dat'):
        PS1catalog(coords[0],coords[1],backpath)
        refcat = pd.read_csv(backpath + "/ref_stars.dat", sep="\t") 
      else:
        print(f'\nReference star catalog in dir={backpath_dirname} already exists. If you wish to make a new one, delete the current "ref_stars.dat" file and rerun the command.\n')
        refcat = pd.read_csv(backpath + "/ref_stars.dat", sep="\t")
    else:
      if not os.path.exists(currpath + '/ref_stars.dat'):
        PS1catalog(coords[0],coords[1],currpath)
        refcat = pd.read_csv(currpath + "/ref_stars.dat", sep="\t") 
      else:
        print(f'\nReference star catalog in dir={currpath_dirname} already exists. If you wish to make a new one, delete the current "ref_stars.dat" file and rerun the command.\n')
        refcat = pd.read_csv(currpath + "/ref_stars.dat", sep="\t")


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
      # w-band conversion taken from Pan-STARRS colour transformations in Table 6. of https://iopscience.iop.org/article/10.1088/0004-637X/750/2/99
      if not np.isnan(row['g']) and not np.isnan(row['r']):
        psw = round(0.04*row['g'] + 0.96*row['r'] - 0.015, 3)
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
    epochs_filter_binned = binFilters(epochs)


    # Load in the skycell objects for each of the epochs
    epochs_skycell_objects = loadSkycellObjects(epochs_filter_binned)
    print(f'Skycell objects loaded. Beginning calibration(s)...\n')


    # Calculate magnitude offsets between skycell objects and reference stars
    offset_objects, outlier_objects = calSkycellOffsets(epochs_skycell_objects, [refcat_g, refcat_r, refcat_i, refcat_z, refcat_y, refcat_w])
    print(f'\nSkycell offsets calibrated successfully.\n')


    # Plot offsets
    offset_collection = {}
    print(f'\nPlotting offsets...\n')
    for mjdbin, filterlist in offset_objects.items():
      for filterbin, skycellobjects in filterlist.items():
        
        # Plot lightcurve comparisons
        fig, ax = plt.subplots()
        fig.set_size_inches(11.7, 8.3, forward=True)
        plt.tight_layout()

        # split out the data into xy axes
        x = [object['CAL.MAG_TRUE'] for object in skycellobjects]
        y = [object['CAL.MAG_OFFSET'] for object in skycellobjects]

        # Plot data as a scatter
        ax.scatter(x, y, s=100, marker='x', color='black', zorder=1)

        # calculate statistics for the offsets
        offest_mean = round(np.nanmean(y),4)
        offest_median = round(np.nanmedian(y),4)
        offest_stdev = round(np.nanstd(y),4)
        
        # Set axes boundaries
        ax.set_xlim([15.0, 22.0])
        ax.set_xticks(np.arange(16, 21.5, 1))
        y_upper_axis_lim = round(offest_mean+(4*offest_stdev), 1)
        y_lower_axis_lim = round(offest_mean-(4*offest_stdev), 1)
        y_axis_range = y_upper_axis_lim-y_lower_axis_lim
        y_axis_ticksize = y_axis_range/5
        ax.set_ylim([y_lower_axis_lim, y_upper_axis_lim + (y_axis_ticksize/4)])
        ax.set_yticks(np.arange(y_lower_axis_lim + (y_axis_ticksize/2), y_upper_axis_lim + (y_axis_ticksize/2), y_axis_ticksize))

        ax.set_xlabel('Reference Star Magnitude', fontsize=26)
        ax.set_ylabel('Observed Magnitude Offset [Target - Reference]', fontsize=26)
        ax.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=24)
        ax.set_title(f'Difference Magnitude Offsets for: MJD {mjdbin} {filterbin}-band', fontsize=28)
        
        ax.text(15.05, +(y_axis_ticksize/3), 'Observed\nFainter', color='r', ha='left', va='center', fontsize=22)
        ax.axhline(linewidth=3, color='r')
        ax.text(15.05, -(y_axis_ticksize/3), 'Observed\nBrighter', color='r', ha='left', va='center', fontsize=22)
        ax.text(18.50, y_upper_axis_lim, f'Mean-offset: {offest_mean}    Median-offset: {offest_median}    Stdev-offset: {offest_stdev}', color='black', ha='center', va='top', fontsize=22, bbox=dict(facecolor='grey', alpha=0.67))

        plt.savefig(f'{curdir}/{mjdbin}_{filterbin}.png', bbox_inches='tight', dpi=600)
        print(f'Saved figure: {curdir}/{mjdbin}_{filterbin}.png.')

        #write offset analysis to offset collection object
        object_id = 'skycell' + '_' + str(mjdbin) + '_' + str(filterbin)
        offset_collection[object_id] = [{
          'subdir' : currpath_dirname,
          'mjd' : mjdbin,
          'filter' : filterbin,
          'offset_mean' : offest_mean,
          'offset_median' : offest_median,
          'offset_stdev' : offest_stdev
        }]

    
    
    # Also write the more extreme offsets to a text file
    writeFileExtremeOffsets(outlier_objects)

    # Also write the more extreme offsets to a text file
    if backpath_dirname != 'Bundles':
      writeFileAllOffsets(offset_collection, backpath)
    else:
      writeFileAllOffsets(offset_collection, currpath)

    print(f'\n\nFinished!\n')
