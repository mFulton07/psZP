from astropy.io import fits
from alive_progress import alive_bar
import matplotlib.pyplot as plt
import numpy as np
import math as m
import os
import re
import csv


def main (bundle_dir):
   
    # Create separate lists of the skycell fits files and reference stars dat files
    # By iterating through the bundle and picking out ".cmf" files for the skycells and ".dat" files for the reference stars
    full_bundle_dir = 'Bundles/'+bundle_dir
    all_files = [files for files in os.listdir('{bundle}'.format(bundle=full_bundle_dir) ) if os.path.isfile(os.path.join('{bundle}'.format(bundle=full_bundle_dir), files) ) ]
    skycells_list = [name for name in all_files if '.cmf' in name]
    refstars_list = [name for name in all_files if '.dat' in name]


    # Create an iterable list of the objects in the skycells
    # By opening up each skycell file one at a time
    # And extracting the relevant info for each psf object
    print('Loading Skycell Objects...')
    skycells_objects_list = []
    mjds_list = []
    cwd = os.getcwd()
    for skycell in skycells_list:
        
        # Extract the header and table data from the skycell file
        fits_image_filename = '{dir}/{bundle}/{file}'.format(dir=cwd, bundle=full_bundle_dir, file=skycell)
        hdul = fits.open(fits_image_filename)
        hdr = hdul[0].header
        data = hdul[1].data
        mjds_list.append(round(hdr['MJD-OBS'],5))
        filterband = hdr['HIERARCH FPA.FILTER'].strip('.012')

        # Create a simplifed record of the psf object
        # And add the object to an internal list containing all skycell objects
        for record in range(len(data)):

            object = {

                'PSF.RA':data['RA_PSF'][record],
                'PSF.DEC':data['DEC_PSF'][record],
                'PSF.CAL_PSF_MAG':float(data['CAL_PSF_MAG'][record]),
                'PSF.INST_PSF_MAG_SIG':float(data['PSF_INST_MAG_SIG'][record]),
                'FPA.MJD':round(hdr['MJD-OBS'],5),
                'FPA.AIRMASS':hdr['AIRMASS'],
                'FPA.FILTER':hdr['HIERARCH FPA.FILTER'].strip('.0'),
                'FPA.EXPTIME':hdr['EXPTIME'],
                'FPA.ZP':hdr['HIERARCH FPA.ZP']
            }
            # Rejecting PSFs...
            # Where the calibration failed
            # Where the calibration mag error is less than 5-sigma (mag_err <= 0.2)
            # Where the calibration mag value is outside of the 15 - 20 range (Too bright and we will see saturation effects, too faint and we will see noise)
            if not m.isnan(object['PSF.CAL_PSF_MAG']):
                if object['PSF.INST_PSF_MAG_SIG'] <= 0.2:
                    if 15 < object['PSF.CAL_PSF_MAG'] < 20:
                        skycells_objects_list.append(object)


    # Create an iterable list of the stars in the reference stars
    # By opening up each reference stars file one at a time
    # And extracting the relevant info for each reference star
    print('Loading Reference Stars...')
    refstars_stars_list = []
    cwd = os.getcwd()
    for refstar in refstars_list:
        
        # Extract the star data (spatial subset in a plain text file with columns [RA DEC MAG MAGERR C1_MAG C2_MAG])
        star_filter = re.search('v0_(.*).dat', refstar).group(1)
        fits_image_filename = '{dir}/{bundle}/{file}'.format(dir=cwd, bundle=full_bundle_dir, file=refstar)
        data = np.loadtxt(fits_image_filename)

        # Create a simplifed record of the reference star
        # And add the reference star to an internal list containing all reference stars
        for record in range(len(data)):

            star_rec = data[record]  
            star = {
                'STAR.RA':star_rec[0],
                'STAR.DEC':star_rec[1],
                'STAR.MAG':star_rec[2],
                'STAR.MAG_ERR':star_rec[3],
                'STAR.FILTER':star_filter
            }
            # Rejecting reference stars...
            # Where the reference mag error is less than 5-sigma (mag_err <= 0.2)
            if star['STAR.MAG_ERR'] <= 0.2:
                    refstars_stars_list.append(star)
    

    # Check if the skycell objects match the position of any of the reference stars (match accuracy is at 0.104 arcsec)
    # Then calculate the MAG_OFFSET between the CAL_PSF_MAG of the objects and MAG of the stars for reference stars with an apparent magnitude between
    # And add those that matched to a new internal crossmatched skycell objects list
    print('Calculating zeropoint offset between skycell objects and reference stars:')
    arcsec_to_degree = 0.0002777777778
    match_acc_arcsec = 0.105
    full_crossmatch_objects_list = []
    outlier_crossmatch_objects_list = []
    with alive_bar(len(skycells_objects_list)) as bar:   
        for object in skycells_objects_list:
            for star in refstars_stars_list:
                match_offset_deg = ((star['STAR.RA'] - object['PSF.RA'])**2 + (star['STAR.DEC'] - object['PSF.DEC'])**2)**0.5
                match_offset_arcsec = (match_offset_deg / arcsec_to_degree)
                if (match_offset_arcsec <= match_acc_arcsec):
                    object['CAL.MAG_OFFSET'] = object['PSF.CAL_PSF_MAG'] - star['STAR.MAG']
                    object['CAL.MAG_TRUE'] = star['STAR.MAG']
                    
                    # If the offset is extreme, also add the object to another internal list for printing out later
                    full_crossmatch_objects_list.append(object)
                    if abs(object['CAL.MAG_OFFSET']) >= 0.5:
                        outlier_crossmatch_objects_list.append(object)
            bar()
    

    # Sigma clip the offsets, removing any that greater than 2-stdev away from the mean
    offsets_list = [offset['CAL.MAG_OFFSET'] for offset in full_crossmatch_objects_list]
    Sigma_clip = 2 * np.nanstd(offsets_list)
    print('Offset limit used in sigma clipping: {clipping}...'.format(clipping=round(Sigma_clip,3)))
    crossmatch_objects_list = [offset for offset in full_crossmatch_objects_list if abs(offset['CAL.MAG_OFFSET']) <= Sigma_clip]


    # Create a list of unqiue mjds for grouping crossmatched skycell objects, keyed using the first mjd recorded for each bin
    mjds_list.sort(key=float, reverse=False)
    binsize_days = 0.00001
    for mjd1 in mjds_list:
        for mjd2 in reversed(mjds_list):
            if 0.0 < abs(mjd1 - mjd2) < binsize_days:
                mjds_list.remove(mjd2)


    # Group the crossmatched skycell objects on mjd
    nightly_crossmatches = {}
    for bin_mjd in mjds_list:
        temp_list = []
        for cal_object in crossmatch_objects_list:
            if abs(cal_object['FPA.MJD'] - bin_mjd) < binsize_days:
                temp_list.append(cal_object)
        nightly_crossmatches[bin_mjd] = temp_list


    # Create the "Outputs" subdirectory if it does not currently exist
    sub_directory = 'Outputs/op_'+bundle_dir
    directory_exists = os.path.exists(sub_directory)
    if not directory_exists:
        os.makedirs(sub_directory)


    # Make a scatter plot of MAG_OFFSET of the objects versus REFCAT_MAG of the reference stars
    print('Creating mag-offsets plots and text files...')
    for bin_mjd, crossmatch_objects in nightly_crossmatches.items():
        fig, ax = plt.subplots()
        x = [offset['CAL.MAG_TRUE'] for offset in crossmatch_objects]
        y = [offset['CAL.MAG_OFFSET'] for offset in crossmatch_objects]
        t = np.unique([offset['FPA.MJD'] for offset in crossmatch_objects])

        offest_mean = round(np.nanmean(y),3)
        offest_median = round(np.nanmedian(y),3)
        offest_stdev = round(np.nanstd(y),3)
        offset_time = round(np.nanmean(t),5)

        ax.scatter(x, y, s=150, marker='x', color='black', zorder=10)
        
        ax.text(14.55, 0.03, 'Observed\nFainter', color='r', ha='left', va='bottom', fontsize=16)
        ax.axhline(linewidth=3, color='r')
        ax.text(14.55, -0.03, 'Observed\nBrighter', color='r', ha='left', va='top', fontsize=16)

        ax.text(17.625, 0.45, 'Mean-offset: {mean}    Median-offset: {median}    Stdev-offset: {stdev}'.format(mean=offest_mean, median=offest_median, stdev=offest_stdev), color='black', ha='center', va='center', fontsize=16, bbox=dict(facecolor='grey', alpha=0.5))

        ax.set_title('Difference Image Offsets for {filter}-band on MJD={time}'.format(filter=filterband, time=offset_time), fontsize=36)
        ax.set_xlabel('REFcat Star AB Mag', fontsize=28)
        ax.set_ylabel('Mag Offset [Skycell - REFcat Star]', fontsize=28)

        ax.set_ylim([-0.5,0.5])
        ax.set_xlim([14.50,20.75])
        ax.tick_params(axis='both', length=6, width=3, labelsize=24)

        fig.set_size_inches(20, 10)
        plt.savefig('{path}/{filter}-{time}.png'.format(path=sub_directory, filter=filterband, time=offset_time))

    
    # Write the more extreme offsets to a text file
    with open('{path}/{filter}_extreme_offsets.txt'.format(path=sub_directory, filter=filterband), 'w', newline='') as output_file:
        keys = outlier_crossmatch_objects_list[0].keys()
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(outlier_crossmatch_objects_list)
    
    
    # Script finished
    print('Finished!')
    print('Please check "{output}" for your final offsets.'.format(output=sub_directory))
    
    return 0





print('#====================================================================================#')
print('      Welcome to the Zeropoint Corrections Tool for Pan-STARRS Difference Images      ')
print(' ')
mode = False
while mode is not True:
    bundle_dir = input('Please enter the name of your data bundle: ')
    
    if os.path.isdir('Bundles/'+bundle_dir) == True:
        mode = True
        main(bundle_dir)   
    else:
        print('Data bundle not found in the "Bundles" directory. Please try again - input is case sensitive...')
print(' ')
print('#====================================================================================#')