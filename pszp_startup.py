from astropy.io import fits
from alive_progress import alive_bar
import matplotlib.pyplot as plt
import numpy as np
import statistics
import os
import re


def main (bundle_dir):
   
    # Create separate lists of the skycell fits files and reference stars fits files
    # By iterating through the bundle and picking out ".cmf" files for the skycells and ".dat" files for the reference stars
    full_bundle_dir = 'Bundles/'+bundle_dir
    all_files = [files for files in os.listdir('{bundle}'.format(bundle=full_bundle_dir) ) if os.path.isfile(os.path.join('{bundle}'.format(bundle=full_bundle_dir), files) ) ]
    skycells_list = [name for name in all_files if '.cmf' in name]
    refstars_list = [name for name in all_files if '.dat' in name]


    # Create an iterable list of the objects in the skycells
    # By opening up each skycell file one at a time
    # And extracting the revlant info for each object record
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
        filterband = hdr['HIERARCH FPA.FILTER'].strip('.0')

        # Create a simplifed version of the object record
        # And add the object to a list containing all skycell objects
        for record in range(len(data['CAL_PSF_MAG'])):

            object = {

                'PSF.RA':round(data['RA_PSF'][record],5),
                'PSF.DEC':round(data['DEC_PSF'][record],5),
                'PSF.CAL_PSF_MAG':float(data['CAL_PSF_MAG'][record]),
                'PSF.CAL_PSF_MAG_SIG':float(data['CAL_PSF_MAG_SIG'][record]),
                'FPA.MJD':round(hdr['MJD-OBS'],5),
                'FPA.BS_RA':hdr['HIERARCH FPA.RA'],
                'FPA.BS_DEC':hdr['HIERARCH FPA.DEC'],
                'FPA.AIRMASS':hdr['AIRMASS'],
                'FPA.FILTER':hdr['HIERARCH FPA.FILTER'].strip('.0'),
                'FPA.EXPTIME':hdr['EXPTIME'],
                'FPA.ZP':hdr['HIERARCH FPA.ZP']
            }
            skycells_objects_list.append(object)


    # Create an iterable list of the stars in the reference stars
    # By opening up each reference stars file one at a time
    # And extracting the revlant info for each star record
    print('Loading Reference Stars...')
    refstars_stars_list = []
    cwd = os.getcwd()
    for refstar in refstars_list:
        
        # Extract the star data (spatial subset in a plain text file with columns [RA DEC MAG MAGERR C1_MAG C2_MAG])
        star_filter = re.search('v0_(.*).dat', refstar).group(1)
        fits_image_filename = '{dir}/{bundle}/{file}'.format(dir=cwd, bundle=full_bundle_dir, file=refstar)
        data = np.loadtxt(fits_image_filename)

        # Create a simplifed version of the star record
        # And add the star to a list containing all reference stars
        for record in range(len(data)):

            star_rec = data[record]  
            star = {
                'STAR.RA':round(star_rec[0],5),
                'STAR.DEC':round(star_rec[1],5),
                'STAR.MAG':star_rec[2],
                'STAR.MAG_ERR':star_rec[3],
                'STAR.FILTER':star_filter
            }
            refstars_stars_list.append(star)
    

    # Check if the skycell objects are at the position of any of the reference stars (match accuracy is at 0.1 arcsec)
    # Then calculate the MAG_OFFSET between the CAL_PSF_MAG of the objects and MAG of the stars for stars between 16.5 and 20.5 mag (Too bright and we will see saturation effects. Too faint and we will see noise in the CAL_PSF_MAG).
    # And add those that matched to a new crossmatched skycell objects list
    print('Calculating zeropoint offset between skycell objects and reference stars:')
    crossmatch_objects_list = []
    with alive_bar(len(skycells_objects_list)) as bar:   
        for object in skycells_objects_list:
            for star in refstars_stars_list:
                if (object['PSF.RA'] == star['STAR.RA'] and object['PSF.DEC'] == star['STAR.DEC'] and 16.5 <= star['STAR.MAG'] <= 20.5):
                    object['CAL.MAG_OFFSET'] = object['PSF.CAL_PSF_MAG'] - star['STAR.MAG']
                    object['CAL.MAG_TRUE'] = star['STAR.MAG']
                    crossmatch_objects_list.append(object)
            bar()


    # Reduce list of mjds to unique binning mjds keyed using the first mjd recorded for each bin
    mjds_list.sort(key=float, reverse=False)
    binsize_days = 0.125
    for mjd1 in mjds_list:
        for mjd2 in reversed(mjds_list):
            if 0.0 < abs(mjd1 - mjd2) < binsize_days:
                mjds_list.remove(mjd2)


    # Bin the crossmatched skycell objects on mjd
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


    # For visiual representation
    # Make a scatter plot of MAG_OFFSET of the objects versus MAG of the stars (true magnitude)
    for bin_mjd, crossmatch_objects in nightly_crossmatches.items():
        fig, ax = plt.subplots()
        x = [offset['CAL.MAG_TRUE'] for offset in crossmatch_objects]
        y = [offset['CAL.MAG_OFFSET'] for offset in crossmatch_objects]

        offest_mean = round(statistics.mean(y),3)
        offest_median = round(statistics.median(y),3)
        offset_time = round(bin_mjd + (binsize_days/2.0),2)

        ax.scatter(x, y, s=150, marker='x', color='black', zorder=10)
        
        ax.text(16.05, 0.03, 'Observed\nFainter', color='r', ha='left', va='bottom', fontsize=16)
        ax.axhline(linewidth=3, color='r')
        ax.text(16.05, -0.04, 'Observed\nBrighter', color='r', ha='left', va='top', fontsize=16)

        ax.text(18.375, 0.75, 'Mean-offset: {mean}    Median-offset: {median}'.format(mean=offest_mean, median=offest_median), color='black', ha='center', va='center', fontsize=16, bbox=dict(facecolor='grey', alpha=0.5))

        ax.set_title('Difference Image Offsets for {filter}-band on MJD={time}'.format(filter=filterband, time=offset_time), fontsize=36)
        ax.set_xlabel('True AB Mag', fontsize=28)
        ax.set_ylabel('Skycell Object Mag Offset', fontsize=28)

        ax.set_ylim([-0.8,0.8])
        ax.set_xlim([16.00,20.75])
        ax.tick_params(axis='both', length=6, width=3, labelsize=24)

        fig.set_size_inches(20, 10)
        plt.savefig('{path}/wband_mjd{time}.png'.format(path=sub_directory, time=offset_time))

    print('Finished!')
    print('Please check "{output}" for your final offsets.'.format(output=sub_directory))
    print(' ')
    

    return 0





print('#====================================================================================#')
print('      Welcome to the Zeropoint Corrections Tool for Pan-STARRS Difference Images')
print(' ')
mode = False
while mode is not True:
    bundle_dir = input('Please enter the name of your data bundle: ')
    
    if os.path.isdir('Bundles/'+bundle_dir) == True:
        mode = True
        main(bundle_dir)   
    else:
        print('Data bundle not found in the "Bundles" directory. Please try again - input is case sensitive...')

