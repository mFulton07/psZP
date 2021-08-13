#pip install astropy
#pip install alive_progress
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
    all_files = [files for files in os.listdir('Bundles/{bundle}'.format(bundle=bundle_dir) ) if os.path.isfile(os.path.join('Bundles/{bundle}'.format(bundle=bundle_dir), files) ) ]
    skycells_list = [name for name in all_files if '.cmf' in name]
    refstars_list = [name for name in all_files if '.dat' in name]


    # Create an iterable list of the objects in the skycells
    # By opening up each skycell file one at a time
    # And extracting the revlant info for each object record
    print('Loading Skycell Objects...')
    skycells_objects_list = []
    unique_mjds_list = []
    cwd = os.getcwd()
    for skycell in skycells_list:
        
        # Extract the header and table data from the skycell file
        fits_image_filename = '{dir}/{bundle}/{file}'.format(dir=cwd, bundle=bundle_dir, file=skycell)
        hdul = fits.open(fits_image_filename)
        hdr = hdul[0].header
        data = hdul[1].data
        unique_mjds_list.append(round(hdr['MJD-OBS'],5))
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
        fits_image_filename = '{dir}/{bundle}/{file}'.format(dir=cwd, bundle=bundle_dir, file=refstar)
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
    print('Calculating photometry offset between skycell objects and reference stars:')
    crossmatch_objects_list = []
    with alive_bar(len(skycells_objects_list)) as bar:   
        for object in skycells_objects_list:
            for star in refstars_stars_list:
                if (object['PSF.RA'] == star['STAR.RA'] and object['PSF.DEC'] == star['STAR.DEC'] and 16.5 <= star['STAR.MAG'] <= 20.5):
                    object['CAL.MAG_OFFSET'] = object['PSF.CAL_PSF_MAG'] - star['STAR.MAG']
                    object['CAL.MAG_TRUE'] = star['STAR.MAG']
                    crossmatch_objects_list.append(object)
            bar()
    print('Calculation Complete!')
    print('Difference image offsets are recorded below and in plot popup:')


    # For each of the crossmatched skycell objects
    # Bin them on unique mjds and calculate the MEAN_MAG_OFFSET in each bin
    # Then print the results to screen
    unique_mjds_list.sort(reverse=True)
    for unique_mjd in unique_mjds_list:
        temp_list = []
        for cal_object in crossmatch_objects_list:
            if cal_object['FPA.MJD'] == unique_mjd:
                temp_list.append(cal_object['CAL.MAG_OFFSET'])
        print(' ')
        print(
            'filter: '+filterband+'\n'+
            'image-mjd: '+str(unique_mjd)+'\n'+
            'mean-mag-offset: '+str(round(statistics.mean(temp_list),3))+'\n'+
            'median-mag-offset: '+str(round(statistics.median(temp_list),3))
        )
    

    # For visiual representation
    # Make a scatter plot of MAG_OFFSET of the objects versus MAG of the stars (true magnitude)
    fig, ax = plt.subplots()
    x = [offset['CAL.MAG_TRUE'] for offset in crossmatch_objects_list]
    y = [offset['CAL.MAG_OFFSET'] for offset in crossmatch_objects_list]

    ax.scatter(x, y, s=150, marker='x', color='black', zorder=10)
    
    ax.text(16.05, 0.03, 'Observed\nFainter', color='r', ha='left', va='bottom', fontsize=16)
    ax.axhline(linewidth=3, color='r')
    ax.text(16.05, -0.04, 'Observed\nBrighter', color='r', ha='left', va='top', fontsize=16)

    ax.text(18.375, 0.75, 'Mean-offset: {mean}    Median-offset: {median}'.format(mean=round(statistics.mean(y),3), median=round(statistics.median(y),3)), color='black', ha='center', va='center', fontsize=16, bbox=dict(facecolor='grey', alpha=0.5))

    ax.set_title('Difference Image Offsets for {filter}-band'.format(filter=filterband), fontsize=36)
    ax.set_xlabel('True AB Mag', fontsize=28)
    ax.set_ylabel('Skycell Object Mag Offset', fontsize=28)

    ax.set_ylim([-0.8,0.8])
    ax.set_xlim([16.00,20.75])
    ax.tick_params(axis='both', length=6, width=3, labelsize=24)

    plt.show()
    

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

