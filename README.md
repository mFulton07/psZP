# psZP
### Please Read the [Zeropoint Corrections Wiki](https://psweb.mp.qub.ac.uk/psat-lv-wiki/index.php/Zeropoint_Corrections_for_Pan-STARRS_difference_images) before using psZP <br />

Welcome to the **Zeropoint Corrections Tool for Pan-STARRS Difference Images**. To get started, follow the below steps.<br /><br />

## Step 1: Cloning the psZP repository
On a terminal, pull down the psZP repo from GitHub into your desired directory:
```
cd  #DESIRED_DIRECTORY#
git clone https://github.com/mFulton07/psZP.git
```
<br /><br />
## Step 2: Creating the Python environment
psZP requires some non-base Python packages to operate.
Create a new conda environment. psZP was originally created in Python V3.7, but any Version 3.7+ should work:
```
conda create --name pszp python=3.7 pip
```
and install the following packages:
```
argparse
requests
alive_progress
astropy
numpy
pandas
matplotlib 
csv
```
<br /><br />
## Step 3: Add your data bundles to the "psZP/Bundles" directory
Your bundles should be in their own sub-directory and contain the _.cmf skycell image_ files. The _.smf image_ files are not required as all the necessary info can be found within the _.cmf image_ files. The _.dat reference stars_ file is no longer required as the code will automatically query PS1-DR2 Stack Archive for reference stars in a 10 arcmin (0.164 degrees) radius of the target position. Note that the bundles can be a mix of filters, but typically should not be a mix of MJDs (unless inter-night stacking is involved). Different nights should be partitioned into further sub-sub-directories, depending on the images that have been used to produce the stacked difference images or flux.
<br /><br /><br />
## Step 4: Run psZP
On a terminal:
1. activate the pszp python environment...
2. cd to the psZP directory...
3. run the psZP startup file and pass in the name of your data bundle sub-directory (or sub-sub-directory if required) along with the Pan-STARRS coordinates of your object...
Example:
```
conda activate pszp
cd  #DESIRED_DIRECTORY#/psZP/
python pszp.py -d MyBundle/AT2020abc/YYYY-MM-DD -c 123.456789 -12.3456789
...
```
Note that psZP can only accept one data bundle sub-directory at a time.
<br /><br /><br />
## Step 5: Check the sub- and sub-sub-directories for the results
When finished, psZP will output the following. The file locations will be specified in the terminal window:
* (Always) For each skycell file, A scatter plot of the Mag Offset versus the True Mag for all the skycell objects in a particular filter for a particular MJD. The figure includes labels for the mean, median and stdev of the offsets plotted.
* (Sometimes) For each skycell file, A comma-deliminated text file of the skycell objects that were extremely offset (>3 stdev from mean) and were rejected from the data before analysing and plotting.
* (Always) A consolidated, comma-deliminated text file containing all offset values determined for the data bundles. This will be saved at the parent directory level if applicable (same location as the reference stars file).
