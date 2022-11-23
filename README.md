# psZP
### Please Read the [Zeropoint Corrections Wiki](https://psweb.mp.qub.ac.uk/psat-lv-wiki/index.php/Zeropoint_Corrections_for_Pan-STARRS_difference_images) before using psZP <br />

Welcome to the **Zeropoint Corrections Tool for Pan-STARRS Difference Images**. To get started follow the below steps.<br /><br />

## Step 1: Cloning the psZP repository
On a terminal, pull down the psZP repo from Github into your desired directory:
```
cd  #DESIRED_DIRECTORY#
git clone https://github.com/mFulton07/psZP.git
```
<br /><br />
## Step 2: Creating the Python environment
psZP requires some non-base python packages to operate.
Create a new conda environment:
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
statistics
csv
```
<br /><br />
## Step 3: Add your data bundles to the "psZP/Bundles" directory
Your bundles should be in their own sub-directory and contain the _.cmf skycell image_ files. The _.smf image_ files are not required as all the necessary info can be found within the _.cmf image_ files. The _.dat reference stars_ file is also not required as the script will automatically query PS1-DR2 Stack Archive for reference stars. Note that the bundles can be a mix of filters and MJDs.
<br /><br /><br />
## Step 4: Run psZP
On a terminal:
1. activate the pszp python environment...
2. cd to the psZP directory...
3. run the psZP startup file and pass in the name of your data bundle sub-directory along with the Pan-STARRS coordinates of your object...
Example:
```
conda activate pszp
cd  #DESIRED_DIRECTORY#/psZP/
python pszp.py -d MyBundle -c 123.456789 -12.3456789
...
```
Note that psZP can only accept one data bundle sub-directory at a time.
<br /><br /><br />
## Step 5: Check "psZP/Outputs" directory for your offsets
When finished, psZP will output:
* (Always) A scatter plot of the Mag Offset versus the True Mag for all the skycell images in a particualr filter on a particualr night. Plot includes labels for the mean, median and stdev of the offsets.
* (Sometimes) A comma deliminated text file of the skycell objects that had the most extreme offsets, if there are any outliers.
