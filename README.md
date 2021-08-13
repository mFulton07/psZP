# psZP
### Please Read the [Zeropoint Corrections Wiki](https://psweb.mp.qub.ac.uk/psat-lv-wiki/index.php/Zeropoint_Corrections_for_Pan-STARRS_difference_images) before using psZP <br />

Welcome to **psZP: Zeropoint Corrections Tool for Pan-STARRS Difference Images**. To get started follow the below steps.<br /><br />

## Step 1: Cloning the psZP repository
On a terminal, pull down the psZP repo into your desired directory:
```
cd  #DIRECTORY_NAME#
git clone https://github.com/mFulton07/psZP.git
```
<br /><br />
## Step 2: Creating the Python environment
psZP requires some non-base python packages to operate.
Create a new conda environment:
```
conda create conda create --name pszp python=3.7
```
and install the following packages:
```
alive_progress
astropy
numpy
matplotlib 
statistics
```
<br /><br />
## Step 3: Add your nightly bundles to the psZP/Bundles directory
Your bundles will be a directory containing the _.cmf & .smf skycell image files_ as well as the _.dat reference stars files_. It is important to note that these **bundles should be filter specific** (combining more than one filter in the bundles will result in errors) but can contain images from multiple nights.
