# psZP
### Please Read the [Zeropoint Corrections Wiki](https://psweb.mp.qub.ac.uk/psat-lv-wiki/index.php/Zeropoint_Corrections_for_Pan-STARRS_difference_images) before using psZP <br />

Welcome to the **Zeropoint Corrections Tool for Pan-STARRS Difference Images**. To get started follow the below steps.<br /><br />

## Step 1: Cloning the psZP repository
On a terminal, pull down the psZP repo into your desired directory:
```
cd  #DESIRED_DIRECTORY#
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
## Step 3: Add your data bundles to the "psZP/Bundles" directory
Your bundles will be a directory containing the _.cmf & .smf skycell image files_ as well as the _.dat reference stars files_. It is important to note that these **bundles should be filter specific** (combining more than one filter in the bundles will result in errors) but can contain skycell images from multiple nights.
<br /><br /><br />
## Step 4: Run psZP
On a terminal:
1. activate the pszp python environment...
2. cd to the psZP directory...
3. run the psZP startup file...
4. follow the psZP prompts provided...
```
conda activate pszp
cd  #DESIRED_DIRECTORY/psZP#
python pszp_startup.py
...
```
Note that psZP can only accept one data bundle at a  time. If you have multiple bundles, simply run psZP multiple times (using a different bundle each time).
<br /><br /><br />
## Step 5: Check "psZP/Outputs" directory for your offsets
When finished, psZP will output per nighlty basis:
* A scatter plot of Skycell Object Mag Offset versus True AB Mag. Plot includes labels for the medain and mean offsets.
