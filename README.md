# FIRES
FIRES (Fits Information Revision and Extraction System) is designed to run through archival and new HERCULES observations to determine whether the FITS header includes correct information, and adds some new useful information.

## Status
FIRES is currently being maintained and hosted up to date here. It will eventually be moved to another github, and this status will be updated to reflect that.
## Installation
Download the .zip and extract to your desired location. You can also clone the git by running
```
cd "your_main_folder"
git clone https://eng-git.canterbury.ac.nz/mwega/fires.git
```
This will make a new folder inside "your_main_folder", containing FIRES.

FIRES was tested in python 3.12.x

OPTIONAL: To set up a new python environment, make sure you have an anaconda installation (anaconda, miniconda, etc.), then run `conda create -n my_env_name python=3.12`
Ensure you then pip install atleast the following:

- numpy
- scipy
- pandas
- astropy
- astroquery (SIMBAD, if you don't want to install all of astroquery)
- os
- glob

Finally, run `conda activate my_env_name` to make sure your new python environment is being used.

## Usage

`fires_run.py` contains the code to run the program. Other files will not do anything if run.

## Functions

