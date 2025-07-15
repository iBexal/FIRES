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
**fires_base.py**
- open_fits(fitsfile):
This opens a given HERCULES fits file, and returns the header and data.
fitsfile is a string of a file location.

- guess_obs_with_exptime(header):
This takes a header, reads through the given exposure length and returns an educated guess on what type of file would be expected to have the given exposure time.
header must be a header returned by open_fits, or one read in by using astropy's fits.open().

- determine_obs_type(header, log=True):
- - This takes a header from open_fits, and checks the observation type using the observation type given by the header, a guess made by guess_obs_with_exptime(), and the object name given by the header.
- - set log=False to disable the check outputs, and only return the observation type
If the guess from the exposure time and the observation type given by the header are the same, the check is passed as a match. 
If the check fails, the guess and the object name are compared. If they match, the check is passed without a match.
If the check if passed, the observation type is returned alongside the check as either a "match", or not.
If all of these checks fail, the check is failed and the observation type is returned alongside the check as "failed".


- get_coords(header, deg=False):
- - This takes a header and returns the sexagesimal ra and dec saved in the comment card of the HERCULES header.
- - If deg=True, it will return the ra and dec in a degrees format.

- def convert_deg_coords(ra_deg, dec_deg):
- - This takes a set of ra and dec in degrees format and converts it into sexagesimal format. 
- - The ra and dec must be `float` or `int`.

- def add_new_coords_to_header(header, log=True):
- - This takes a header and adds new cards based on the sexigesimal coordinates saved in the comment of the header, and coordinates from SIMBAD. The updated header is returned.
- - Set log=False to turn off the log output, and return only the new header
If the exposure type is a `Stellar`, SIMBAD will be queried using the object name to generate coordinates. If these coordinates are sufficiently seperated from the coordinates saved in the header, this result is also returned alongside the new header. If the coordinates are further north than +15 degrees, the SIMBAD coordinates are noted as "Unreasonable".

- def save_new_fits(header, data, folder=None):
This takes a FITS header, respective data, and a folder location (as a `string`), and saves a new FITS file in a subfolder of `folder` with a given filename structure.
If `folder=None` (i.e. it is not defined), the subfolder is created in the same folder as `fires_base.py`.
The filename structure is as follows: `HYYYYMMDD-Targetname-Exptime-RNo.fit`; where YYYYMMDD is the date of observation from the header, Targetname is the object name from the header, Exptime is the exposure time from the header , and RNo being the running number of exposures for the given object name and given date.
The exposure time is saved as four digits long, with frontrunning 0s if required. `p` represents a decimal point if required, with only one decimal place being saved.
The running number is calculated by the number of files already saved in the subfolder, and is not respective of actual observation time or real numbers of exposures taken.

