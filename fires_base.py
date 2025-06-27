import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.coordinates import Angle
import astropy.units as u
import os
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord


def open_fits(fitsfile):
    # OPEN HERCULES FIT(S) FILE, return the header and image data
    with fits.open(fitsfile) as hdulist:
        data = hdulist[0].data
        header = hdulist[0].header
    return header, data

def guess_obs_with_exptime(header):
    # This returns 5 different guesses for the type of observation, based on the HERCEXPT header and the exposure time.
    # If the HERCEXPT header is present, it will return the type of observation based on the exposure time: White Light, Thorium Light, Stellar, or Stellar short (if the exposure time is suspiciously short for a stellar observation).
    # If all checks fail, it will return 'Unknown'.
    obs_type_guess = 'Unknown'  # Default guess
    
    exptime = header['EXPTIME']
    if not isinstance(exptime, float):
        print(f'EXPTIME is not a number: {exptime}, type: {type(exptime)}')
        obs_type_guess = 'Failed'  # If EXPTIME is not a number, return 'Unknown'
        return obs_type_guess  # Return 'Unknown' if EXPTIME is not a number
    if exptime < 2.5 and exptime > 0:
        obs_type_guess = ' White Light'
    elif exptime < 12.5 and exptime >= 2.5:
        obs_type_guess = ' Thorium Light'
    elif exptime >= 25:
        obs_type_guess = ' Stellar'
    elif exptime < 25 and exptime >= 12.5:
        # Stellars shouldnt be less than 30 seconds due to camera shutter, but physically it is possible
        obs_type_guess = ' Stellar short'
    else:
        # if somehow all these checks fail, return 'Unknown'
        obs_type_guess = ' Unknown'
    return obs_type_guess

def determine_obs_type(header):
    thorium_names = ['thar', 'thorium', 'thar lamp', 'throium', 'thorium lamp', 'arc']
    failed = False
    match = True
    # Determine the type of observation based on the HERCEXPT header and exposure time.
    exp_type = header['HERCEXPT'] # HERC Exposure Type
    obj_name = header['OBJECT'] # Object Name
    obs_type_guess = guess_obs_with_exptime(header) #Use EXPTIME Exposure TIme
    # Check if exp_type and guess match
    if exp_type.lower().strip() == obs_type_guess.lower().strip():
        obs_type = exp_type
    else:
        print(f'Observation type guess {obs_type_guess} does not match HERCEXPT {exp_type} for {obj_name}.')
        # Do more looking
        # Check object name
        if 'hd_' in obj_name.lower() and obs_type_guess == ' Stellar':
            obs_type = ' Stellar'
            match = False
        elif obj_name.lower() in thorium_names and obs_type_guess == ' Thorium Light':
            obs_type = ' Thorium Light'
            match = False
        else:
            # If no guesses match with normal expected values, return the guess alongside note checks failed
            failed = True
            obs_type = obs_type_guess
    if failed:
        failed = 'Failed'
        print(f'Observation type guess failed for {obj_name} with exposure time {header["EXPTIME"]}. Guess: {obs_type_guess}, HERCEXPT: {exp_type}')
        return obs_type, failed
    return obs_type, match

# gets the saved coordinates from the comment card in the header (only if HERCEXPT is stellar)
def get_coords(header, deg=False):
    assert header['HERCEXPT']==' Stellar', 'Not a stellar observation'
    assert deg in [True, False], 'deg must be True or False'
    if not deg:
        comment = header['COMMENT']
        # get the right comment
        ra = comment[0]
        dec = comment[1]
        return ra, dec
    elif deg:
        ra_deg = header['POSTN-RA']
        dec_deg = header['POSTN-DE']
        return ra_deg, dec_deg

def convert_deg_coords(ra_deg, dec_deg):
    # Convert RA from degrees to HMS (RA is measured in hours)
    ra_hms = Angle(ra_deg * u.deg).to_string(unit=u.hour, sep=':', precision=2)
    # Convert DEC from degrees to DMS
    dec_dms = Angle(dec_deg * u.deg).to_string(unit=u.deg, sep=':', alwayssign=True, precision=2)
    return ra_hms, dec_dms

def add_new_coords_to_header(header):
    # Set up Simbad query
    simbad = Simbad()
    simbad.TIMEOUT = 500  # Increase timeout for Simbad queries
    simbad.add_votable_fields('ra', 'dec')
    # Add a key-value pair to the header
    star_name = header['OBJECT']
    ra, dec = get_coords(header, deg=True)
    ra_hms, dec_dms = get_coords(header)  # Get the coordinates in HMS and DMS format
    # Add RA and DEC to the header in HMS and DMS format not in comment...
    header['RA'] = ra_hms
    header['DEC'] = dec_dms

    # Go to SIMBAD with the star name and get the coordinates
    result = simbad.query_object(star_name)
    if result is None: # Check if the request was successful
        print(f'No results found for {star_name} in SIMBAD.')
        return header # change nothing if request failed
    
    ra_simbad = result[0]['ra'] # extract RA and DEC from the result
    dec_simbad = result[0]['dec']

    # Check if the coordinates are similar to the ones in the header
    # Assume within ~ 2 arcminute is close enough
    seperation = 'Small'
    header_coords = SkyCoord(ra, dec, unit=(u.deg, u.deg))
    simbad_coords = SkyCoord(ra_simbad, dec_simbad, unit=(u.deg, u.deg)) 

    sep = header_coords.separation(simbad_coords)  # Calculate separation btween header coordinates and SIMBAD coordinates
    print(sep)

    if sep.arcmin > 2:
        seperation = 'Large'  # If the separation is greater than 2 arcminute, we consider it not similar
        # Check to see if the simbad coordinates are reasonable for the 1m telescope
        if dec_simbad > 15:
            separation = 'Unreasonable'
            print(f'Warning: SIMBAD coordinates for {star_name} are not reasonable for the 1m telescope. Difference from header: {sep.deg:.2f} degrees')
            return header, separation

    ra_simbad_hms, dec_simbad_dms = convert_deg_coords(ra_simbad, dec_simbad) # Convert to HMS and DMS format (: delimited)

    # Add the coordinates to the header
    header['SMBD_RA'] = ra_simbad_hms
    header['SMBD_DEC'] = dec_simbad_dms
    # print(ra_simbad, dec_simbad)
    return header, seperation

def save_new_fits(header, data, folder=None):
    # Set up required fields for the new filename
    obs_date = header['DATE'].strip().replace('-', '')  # Format the date for the filename
    target_name = header['OBJECT'].strip().replace(' ', '_').replace('-','_')  # Format the target name for the filename
    exposure_time = header['EXPTIME']
    if float(exposure_time).is_integer():
        exposure_time = float(exposure_time)  # Convert to float if it is an integer

    exposure_time = f"{header['EXPTIME']:.1f}".replace('.', 'p').zfill(4)  # Format the exposure time for the filename (1dp, padded to 4 characters if needed), p representing the decimal point
    if exposure_time[-1] == '0' and exposure_time[-2]=='p':  # If the last two characters are 0 and p, remove them (no need to keep decimal point if it is a whole number)
        exposure_time = exposure_time[:-2]
    if len(exposure_time) > 4:
        exposure_time = exposure_time[:4]
        if exposure_time[-1] == 'p': # if we truncate and leave the p, no need to keep it
            exposure_time = exposure_time[:-1]

    # Save the new header and data to a new fits file
    if folder is None: # default to subfolder in folder of this script
        current_folder = os.path.dirname(os.path.abspath(__file__))
        output_folder = os.path.join(current_folder, 'checked_files')
    else:
        output_folder = os.path.join(folder, 'checked_files')  # Use provided folder, with a subfolder for new files
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    running_number = 1
    # Check number of files in the output folder with the same objectname to determine the running number
    existing_files = [f for f in os.listdir(output_folder) if f.startswith(f'H{obs_date}-{target_name}-')]
    if existing_files:
        running_number += len(existing_files)  # Increment running number based on existing files
    running_number = str(running_number).zfill(3)  # Pad the running number to 3 characters

    filename_structure = f'H{obs_date}-{target_name}-{exposure_time}-{running_number}.fit' # Set up the filename structure as HYYYYMMDD-Targetname-Exptime-RNo.fit
    if len(filename_structure) > 30:
        print(f'Warning: Filename {filename_structure} is too long ({len(filename_structure)} characters...).')

    output_file = os.path.join(output_folder, filename_structure)
    
    # Create a new HDU with the updated header and data
    hdu = fits.PrimaryHDU(data=data, header=header)
    hdu.writeto(output_file, overwrite=True)
    
    print(f'Saved updated FITS file to {output_file}')