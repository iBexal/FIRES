import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
import os

from fires_base import *
from order_tracing import *

# def __main__():
current_folder = os.path.dirname(os.path.abspath(__file__))
template_folder = os.path.join(current_folder, 'template_files')

files = ['test_data/correct/J0763016.fit', 'test_data/correct/J0763017.fit', 'test_data/correct/J0763018.fit', 'test_data/correct/J0763019.fit', 'test_data/correct/J0763028.fit', 'test_data/correct/J0763029.fit', 'test_data/correct/J0763042.fit']

checks = pd.DataFrame(columns=['File', 'HERCEXPT', 'Guess', 'Match', 'Seperation'])

for i, file in enumerate(files):
    header, data = open_fits(file)
    if '42' in file:
        header['HERCEXPT'] = 'White Light'  # Force a known type for testing
    obs_type, match = determine_obs_type(header)
    if match == 'Failed':
        # We need to go deeper
        data = data[:, 3300:3600] #for speed, plus allows for `len(tellurics) <= 3` to be a valid choice
        orders = order_tracing(header, data)

        fluxes = get_flux_from_orders(data, orders)
        norm_fluxes = cut_order_edge(norm_orders(fluxes))
        peaks, obs_type = find_tellurics_exp_type(norm_fluxes)
        match = 'order_traced'


    header['HERCEXPT'] = obs_type  # Update header with the determined observation type
    if 'stellar' in obs_type.lower():
        header = add_new_coords_to_header(header)
        coords = get_coords(header)
        print(coords, header['SMBD_RA'], header['SMBD_DEC'])

        checks.loc[i]  = [file,
        header['HERCEXPT'],
        obs_type,
        match
    ]
print(checks)
