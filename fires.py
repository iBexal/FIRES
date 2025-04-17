import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from skimage.metrics import structural_similarity as ssim
import os
import scipy as sp
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter, maximum_filter1d, median_filter
from copy import deepcopy
# from skimage.morphology import watershed
# from skimage.feature import peak_local_max

current_folder = os.path.dirname(os.path.abspath(__file__))
template_folder = os.path.join(current_folder, 'template_files')

def open_fits(fitsfile):
    hdulist = fits.open(fitsfile)
    data = hdulist[0].data
    header = hdulist[0].header
    hdulist.close()
    return header, data

def check_obs_type(header):
    # Open HERCEXPT card to check the type of observation the file thinks it is
    if header['HERCEXPT']==' White Light':
        return 'white'
    elif header['HERCEXPT']==' Thorium Light':
        return 'thorium'
    elif header['HERCEXPT']==' Stellar':
        return 'stellar'
    else:
        return None

# gets the saved coordinates from the comment card in the header (only if HERCEXPT is stellar)
def get_coords(header):
    assert header['HERCEXPT']==' Stellar', 'Not a stellar observation'
    comment = header['COMMENT']
    # get the right comment
    ra = comment[0]
    dec = comment[1]
    return ra, dec

def order_tracing(header, data):
    # get obs type
    obs_type = check_obs_type(header)
    if obs_type == 'thorium':
        print('Thorium, cannot order trace')
        return None
    
    # data = data[350:-350,500:]
    
    reference_column = int(data.shape[0] / 2) # Reference column for tracing (middle of the image)
    # spatial_profile = np.median(data, axis=0)  # Median collapse along y-axis
    spatial_profile = data[reference_column,:] # Collapse along y-axis

    # Detect peaks (orders)
    # Apply Gaussian smoothing to reduce noise and enhance weak orders
    smoothed_profile = gaussian_filter(spatial_profile, sigma=5)
    
    # sobel filter
    smooth_data = data
    smoothed_image = gaussian_filter(smooth_data, sigma=2)
    

    # Determine a dynamic prominence threshold (adjust if necessary)
    peak_prominence = (np.max(smoothed_profile) - np.min(smoothed_profile)) * 0.01  # 2% of peak range
    min_height = np.max(smoothed_profile) * 0.007  # Ignore peaks below 1% of the highest peak

    # Detect peaks with prominence filtering
    peaks, properties = find_peaks(
        smoothed_profile, 
        height=min_height,  # Ignore weak peaks
        distance=20, 
        prominence=peak_prominence  # Require peaks to stand out
    )

    order_traces = []

    # Trace each order **up and down** from detected peak
    for peak in peaks:
        x_positions = [peak] # peak is a given COLUMN
        y_positions = [reference_column]

        # Trace downwards
        fluxes = []
        for y in range(reference_column + 1, data.shape[0]):
            search_range = 10
            x_min = max(0, x_positions[-1] - search_range)
            x_max = min(data.shape[1] - 1, x_positions[-1] + search_range)
            
            row_slice = smoothed_image[y, x_min:x_max]
            if np.any(row_slice):
                new_x = np.argmax(row_slice) + x_min
            if abs(new_x - x_positions[-1]) <= search_range:  # Ensure continuity
                x_positions.append(new_x)
                y_positions.append(y)
            else:
                break
        
        # Reverse the lists to start tracing upwards from the original peak
        x_positions.reverse()
        y_positions.reverse()

        # Trace upwards
        for y in range(reference_column - 1, 0, -1):
            search_range = 10
            x_min = max(0, x_positions[-1] - search_range)
            x_max = min(data.shape[1] - 1, x_positions[-1] + search_range)
            
            row_slice = smoothed_image[y, x_min:x_max]
            if np.any(row_slice):  # Avoid empty slices
                new_x = np.argmax(row_slice) + x_min
            if abs(new_x - x_positions[-1]) <= search_range:  # Ensure continuity
                x_positions.append(new_x)
                y_positions.append(y)
            else:
                break

        # Reverse the lists back to their original order
        x_positions.reverse()
        y_positions.reverse()

        # Store the traced order
        order_traces.append((x_positions, y_positions))
    return order_traces

def get_flux_from_orders(data, order_traces):
    # get obs type
    # obs_type = check_obs_type(header)
    
    # Get the flux from the orders
    fluxes = []
    for order in order_traces:
        column_positions, row_positions = order
        flux = []
        # fit gaussian across each order per column and extract flux weighted by the gaussian
        
        for i in range(-3,3):
            flux.append(data[row_positions, np.array(column_positions)+i])
        # average the fluxes, weighted by a gaussian shape
        flux = np.array(flux)
        flux = np.mean(flux, axis=0)
        
        # reverse the order
        flux = flux[::-1]
        fluxes.append(flux)
    return fluxes


def norm_orders(order_fluxes):
    # normalise the fluxes
    norm_fluxes = []
    for i, fluxes in enumerate(order_fluxes):
        fluxes = fluxes / np.max(fluxes)
        
        order_new = deepcopy(fluxes)
        order_filter = deepcopy(fluxes)
        # median filter
        order_filter = median_filter(order_filter, size=10)
        # max filter
        order_filter = maximum_filter1d(order_filter, size=40, mode='constant')

        # smooth the order_filter
        order_filter = gaussian_filter(order_filter, sigma=15)
        order_new = order_new/order_filter
        
        norm_fluxes.append(order_new)
    return np.array(norm_fluxes)

def cut_order_edge(orders):
    # cut the edges of the orders
    new_orders = []
    for i, order in enumerate(orders):
        length = len(order)
        # take only middle 50% of order
        new_orders.append(order[int(length/4):int(length*3/4)])
    return new_orders

def find_tellurics(orders):
    # find the tellurics in the orders
    tellurics = []
    no_peaks = 0
    for i, order in enumerate(orders):
        # find the peaks
        # flip order so absorption lines are up
        order = 1-order
        peaks, _ = find_peaks(order, height=0.8)
        # append the peaks to the tellurics list
        tellurics.append(peaks)
        if len(peaks) > 0:
            no_peaks += 1
    
    obj_type_guess = 'stellar'
    if no_peaks == 0:
        obj_type_guess = 'white'
    return tellurics, obj_type_guess


# header, data = open_fits('test_data/correct/J0763019.fit')
files = ['16', '17', '18', '19', '28', '29', '42']
for i in files:
    header, data = open_fits(f'test_data/correct/J07630{i}.fit')

    data = data[:, 3300:3600]

    orders = order_tracing(header, data)
    if orders is None:
        continue
    fluxes = get_flux_from_orders(data, orders)
    norm_fluxes = cut_order_edge(norm_orders(fluxes))
    peaks = find_tellurics(norm_fluxes)


    print(peaks[1], header['HERCEXPT'])


plt.figure()


add=0
for i, order in enumerate(norm_fluxes):
    xs = np.array(range(len(order)))
    plt.plot(xs+add, order)
    add += xs[-1]

# plt.xlim(10000,12000)
# plt.figure()
# plt.imshow(data, cmap='gray', vmin=np.percentile(data, 1), vmax=np.percentile(data, 99))
# plt.figure()
# plt.imshow(white_data, cmap='gray', vmin=np.percentile(white_data, 1), vmax=np.percentile(white_data, 99))
plt.show()