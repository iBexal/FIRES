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


def order_tracing(header, data, extract=False):
    # get obs type

    reference_column = int(data.shape[0] / 2) # Reference column for tracing (middle of the image)

    spatial_profile = data[reference_column,:] # Collapse along y-axis

    # Detect peaks (orders)
    # Apply Gaussian smoothing to reduce noise and enhance weak orders
    smoothed_profile = gaussian_filter(spatial_profile, sigma=5)
    
    # IDK why this is duplicate of above (but less aggressive) but it works
    smooth_data = data
    smoothed_image = gaussian_filter(smooth_data, sigma=2)
    

    # Determine a dynamic prominence threshold (adjust if necessary)
    peak_prominence = (np.max(smoothed_profile) - np.min(smoothed_profile)) * 0.01  # 2% of peak range
    min_height = np.max(smoothed_profile) * 0.007  # Ignore peaks below 0.7% of the highest peak

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
    # Get the flux from the orders
    fluxes = []
    for order in order_traces:
        column_positions, row_positions = order
        fluxs = []
        
        for i in range(-3,4):
            fluxs.append(data[row_positions, np.array(column_positions)+i])
        
        fluxs = np.array(fluxs)
        # TODO: fit gaussian shape across each order per column and extract flux weighted by the gaussian
        # for now, use some guess weights
        weights = np.array([0.05, 0.05, 0.15, 0.5, 0.15, 0.05, 0.05])
        # average the fluxes using the weights
        flux = np.mean(fluxs, axis=0)
        # flux = np.average(fluxs, axis=0, weights=weights)
        
        # reverse the order as the image is backwards
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

        # smooth the traced model aggresively
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

def find_tellurics_exp_type(orders):
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
    
    obj_type_guess = ' Stellar'
    if no_peaks == 0: #if no peaks are found, it must be a white
        obj_type_guess = ' White Light'
    
    if obj_type_guess == ' Stellar':
        if len(tellurics) <= 3: #if few orders are succesfully traced, it is likely a thorium as stellars tend to not struggle with order tracing in the red
            # TODO: Make this more robust (<= 3 is a bit arbitrary, should try to use size of image or something)
            obj_type_guess = ' Thorium Light'
    return tellurics, obj_type_guess