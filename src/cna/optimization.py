#!/usr/bin/env python3

import scipy.stats
import scipy.optimize
import scipy.signal
import random
import numpy as np
import logging
import plotly.graph_objects as go
import plotly

from src.utils.statistics import weighted_means

logger = logging.getLogger()

# parameters for modeling
# N_GAUSS = 4
HIST_RANGE = 0, 500

def smooth_triangle(data, degree):
    triangle=np.concatenate((np.arange(degree + 1), np.arange(degree)[::-1])) # up then down
    smoothed=[]

    for i in range(degree, len(data) - degree * 2):
        point=data[i:i + len(triangle)] * triangle
        smoothed.append(np.sum(point)/np.sum(triangle))
    if len(data) > degree:
        # Handle boundaries
        smoothed=[smoothed[0]]*int(degree + degree/2) + smoothed
        while len(smoothed) < len(data):
            smoothed.append(smoothed[-1])
        return smoothed
    else:
        return data

# def parse_segments_cov(filename):
#     cov = []
#     for line in open(filename, "r"):
#         cov.extend(map(float, line.strip().split(",")))
#     return cov


def peak_detection_optimization(args, input_segments, input_weights, tumor_cov):

    NUM_SAMPLES = 1000

    observed = input_segments
    weights = input_weights
    hist_mean = weighted_means(input_segments, input_weights)

    # computing observed histogram + normalization
    hist_bins = range(HIST_RANGE[0], HIST_RANGE[1] + 1)
    hist_scale = NUM_SAMPLES / 50
    observed_hist = np.histogram(observed, weights=weights, bins=hist_bins)[0]
    observed_hist = observed_hist / hist_scale

    if args.first_copy == 0:
        # Autocorrelation method
        corr = scipy.signal.correlate(observed_hist, observed_hist, mode="full")
        corr = corr[corr.size // 2:]  # autocorrelation, only positive shift, so getting right half of the array

        minima = scipy.signal.argrelmin(corr)
        #check that at least one local minimum smaller than mean coverage  exists -> self-correlation makes sense
        if len(minima[0]) > 0 and minima[0][0] < hist_mean:
            first_min = minima[0][0]
            corr_max = np.argmax(corr[first_min:]) + first_min

        #no good local minima -> just one coverage peak, likely balanced kariotype. set CN=1 to the peak value
        else:
            first_min = -1
            corr_max = np.argmax(observed_hist)

        # logger.info("First minimum %s", first_min)
        # logger.info("Max correlation peak %s", corr_max)

        # testing if 1/2 of the highest peak is also a peak. compare it with 1/4 and 3/4 (which should be valleys)
        half_peak = corr_max // 2  #change to 1 if dodn't want to consider neraset as half peak
        valley_left, valley_right = corr_max // 4, corr_max * 3 // 4
        is_half_peak = corr[half_peak] > max(corr[valley_left], corr[valley_right])
        # print(half_peak, valley_left, valley_right)
        # logger.info("Half peak: %s", is_half_peak)
    else:
        corr = observed_hist
        is_half_peak = False
        corr_max = args.first_copy

    if args.consider_wgd:
        is_half_peak = True

    single_copy_cov = corr_max + 0.012
    final_peaks_half = 0
    single_copy_cov_half = 0
    if is_half_peak:
        single_copy_cov_half = half_peak + 0.012

    #logger.info("Estimated single copy coverage: %s", single_copy_cov)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=list(range(len(corr))), y=corr, mode='lines', name='correlation peaks'))
    if is_half_peak:
        fig.add_shape(type="line", x0=half_peak, x1=half_peak, y0=0, y1=corr[half_peak], line=dict(color="yellow", width=2))
        fig.add_shape(type="line", x0=corr_max, x1=corr_max, y0=0, y1=corr[corr_max], line=dict(color="red", width=2))
    else:
        fig.add_shape(type="line", x0=corr_max, x1=corr_max, y0=0, y1=corr[corr_max], line=dict(color="red", width=2))
    plotly.offline.plot(fig, filename=args.out_dir_plots + '/' + args.genome_name +'_optimized_peak.html', auto_open=False)

    last_copy_state = int(max(observed) // single_copy_cov + 1)
    final_peaks = [i * single_copy_cov for i in range(0, last_copy_state + 1)]

    if is_half_peak:
        last_copy_state_half = int(max(observed) // single_copy_cov_half + 1)
        final_peaks_half = [i * single_copy_cov_half for i in range(0, last_copy_state_half + 1)]

    final_peaks_subclonal = [single_copy_cov//2] + [i * single_copy_cov+(single_copy_cov//2) for i in range(1, last_copy_state)]

    return first_min, final_peaks, is_half_peak, final_peaks_half, final_peaks_subclonal, np.arange(0, 500), observed_hist, single_copy_cov, single_copy_cov_half
