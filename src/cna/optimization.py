#!/usr/bin/env python3

import scipy.stats
import scipy.optimize
import scipy.signal
import random
import numpy as np
import logging
import plotly.graph_objects as go
import plotly

import matplotlib.pyplot as plt


def weighted_means(vals, weights):
    if not len(vals) == len(weights) or sum(weights) == 0:
        return 0
    weighted_vals = []
    vals_n_weights = [(vals[i], weights[i]) for i in range(0, len(weights))]
    for tup in vals_n_weights:
        weighted_vals.append(round(tup[0] * tup[1], 3))
    return sum(weighted_vals)/sum(weights)
import sys


logger = logging.getLogger()

# parameters for modeling
# N_GAUSS = 4

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


def gaussian_mixture(x, peak_positions, sigma=1.0, amplitudes=None):
    if amplitudes is None:
        amplitudes = np.ones(len(peak_positions))

    if np.isscalar(sigma):
        sigma = np.full(len(peak_positions), sigma)

    mixture = np.zeros_like(x, dtype=float)

    for pos, amp, sig in zip(peak_positions, amplitudes, sigma):
        mixture += amp * scipy.stats.norm.pdf(x, loc=pos, scale=sig)

    return mixture


class ProfileException(Exception):
    pass


def cn_one_inference(input_segments, input_weights, cov_ploidy):
    observed = input_segments
    weights = input_weights
    hist_mean = weighted_means(input_segments, input_weights)

    hist_max = np.quantile(input_segments, 0.50) * 2
    hist_range = (1, hist_max)
    HIST_RATE = hist_max / 100
    print("Hist range", hist_range, "rate", HIST_RATE)

    # computing observed histogram + normalization
    hist_bins = np.arange(hist_range[0], hist_range[1] + 1, HIST_RATE)
    observed_hist = np.histogram(observed, weights=weights, bins=hist_bins)[0]

    # Autocorrelation method
    corr = scipy.signal.correlate(observed_hist, observed_hist, mode="full")
    corr = corr[corr.size // 2:]  # autocorrelation, only positive shift, so getting right half of the array

    #smoothing correlation a bit
    smooth_corr = scipy.signal.savgol_filter(corr, window_length=5, polyorder=3)

    peaks_mixture = None
    first_min = 0
    max_corr_peak = np.argmax(smooth_corr)
    scaled_peaks = []

    try:
        #check that at least one local minimum smaller than mean coverage  exists -> self-correlation makes sense
        #if yes, will find peaks after th minimum. if not, can't infer much from this profile
        minima = scipy.signal.argrelmin(smooth_corr)
        peaks, stats = scipy.signal.find_peaks(smooth_corr, distance=(5 / HIST_RATE),
                                               width=1 / HIST_RATE, rel_height=1.0)

        if len(minima[0]) == 0 or minima[0][0] >= hist_mean:
            cn_one = np.argmax(observed_hist) * HIST_RATE / cov_ploidy
            raise ProfileException("No local minima, likely a single peak. Setting CN=1 to histogram maximum")

        #making sure all peaks are past minima, although it always should be the case
        first_min = minima[0][0]
        peaks = [p for p in peaks if p > first_min]
        scaled_peaks = [round(hist_bins[p], 1) for p in peaks]
        print("Correlation peaks", scaled_peaks)

        if len(peaks) == 0:
            cn_one = np.argmax(observed_hist) * HIST_RATE / cov_ploidy
            raise ProfileException("No histogram peaks found, defaulting to histogram maximum")

        #Have a first local mimima and at least one peak in the correlation spectra
        max_corr_peak = max(peaks, key=lambda p: smooth_corr[p])
        print("First minimum", first_min * HIST_RATE)
        print("Max correlation peak", max_corr_peak * HIST_RATE)

        #filtering very small peaks
        HEIGHT_RATE = 0.3
        peaks = [p for p in peaks if smooth_corr[p] >= smooth_corr[max_corr_peak] * HEIGHT_RATE]

        if len(peaks) == 1:
            cn_one = max_corr_peak * HIST_RATE
            raise ProfileException("Just one correlation peak, using it for CN=1")

        #now, we have 2+ peaks, need to find the distance between them
        peaks_mixture = gaussian_mixture(hist_bins[:-1], scaled_peaks, sigma=1 * HIST_RATE,
                                         amplitudes=[smooth_corr[int(p / HIST_RATE)] for p in scaled_peaks])

        gauss_corr = scipy.signal.correlate(peaks_mixture, peaks_mixture, mode="full")[len(peaks_mixture) - 1:]
        peaks_gauss, stats_gauss = scipy.signal.find_peaks(gauss_corr)
        cn_one = hist_bins[peaks_gauss[0]]

    except ProfileException as e:
        print(e)
        pass

    print("Haploid CN=1 estimate", cn_one)

    #some visualization
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1)
    ax1.bar(hist_bins[:-1], observed_hist, width=HIST_RATE, edgecolor='black')

    #ax2.plot(hist_bins[:-1], corr)
    #ax2.set_ylim(0, corr[corr_max] * 2)
    ax2.plot(hist_bins[:-1], smooth_corr)
    ax2.set_ylim(0, smooth_corr[max_corr_peak] * 2)
    ax2.plot([first_min * HIST_RATE, first_min * HIST_RATE],
             [0, smooth_corr[first_min]], color="red", linewidth=2)
    for p in scaled_peaks:
        top = smooth_corr[int(p / HIST_RATE)]
        ax2.plot([p, p], [0, top], color="green", linewidth=2)

    if peaks_mixture is not None:
        ax3.plot(hist_bins[:-1], peaks_mixture)

        ax4.plot(hist_bins[:-1], gauss_corr)
        ax4.plot([cn_one, cn_one], [0, gauss_corr[peaks_gauss[0]]], color="red", linewidth=2)

    plt.show()

    return cn_one, first_min * HIST_RATE


def peak_detection_optimization(args, input_segments, input_weights, tumor_cov):
    COV_PLOIDY = 2  #change to 1 if using haploid coverage as input
    cn_one, first_min, max_coverage = cn_one_inference(input_segments, input_weights, COV_PLOIDY)

    #overwriting with command line arguments if needed
    is_half_peak = False
    if args.consider_wgd:
        print("CN=1 overwritten by command line argument")
        is_half_peak = True

    if args.first_copy != 0:
        print("CN=1 overwritten by command line argument")
        is_half_peak = False
        cn_one = args.first_copy
    ###

    single_copy_cov = cn_one + 0.012
    final_peaks_half = 0
    single_copy_cov_half = 0

    if is_half_peak:
        single_copy_cov_half = half_peak + 0.012

    last_copy_state = int(max(input_segments) // single_copy_cov + 1)
    final_peaks = [i * single_copy_cov for i in range(0, last_copy_state + 1)]

    if is_half_peak:
        last_copy_state_half = int(max(input_segments) // single_copy_cov_half + 1)
        final_peaks_half = [i * single_copy_cov_half for i in range(0, last_copy_state_half + 1)]

    final_peaks_subclonal = [single_copy_cov//2] + [i * single_copy_cov+(single_copy_cov//2) for i in range(1, last_copy_state)]

    return (scaled_first_min, final_peaks, is_half_peak, final_peaks_half, final_peaks_subclonal,
            np.arange(0, 500), observed_hist, single_copy_cov, single_copy_cov_half)


class Dummy:
    pass

if __name__ == "__main__":
    segments = []
    weights = []

    """
    #parsing form optimization_segments.csv
    with open(sys.argv[1], "r") as fin:
        for line in fin:
            fields = line.strip().split()
            segments.append(float(fields[0]))
            weights.append(int(fields[1]))
            #weights.append(1)
    """

    #parsing from coverage.ps
    WINDOW = 100    #5mb with 50kb bins
    with open(sys.argv[1], "r") as fin:
        hp1 = []
        hp2 = []
        hp3 = []
        for line in fin:
            fields = line.split("\t")
            #min_cov, max_cov = sorted([float(fields[3]), float(fields[4])])
            min_cov, max_cov = float(fields[3]), float(fields[4])
            hp1.append(min_cov)
            hp2.append(max_cov)
            hp3.append(min_cov + max_cov)
            #hp1.append(float(fields[3]) + float(fields[4]))
            if len(hp1) >= WINDOW:
                #segments.append(np.median(hp1))
                #segments.append(np.median(hp2))
                segments.append(np.median(hp3))
                hp1 = []
                hp2 = []
                hp3 = []
    weights = [1] * len(segments)

    args = Dummy()
    args.first_copy = 0
    args.consider_wgd = False
    args.out_dir_plots = "."
    args.genome_name = "Sample"
    peak_detection_optimization(args, segments, weights, None)
