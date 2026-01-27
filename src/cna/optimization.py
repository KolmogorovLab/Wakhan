#!/usr/bin/env python3

import scipy.stats
import scipy.optimize
import scipy.signal
import scipy.spatial
import random
import numpy as np
import logging
import plotly.graph_objects as go
import plotly

import ruptures as rpt
from collections import defaultdict
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


"""
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
"""


"""
def gaussian_mixture(x, peak_positions, sigma=1.0, amplitudes=None):
    if amplitudes is None:
        amplitudes = np.ones(len(peak_positions))

    if np.isscalar(sigma):
        sigma = np.full(len(peak_positions), sigma)

    mixture = np.zeros_like(x, dtype=float)

    for pos, amp, sig in zip(peak_positions, amplitudes, sigma):
        mixture += amp * scipy.stats.norm.pdf(x, loc=pos, scale=sig)

    return mixture
"""


class ProfileException(Exception):
    pass


def cn_one_inference(input_segments, input_weights, cov_ploidy):
    observed = input_segments
    weights = input_weights

    # computing observed coverage histogram
    hist_max = np.quantile(input_segments, 0.50) * 2
    hist_range = (0, hist_max)
    hist_mean = weighted_means(input_segments, input_weights)
    HIST_RATE = hist_max / 100
    print("Hist range", hist_range, "rate", HIST_RATE)

    hist_bins = np.arange(hist_range[0], hist_range[1] + 1, HIST_RATE)
    observed_hist = np.histogram(observed, weights=weights, bins=hist_bins)[0]

    #gentle smoothing
    raw_hist = np.copy(observed_hist)
    observed_hist = scipy.signal.savgol_filter(observed_hist, window_length=5, polyorder=3)

    #computing autocorrelation to find modes
    corr = scipy.signal.correlate(observed_hist, observed_hist, mode="full")
    corr = corr[corr.size // 2:]  # autocorrelation, only positive shift, so getting right half of the array

    #smoothing correlation a bit
    smooth_corr = scipy.signal.savgol_filter(corr, window_length=5, polyorder=3)

    #peaks_mixture = None
    first_min = 0
    max_corr_peak = np.argmax(smooth_corr)
    scaled_peaks = []

    try:
        #check that at least one local minimum smaller than mean coverage  exists -> self-correlation makes sense
        #if yes, will find peaks after th minimum. if not, can't infer much from this profile
        minima = scipy.signal.argrelmin(smooth_corr)
        peaks, stats = scipy.signal.find_peaks(smooth_corr, distance=5,
                                               width=3, rel_height=1.0)

        if len(minima[0]) == 0 or minima[0][0] >= hist_mean:
            cn_one = np.argmax(observed_hist) * HIST_RATE / cov_ploidy
            raise ProfileException("No local minima, likely a single peak. Setting CN=1 to histogram maximum")

        #making sure all peaks are past minima, although it always should be the case
        first_min = minima[0][0]
        peaks = [p for p in peaks if p > first_min]

        if len(peaks) == 0:
            cn_one = np.argmax(observed_hist) * HIST_RATE / cov_ploidy
            raise ProfileException("No histogram peaks found, defaulting to histogram maximum")

        #Have a first local mimima and at least one peak in the correlation spectra
        max_corr_peak = max(peaks, key=lambda p: smooth_corr[p])

        #filtering very small peaks
        HEIGHT_RATE = 0.3
        peaks = [p for p in peaks if smooth_corr[p] >= smooth_corr[max_corr_peak] * HEIGHT_RATE]
        scaled_peaks = [round(p * HIST_RATE, 2) for p in peaks]

        print("First minimum", first_min * HIST_RATE)
        print("Correlation peaks", scaled_peaks)
        print("Max correlation peak", max_corr_peak * HIST_RATE)

        cn_one = scaled_peaks[0]

        """
        if len(peaks) == 1:
            cn_one = max_corr_peak * HIST_RATE
            raise ProfileException("Just one correlation peak, using it for CN=1")

        #now, we have 2+ peaks, need to find the distance between them
        peaks_mixture = gaussian_mixture(hist_bins[:-1], scaled_peaks, sigma=(1 * HIST_RATE),
                                         amplitudes=[smooth_corr[int(p / HIST_RATE)] for p in scaled_peaks])

        gauss_corr = scipy.signal.correlate(peaks_mixture, peaks_mixture, mode="full")[len(peaks_mixture) - 1:]
        peaks_gauss, stats_gauss = scipy.signal.find_peaks(gauss_corr)
        if len(peaks_gauss) == 0:
            cn_one = np.argmax(observed_hist) * HIST_RATE / cov_ploidy
            peaks_gauss = [max(gauss_corr)]
            raise ProfileException("No histogram peaks found, defaulting to histogram maximum")
        cn_one = peaks_gauss[0] * HIST_RATE
        """

    except ProfileException as e:
        print(e)
        pass

    print("Haploid CN=1 estimate", cn_one)

    #some visualization
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, layout='constrained')
    ax1.bar(hist_bins[:-1], raw_hist, width=HIST_RATE, edgecolor='black')
    ax1.set_title("Segment coverage histogram")

    ax2.bar(hist_bins[:-1], observed_hist, width=HIST_RATE, edgecolor='black')
    ax2.set_title("Smoothed segment coverage histogram")

    corr_top = max(smooth_corr[max_corr_peak], smooth_corr[first_min]) * 2
    ax3.plot(hist_bins[:-1], corr)
    ax3.set_ylim(0, corr_top)
    ax3.set_title("Histogram self-correlation")

    ax4.plot(hist_bins[:-1], smooth_corr)
    ax4.set_ylim(0, corr_top)
    ax4.set_title("Smoothed self-correlation")
    ax4.plot([first_min * HIST_RATE, first_min * HIST_RATE],
             [0, smooth_corr[first_min]], color="red", linewidth=2)
    for p in scaled_peaks:
        top = smooth_corr[int(p / HIST_RATE)]
        ax4.plot([p, p], [0, top], color="green", linewidth=2)

    plt.show()

    return cn_one, first_min * HIST_RATE, hist_bins, observed_hist


"""
def _custom_corr(array):
    corr = []
    for i in range(1, len(array) - 1):
        #shifted = np.roll(array, i)
        #shifted[:i] = 0
        part_a = array[:-i]
        part_b = array[i:]

        #corr.append(scipy.spatial.distance.braycurtis(part_a, part_b))
        pearson = scipy.stats.pearsonr(part_a, part_b)
        print(pearson)
        if not np.isnan(pearson.statistic):
            corr.append(1 + pearson.statistic)
        else:
            corr.append(0)
    corr = [0] + corr + [0]
    print(corr)
    return np.array(corr)
"""


def peak_detection_optimization(args, input_segments, input_weights, tumor_cov):
    COV_PLOIDY = 2  #change to 1 if using haploid coverage as input
    cn_one, first_min, hist_x, hist_y = cn_one_inference(input_segments, input_weights, COV_PLOIDY)

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

    return (first_min, final_peaks, is_half_peak, final_peaks_half, final_peaks_subclonal,
            hist_x, hist_y, single_copy_cov, single_copy_cov_half)


class Dummy:
    pass


def parse_coverage_bed(filename, avg_window, phased):
    """
    Reading Wakhan bed with phase-corrected binned coverage
    """
    bin_size = None

    segments = []
    weights = []

    raw_segments = []
    with open(filename, "r") as fin:
        hp1 = []
        hp2 = []
        hp3 = []
        for line in fin:
            fields = line.split("\t")
            if bin_size is None:
                bin_size = int(fields[2]) - int(fields[1])
                print("Bin size:", bin_size)

            min_cov, max_cov = float(fields[3]), float(fields[4])   #phased
            if min_cov + max_cov > 0:
                hp1.append(min_cov)
                hp2.append(max_cov)
                hp3.append(min_cov + max_cov)
            """
            if float(fields[3]) > 0:
                hp1.append(0)
                hp2.append(0)
                hp3.append(float(fields[3]))
            """

            if len(hp1) >= round(avg_window / bin_size):
                if phased:
                    raw_segments.append(hp1)
                    raw_segments.append(hp2)
                else:
                    raw_segments.append(hp3)
                hp1 = []
                hp2 = []
                hp3 = []

    std_devs = []
    for seg in raw_segments:
        std_devs.append(np.std(seg))
    max_dev = np.quantile(std_devs, 0.90)
    print("Std median:", np.median(std_devs))
    print("Std Q90", max_dev)

    for seg in raw_segments:
        if np.std(seg) < max_dev:
            segments.append(np.median(seg))

    #print(len(raw_segments))
    #for seg in raw_segments:
    #    print(np.median(seg), np.std(seg))

    weights = [1] * len(segments)
    return segments, weights


def parse_coverage_bed_cpd(filename, phased):
    bin_size = None
    cov_by_chrom = defaultdict(list)

    with open(filename, "r") as fin:
        for line in fin:
            fields = line.split("\t")

            if bin_size is None:
                bin_size = int(fields[2]) - int(fields[1])
                print("Bin size:", bin_size)

            min_cov, max_cov = float(fields[3]), float(fields[4])   #phased
            if min_cov + max_cov > 0:
                cov_by_chrom[fields[0]].append(min_cov + max_cov)

    median_cov = np.median(sum(cov_by_chrom.values(), []))
    #print("Median coverage:", median_cov)

    MIN_SEG_LEN_BP = 5000000
    MIN_JUMP_BP = 1000000
    MIN_SEG_LEN = MIN_SEG_LEN_BP // bin_size

    chroms = ["chr" + str(i) for i in range(1, 23)]
    fig, subplots = plt.subplots(len(chroms) // 2, 2, sharey=False)

    all_segments = []
    for i, ch in enumerate(chroms):
        sp = subplots[i // 2][i % 2]
        coverage = np.array(cov_by_chrom[ch])
        y_max = np.quantile(coverage, 0.90)

        pen_auto = np.log(len(coverage)) * median_cov ** 2 / 5
        print(ch + " ", end="", flush=True)

        algo_pelt = rpt.Pelt(model="l2", min_size=MIN_SEG_LEN,
                             jump=(MIN_JUMP_BP // bin_size))
        algo_pelt.fit(coverage)
        breakpoints = algo_pelt.predict(pen=pen_auto)

        breakpoints = [0] + breakpoints + [len(coverage)]
        for start, end in zip(breakpoints[:-1], breakpoints[1:]):
            segment = coverage[start : end]
            if len(segment) > MIN_SEG_LEN:
                all_segments.append(segment)

        ### visualization
        sp.plot(coverage)
        for bp in breakpoints:
            sp.plot([bp, bp], [0, y_max], color="red")
        sp.set_ylim(0, y_max)
        sp.set_ylabel(ch)
        ###

    print("")
    print("Total cpd fragments:", len(all_segments))
    plt.show()

    """
    std_devs = []
    for seg in all_segments:
        std_devs.append(np.std(seg))
    max_dev = np.quantile(std_devs, 0.90)
    print("Std median:", np.median(std_devs))
    print("Std Q90", max_dev)
    """
    max_dev = 9999  #filtering disabled

    filtered_segments = []
    filtered_weights = []
    for seg in all_segments:
        if np.std(seg) < max_dev:
            filtered_segments.append(np.median(seg))
            filtered_weights.append(np.log(len(seg)))
            #filtered_weights.append(len(seg))

    return filtered_segments, filtered_weights


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
    WINDOW = 5000000    #5Mb
    segments, weights = parse_coverage_bed_cpd(sys.argv[1], None)
    #segments, weights = parse_coverage_bed(sys.argv[1], WINDOW, phased=False)

    args = Dummy()
    args.first_copy = 0
    args.consider_wgd = False
    args.out_dir_plots = "."
    args.genome_name = "Sample"
    peak_detection_optimization(args, segments, weights, None)
