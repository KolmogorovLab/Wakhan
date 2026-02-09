#!/usr/bin/env python3

import scipy.signal
import random
import numpy as np
import logging
import os

import ruptures as rpt
from collections import defaultdict
import matplotlib.pyplot as plt


logger = logging.getLogger()
logging.getLogger('matplotlib').setLevel(logging.WARNING)


class ProfileException(Exception):
    pass


def _weighted_means(vals, weights):
    if not len(vals) == len(weights) or sum(weights) == 0:
        return 0
    weighted_vals = []
    vals_n_weights = [(vals[i], weights[i]) for i in range(0, len(weights))]
    for tup in vals_n_weights:
        weighted_vals.append(round(tup[0] * tup[1], 3))
    return sum(weighted_vals)/sum(weights)


def cn_one_inference(input_segments, input_weights, phased, plot_path):
    """
    Input: segments set with given coverage and length. Haploid coverage should use cov_ploidy=1, diploid=2
    """
    observed = input_segments
    weights = input_weights
    cov_ploidy = 1 if phased else 2

    # computing observed coverage histogram
    hist_max = np.quantile(input_segments, 0.50) * 2
    hist_range = (0, hist_max)
    hist_mean = _weighted_means(input_segments, input_weights)
    HIST_RATE = hist_max / 100

    logger.info("CN=1 inference")
    logger.info("Hist range {}, rate: {}".format(hist_range, HIST_RATE))

    hist_bins = np.arange(hist_range[0], hist_range[1] + 1, HIST_RATE)
    observed_hist = np.histogram(observed, weights=weights, bins=hist_bins)[0]
    observed_hist[0] = 0

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
    max_peak_id = None
    scaled_peaks = []
    cn_one_alt = None

    try:
        #check that at least one local minimum smaller than mean coverage  exists -> self-correlation makes sense
        #if yes, will find peaks after th minimum. if not, can't infer much from this profile
        minima = scipy.signal.argrelmin(smooth_corr)
        peaks, peak_prop = scipy.signal.find_peaks(smooth_corr, distance=5,
                                                   width=1, rel_height=1.0)

        if len(minima[0]) == 0 or minima[0][0] * HIST_RATE >= hist_mean:
            cn_one = np.argmax(observed_hist) * HIST_RATE / cov_ploidy
            raise ProfileException("No local minima, likely a single peak. Setting CN=1 to histogram maximum")

        MIN_COV_TO_CORR = 0.25
        min_corr = np.argmax(observed_hist) * MIN_COV_TO_CORR
        logger.info("Minimum correlation %4.2f", min_corr * HIST_RATE)

        #making sure all peaks are past minima and minimum correlation (that depends on coverage)
        first_min = minima[0][0]
        logger.info("First minimum {}".format(first_min * HIST_RATE))
        filtered_peak_ids = [i for i in range(len(peaks)) if peaks[i] > max(first_min, min_corr)]

        #Have a first local mimima and at least one peak in the correlation spectra
        logger.info("Initial peaks: {}".format([round(peaks[i] * HIST_RATE, 2) for i in filtered_peak_ids]))

        #Only keep peaks with relative prominence >0.25
        MIN_REL_PROM = 0.25
        filtered_peak_ids = [i for i in filtered_peak_ids if peak_prop["prominences"][i] / smooth_corr[peaks[i]] > MIN_REL_PROM]
        logger.info("Peaks w/ prom > 0.25: {}".format([round(peaks[i] * HIST_RATE, 2) for i in filtered_peak_ids]))

        if len(filtered_peak_ids) == 0:
            cn_one = np.argmax(observed_hist) * HIST_RATE / cov_ploidy
            raise ProfileException("No histogram peaks found, defaulting to histogram maximum")

        max_peak_id = max(filtered_peak_ids, key=lambda p: smooth_corr[peaks[p]])
        logger.info("Max peak: {}".format(peaks[max_peak_id] * HIST_RATE))

        #Now, filter peaks that are too small, compared to the top peak
        MIN_TO_PEAK = 0.5
        MIN_TO_MAX = 0.2
        min_height = max(smooth_corr[peaks[max_peak_id]] * MIN_TO_PEAK,
                         max(smooth_corr[first_min:]) * MIN_TO_MAX)
        #min_height = max(smooth_corr[peaks[max_peak_id]], smooth_corr[first_min]) * MIN_HEIGHT_TO_TOP
        #min_height = smooth_corr[peaks[max_peak_id]] * MIN_HEIGHT_TO_TOP
        #min_height = max(smooth_corr[first_min:]) * MIN_HEIGHT_TO_TOP


        filtered_peak_ids = [i for i in filtered_peak_ids if smooth_corr[peaks[i]] > min_height]
        logger.info("Peaks >0.5 of top/saddle: {}".format([round(peaks[i] * HIST_RATE, 2) for i in filtered_peak_ids]))

        if len(filtered_peak_ids) == 0:
            cn_one = np.argmax(observed_hist) * HIST_RATE / cov_ploidy
            max_peak_id = None
            raise ProfileException("No histogram peaks found, defaulting to histogram maximum")

        scaled_peaks = [round(peaks[p] * HIST_RATE, 2) for p in filtered_peak_ids]

        if max_peak_id == filtered_peak_ids[0]:
            cn_one = peaks[filtered_peak_ids[0]] * HIST_RATE
        else:
            cn_one = peaks[max_peak_id] * HIST_RATE
            cn_one_alt = peaks[filtered_peak_ids[0]] * HIST_RATE

    except ProfileException as e:
        logger.info(e)
        pass

    logger.info("Best CN=1 estimate: {}".format(cn_one))
    logger.info("Alt CN=1 estimate: {}".format(cn_one_alt))

    #some visualization
    if plot_path is not None:
        fig, (ax1, ax2) = plt.subplots(2, 1, layout='constrained')
        ax1.bar(hist_bins[:-1], observed_hist, width=HIST_RATE, edgecolor='black')
        ax1.set_title("Segment coverage histogram")

        #if max_peak_id is not None:
        #    corr_top = max(smooth_corr[peaks[max_peak_id]], smooth_corr[first_min]) * 2
        #else:
        #    corr_top = smooth_corr[first_min] * 2
        corr_top = max(smooth_corr[first_min:]) * 2

        ax2.plot(hist_bins[:-1], smooth_corr)
        ax2.set_ylim(0, corr_top)
        ax2.set_title("Self-correlation")
        ax2.plot([first_min * HIST_RATE, first_min * HIST_RATE],
                 [0, corr_top], color="red", linewidth=2)
        ax2.plot([min_corr * HIST_RATE, min_corr * HIST_RATE],
                 [0, corr_top], color="orange", linewidth=2)
        for p in scaled_peaks:
            top = smooth_corr[int(p / HIST_RATE)]
            ax2.plot([p, p], [0, top], color="green", linewidth=2)

        if plot_path == "show":
            plt.show()
        else:
            plt.savefig(plot_path, dpi=300)
    ###

    return cn_one, cn_one_alt, first_min * HIST_RATE, hist_bins, observed_hist


def peak_detection_optimization(args, input_segments, input_weights, phased):
    """
    Older functions that is called from main
    """
    cn_peaks_plot = os.path.join(args.out_dir_plots, 'coverage_data', 'cn_peaks.png')
    cn_one, cn_one_alt, first_min, hist_x, hist_y = cn_one_inference(input_segments, input_weights, phased, plot_path=cn_peaks_plot)
    is_half_peak = (cn_one_alt is not None)

    #overwriting with command line arguments if needed
    #if args.consider_wgd:
    #    logger.info("CN=1 overwritten by command line argument")
    #    is_half_peak = True

    if args.first_copy != 0:
        logger.info("CN=1 overwritten by command line argument")
        is_half_peak = False
        cn_one = args.first_copy
    ###

    single_copy_cov = cn_one
    final_peaks_half = 0
    single_copy_cov_half = 0
    if is_half_peak:
        single_copy_cov_half = cn_one_alt

    last_copy_state = int(max(input_segments) // single_copy_cov + 1)
    final_peaks = [i * single_copy_cov for i in range(0, last_copy_state + 1)]

    if is_half_peak:
        last_copy_state_half = int(max(input_segments) // single_copy_cov_half + 1)
        final_peaks_half = [i * single_copy_cov_half for i in range(0, last_copy_state_half + 1)]

    final_peaks_subclonal = [single_copy_cov//2] + [i * single_copy_cov+(single_copy_cov//2) for i in range(1, last_copy_state)]

    return (final_peaks, is_half_peak, final_peaks_half, final_peaks_subclonal,
            hist_x, hist_y, single_copy_cov, single_copy_cov_half)


def parse_coverage_bed_windows(filename, phased):
    """
    Reading Wakhan bed with phase-corrected binned coverage, average by larger window
    """
    WINDOW = 5000000    #5Mb
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
                logger.info("CPD bin size: {}".format(bin_size))

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

            if len(hp1) >= round(WINDOW / bin_size):
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
    logger.info("Std median: {}".format(np.median(std_devs)))
    logger.info("Std Q90: {}".format(max_dev))

    for seg in raw_segments:
        if np.std(seg) < max_dev:
            segments.append(np.median(seg))

    weights = [1] * len(segments)
    return segments, weights


def parse_coverage_bed_cpd(filename, phased, plot_path):
    """
    Parses Wakhan phased_coverage file, fragments with CPD and returns segments coverage and weights
    """

    MIN_SEG_LEN_BP = 5000000
    MIN_JUMP_BP = 1000000

    bin_size = None
    cov_by_chrom = defaultdict(list)
    cov_by_chrom_hp1 = defaultdict(list)
    cov_by_chrom_hp2 = defaultdict(list)

    with open(filename, "r") as fin:
        for line in fin:
            fields = line.split("\t")

            if bin_size is None:
                bin_size = int(fields[2]) - int(fields[1])
                logger.info("CPD bin size: {}".format(bin_size))

            min_cov, max_cov = float(fields[3]), float(fields[4])   #phased
            #min_cov, max_cov = float(fields[3]), 0
            if min_cov + max_cov > 0:
                cov_by_chrom_hp1[fields[0]].append(min_cov)
                cov_by_chrom_hp2[fields[0]].append(max_cov)

    for chr in cov_by_chrom_hp1:
        if not phased:
            cov_by_chrom[chr] = [x + y for (x, y) in zip(cov_by_chrom_hp1[chr], cov_by_chrom_hp2[chr])]
        else:
            cov_by_chrom[chr] = cov_by_chrom_hp1[chr] + cov_by_chrom_hp2[chr]

    median_cov = np.median(sum(cov_by_chrom.values(), []))
    MIN_SEG_LEN = MIN_SEG_LEN_BP // bin_size
    chroms = ["chr" + str(i) for i in range(1, 23)]

    if plot_path is not None:
        fig, subplots = plt.subplots(len(chroms) // 2, 2, sharey=False)

    logger.info("Segmenting input coverage")
    all_segments = []
    for i, ch in enumerate(chroms):
        coverage = np.array(cov_by_chrom[ch])
        y_max = np.quantile(coverage, 0.90)

        pen_auto = np.log(len(coverage)) * median_cov ** 2 / 1
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
        if plot_path is not None:
            sp = subplots[i // 2][i % 2]
            sp.plot(coverage)
            for bp in breakpoints:
                sp.plot([bp, bp], [0, y_max], color="red")
            sp.set_ylim(0, y_max)
            sp.set_ylabel(ch)

            for bp1, bp2 in zip(breakpoints[:-1], breakpoints[1:]):
                segment = coverage[bp1 : bp2]
                if len(segment) > 1:
                    seg_cov = np.median(segment)
                    sp.plot([bp1, bp2], [seg_cov, seg_cov], color="red")
        ###

    logger.info("Total cpd fragments: {}".format(len(all_segments)))

    if plot_path is not None:
        fig.set_size_inches(12, 8)
        if plot_path == "show":
            plt.show()
        else:
            plt.savefig(plot_path, dpi=300)

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
            #filtered_weights.append(np.log(len(seg)))
            filtered_weights.append(len(seg))
            #filtered_weights.append(len(seg) ** 2)

    return filtered_segments, filtered_weights


#entry point for debugging
if __name__ == "__main__":
    import sys
    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: " "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: " "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    console_log.setLevel(logging.INFO)
    logger.setLevel(logging.INFO)
    logger.addHandler(console_log)

    """
    segments = []
    weights = []
    #parsing form optimization_segments.csv
    with open(sys.argv[1], "r") as fin:
        for line in fin:
            fields = line.strip().split()
            segments.append(float(fields[0]))
            weights.append(int(fields[1]))
            #weights.append(1)
    """

    PHASED = True
    #segments, weights = parse_coverage_bed_windows(sys.argv[1], phased=PHASED)
    segments, weights = parse_coverage_bed_cpd(sys.argv[1], phased=PHASED, plot_path="show")
    cn_one_inference(segments, weights, phased=PHASED, plot_path="show")
