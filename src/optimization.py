#!/usr/bin/env python3

import sys
import scipy.stats
import scipy.optimize
import scipy.signal
#from matplotlib import pyplot as plt
import random
import numpy as np

# parameters for modeling
# N_GAUSS = 4
HIST_RANGE = 0, 100

def merge_similar_elements(data, threshold):
  """
  This function takes a list of data and a threshold.
  It groups elements whose difference is less than the threshold and returns a list
  containing the mean of each group.

  Args:
      data: A list of numerical data.
      threshold: The maximum difference allowed for elements to be considered similar.

  Returns:
      A list containing the mean of each group of similar elements.
  """

  if not data:
    return []

  merged_data = []
  current_group = []
  current_value = data[0]

  for value in data[1:]:
    if abs(value - current_value) <= threshold:
      current_group.append(value)
    else:
      # Add the mean of the previous group
      merged_data.append(sum(current_group) / len(current_group))
      current_group = [value]
      current_value = value

  # Add the mean of the last group
  merged_data.append(sum(current_group) / len(current_group))

  return merged_data

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

def simulate(sample_size):
    """
    Simulates sum of X gaussians
    """
    MEANS = [10, 25, 25, 55]
    SIG = 3

    distr = [scipy.stats.norm(loc=m, scale=SIG) for m in MEANS]
    sampled = []
    for _ in range(sample_size):
        distr_id = random.randint(0, len(MEANS) - 1)
        sampled.append(distr[distr_id].rvs())

    return sampled


def model_distribution(offset, peak_dist, sig, abundance):
    """
    Generates theoretical distribution function, sum of multiple gaussians.
    The return type is still a function with one input agrument (x)
    """
    # abundance = [abundance[0]] * len(abundance)
    # peak_dist = 15
    norms = [scipy.stats.norm(loc=(offset + i * peak_dist), scale=sig) for i in range(len(abundance))]
    return lambda x: sum([abundance[i] * norms[i].pdf(x) for i in range(len(abundance))])


def fit_optimization(observed_hist, num_gauss, guess_offset, guess_distance, guess_sig):
    """
    Fits the observed historgram with given initial guesses
    """

    def loss_function(offset, peak_dist, sig, abundance):
        theor_fixed = model_distribution(offset, peak_dist, sig, abundance)
        return 10 * sum([(observed_hist[x] - theor_fixed(x)) ** 2
                         for x in range(*HIST_RANGE)])

    def optimization_wrapper(x):
        return loss_function(x[0], x[1], x[2], x[3:])

    num_samples = sum(observed_hist)
    x0 = [guess_offset, guess_distance, guess_sig] + [num_samples / num_gauss] * num_gauss
    x0 = np.array(x0)

    # positive bounds for all variablex, x > 0
    bounds = [(0, None) for x in x0]

    # constaint: offset < peak_dist
    # cons_matrix = [-1, 1, 0] + [0] * num_gauss
    # constraint = scipy.optimize.LinearConstraint(np.array(cons_matrix), lb=0, ub="inf")

    res = scipy.optimize.minimize(optimization_wrapper, x0, method="Powell", bounds=bounds)
    # res = scipy.optimize.minimize(optimization_wrapper, x0, method="nelder-mead", bounds=bounds)

    # print(res)
    return res


def parse_segments_cov(filename):
    cov = []
    for line in open(filename, "r"):
        cov.extend(map(float, line.strip().split(",")))
    return cov


def parse_segments_cov_bins(filename):
    cov = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            elements = line.strip().split("\t")
            hp1 = float(elements[3])
            hp2 = float(elements[4])
            cov.append(hp1)
            cov.append(hp2)
    return cov


def peak_detection_optimization(input_segments, input_weights):

    NUM_SAMPLES = 1000
    # observed = simulate(NUM_SAMPLES)
    # observed = parse_segments_cov(sys.argv[1])
    observed = input_segments#parse_segments_cov_bins(sys.argv[1])
    weights = input_weights
    # observed = [x1 for x1 in observed if x1 > 250 and x1 < 300]
    #NUM_SAMPLES = len(observed)
    print(observed)

    # computing observed histogram + normalization
    hist_bins = range(HIST_RANGE[0], HIST_RANGE[1] + 1)
    hist_scale = NUM_SAMPLES / 50
    observed_hist = np.histogram(observed, weights=weights, bins=hist_bins)[0]
    observed_hist = observed_hist / hist_scale

    xx = range(*HIST_RANGE)
    #plt.bar(xx, observed_hist, color="blue")

    # peak peaking (independent from optimization approach)
    peaks = scipy.signal.find_peaks_cwt(observed_hist, range(1, 5)).tolist()
    print("Peaks", peaks)
    peaks =  merge_similar_elements([peaks[0]] + peaks, 3)
    print("Peaks", peaks)

    # https://stackoverflow.com/questions/25571260/scipy-signal-find-peaks-cwt-not-finding-the-peaks-accurately
    from scipy.ndimage.filters import gaussian_filter1d
    from scipy import signal
    # dataFiltered = gaussian_filter1d(observed_hist, sigma=5)
    # tMax = signal.argrelmax(dataFiltered)[0].tolist()
    # print('tMax:', tMax)

    pair_distance = [peaks[i] - peaks[i - 1] for i in range(1, len(peaks))]
    est_gauss_dist = np.median(pair_distance)
    print("Estimated peak distance", est_gauss_dist)
    est_offset = peaks[0]
    print("Estimated offset", est_offset)
    est_sigma = 3

    # if peaks[-1] > tMax[-1] + est_gauss_dist:
    #     tMax = tMax + [peaks[-1]]
    # if tMax[0] > peaks[0] + est_gauss_dist:
    #     tMax = [peaks[0]] + tMax
    # print('tMax:', tMax)
    #
    # add_extra = []
    # for i in range(len(tMax) - 1):
    #     if tMax[i + 1] - tMax[i] > est_gauss_dist * 2:
    #         add_extra.append((tMax[i] + tMax[i + 1]) / 2)
    #
    # for i, inner in enumerate(add_extra):
    #     tMax.append(inner)
    # peaks_ = sorted(tMax)
    # print('tMax:', peaks_)

    # numerical optimization
    #fr = fit_optimization(observed_hist, len(peaks), est_offset, est_gauss_dist, est_sigma)
    #distr_test = model_distribution(fr.x[0], fr.x[1], fr.x[2], fr.x[3:])
    #yy = [distr_test(x) for x in xx]

    # peaks = scipy.signal.find_peaks_cwt(yy, range(1, 5)).tolist()
    # if est_offset > est_gauss_dist//2:
    #     peaks = [est_offset-est_offset] + peaks

    #plt.plot(xx, yy, "r")
    # plt.show()
    #plt.savefig("test.pdf", format="pdf", bbox_inches="tight")

    if est_offset > 3:#est_offset > est_gauss_dist//2:
        final_peaks = [0] + [peaks[0]] + [i * est_gauss_dist for i in range(1, len(peaks) + 13)]
    else:
        final_peaks = [peaks[0]] + [i*est_gauss_dist for i in range(1, len(peaks) + 13)]


    # if est_offset > 3:
    #     final_peaks = [0] + peaks + [i * est_gauss_dist for i in range(len(peaks), len(peaks) + 3)]
    # else:
    #     final_peaks = peaks + [i * est_gauss_dist for i in range(len(peaks), len(peaks) + 3)]

    #peaks = scipy.signal.find_peaks_cwt(yy, range(1, 5)).tolist()
    #print("optimized peaks", peaks)

    # Autocorrelation method
    corr = scipy.signal.correlate(observed_hist, observed_hist, mode="full")
    corr = corr[corr.size // 2:]  # autocorrelation, only positive shift, so getting right half of the array
    first_min = scipy.signal.argrelmin(corr)[0][0]
    corr_max = np.argmax(corr[first_min:]) + first_min
    print("First minimum", first_min)
    print("Max correlation peak", corr_max)

    # testing if 1/2 of the highest peak is also a peak. compare it with 1/4 and 3/4 (which should be valleys)
    half_peak = corr_max // 2
    valley_left, valley_right = corr_max // 4, corr_max * 3 // 4
    is_half_peak = corr[half_peak] > max(corr[valley_left], corr[valley_right])
    # print(half_peak, valley_left, valley_right)
    print("Half peak:", is_half_peak)

    if not is_half_peak:
        single_copy_cov = corr_max
    else:
        single_copy_cov = half_peak
    print("Estimated single copy coverage:", single_copy_cov)

    last_copy_state = int(max(observed) // single_copy_cov + 1)
    final_peaks = [i * single_copy_cov for i in range(0, last_copy_state)]

    final_peaks_subclonal = [single_copy_cov//2] + [i * single_copy_cov+(single_copy_cov//2) for i in range(1, last_copy_state)]

    return final_peaks, final_peaks_subclonal, np.arange(0, 500), observed_hist