import numpy as np
#import talib as ta
import logging

logger = logging.getLogger()

#def talib_sma(array, period):
#    return ta.MA(array, period)

def smoothing(haplotype_1_values, haplotype_2_values, unphased_reads_values = None, conv_window_size = 5):
    if not unphased_reads_values == None:
        unphased_reads_values = smooth_triangle(unphased_reads_values, conv_window_size)

    haplotype_1_values = smooth_triangle(haplotype_1_values, conv_window_size)
    haplotype_2_values = smooth_triangle(haplotype_2_values, conv_window_size)
    return haplotype_1_values, haplotype_2_values, unphased_reads_values

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