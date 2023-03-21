import numpy as np
#import talib as ta

#def talib_sma(array, period):
#    return ta.MA(array, period)
def convolve_sma(x, N):
    return np.convolve(x, np.ones((N,))/N)[(N-1):]
def smoothing(unphased_reads_values, haplotype_1_values, haplotype_2_values, conv_window_size):
    return convolve_sma(unphased_reads_values, conv_window_size), \
        convolve_sma(haplotype_1_values, conv_window_size), \
        convolve_sma(haplotype_2_values, conv_window_size)
