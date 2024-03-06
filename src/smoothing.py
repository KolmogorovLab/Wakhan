import numpy as np
#import talib as ta

#def talib_sma(array, period):
#    return ta.MA(array, period)
def convolve_sma(x, N):
    return np.convolve(x, np.ones((N,))/N)[(N-1):]

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
    # Handle boundaries
    smoothed=[smoothed[0]]*int(degree + degree/2) + smoothed
    while len(smoothed) < len(data):
        smoothed.append(smoothed[-1])
    return smoothed

def smooth_gaussian(list, degree):
    window = degree*2-1
    import numpy
    weight = numpy.array([1.0]*window)
    weightGauss = []
    for i in range(window):
        i = i-degree+1
        frac = i/float(window)
        gauss = 1/(numpy.exp((4*(frac))**2))
        weightGauss.append(gauss)
    weight = numpy.array(weightGauss)*weight
    smoothed = [0.0]*(len(list)-window)
    for i in range(len(smoothed)):
        smoothed[i] = sum(numpy.array(list[i:i+window])*weight)/sum(weight)
    return smoothed