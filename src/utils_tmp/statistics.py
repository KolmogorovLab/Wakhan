import numpy as np


def remove_outliers_iqr(data, threshold=5):
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    lower_bound = q1 - threshold * iqr
    upper_bound = q3 + threshold * iqr
    outliers_removed = data[(data >= lower_bound) & (data <= upper_bound)]
    return outliers_removed


def weighted_means(vals, weights):
    if not len(vals) == len(weights) or sum(weights) == 0:
        return 0
    weighted_vals = []
    vals_n_weights = [(vals[i], weights[i]) for i in range(0, len(weights))]
    for tup in vals_n_weights:
        weighted_vals.append(round(tup[0] * tup[1], 3))
    return sum(weighted_vals)/sum(weights)


def add_intermediaries(numbers, max_difference):

    if len(numbers) < 2:
        return numbers

    result = [numbers[0]]
    for i in range(1, len(numbers)):
        current_diff = numbers[i] - result[-1]
        if current_diff > max_difference:
            num_intermediaries = current_diff // max_difference
            for j in range(1, num_intermediaries + 1):
                result.append(result[-1] + j + max_difference)
        result.append(numbers[i])

    return result


def adjust_extreme_outliers(hp_data):
    for j, val in enumerate(hp_data):
        if val > 500 and j > 0 and j < len(hp_data):
            hp_data[j] = hp_data[j-1]
        elif val > 500 and j == 0:
            hp_data[j] = hp_data[j+1]
    return hp_data


def find_peak_median_without_outliers(data):
    # from scipy.signal import find_peaks
    # peaks, _ = find_peaks(data)
    # peak_values = [data[i] for i in peaks]
    # peak_index = peaks[np.argmax(peak_values)]
    # peak_value = data[peak_index]

    if len(data) == 0: #second time almost big bp boundries based PS will be one and will be set to 0, change condition
        return 0
    else:
        return statistics.median(data)


    #if len(data) < 5:
    #    return statistics.median(data)

    # 2. Remove outliers using IQR method
    # q1 = np.percentile(data, 25)
    # q3 = np.percentile(data, 75)
    # iqr = q3 - q1
    # lower_bound = q1 - 1.5 * iqr
    # upper_bound = q3 + 1.5 * iqr
    # filtered_data = data[(data >= lower_bound) & (data <= upper_bound)]
    def remove_outliers_modified_z(data, peak_value, threshold=3.5):
        data = np.array(data)
        median = peak_value #np.median(data)
        mad = np.median(np.abs(data - median))  # Median Absolute Deviation
        if mad == 0:
            return data  # No variation
        modified_z_scores = 0.6745 * (data - median) / mad
        return data[np.abs(modified_z_scores) <= threshold]

    filtered_data = remove_outliers_modified_z(data, peak_value)
    # 3. Compute median of filtered data
    return np.median(filtered_data)
