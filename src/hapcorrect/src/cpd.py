import numpy as np
import ruptures as rpt
from random import randint
import statistics

def cpd_positions_means(haplotype1_means, haplotype2_means, arguments):
    snps_haplotype1_mean, snps_haplotype1_len, snps_haplotype1_pos = change_point_detection_algo(arguments['bin_size'],
                                                                                                 haplotype1_means)
    hp1_pos_start = []
    hp1_pos_end = []
    for i in range(len(snps_haplotype1_pos) - 1):
        hp1_pos_start.append(snps_haplotype1_pos[i] if snps_haplotype1_pos[i] < 1 else snps_haplotype1_pos[i] + 1)
        hp1_pos_end.append(snps_haplotype1_pos[i + 1])

    snps_haplotype2_mean, snps_haplotype2_len, snps_haplotype2_pos = change_point_detection_algo(arguments['bin_size'],
                                                                                                 haplotype2_means)
    hp2_pos_start = []
    hp2_pos_end = []
    for i in range(len(snps_haplotype2_pos) - 1):
        hp2_pos_start.append(snps_haplotype2_pos[i] if snps_haplotype2_pos[i] < 1 else snps_haplotype2_pos[i] + 1)
        hp2_pos_end.append(snps_haplotype2_pos[i + 1])

    return snps_haplotype1_mean, hp1_pos_start, hp1_pos_end, snps_haplotype2_mean, hp2_pos_start, hp2_pos_end

def change_point_detection_algo(bin_size, haplotype_means):
    data = np.array(haplotype_means, dtype='uint8')  # numpy.clip(haplotype1_means, a_min=1, a_max=300)
    algo = rpt.Pelt(model="rbf", jump=25).fit(data)
    result = algo.predict(pen=10)
    change_points = [i for i in result if i <= len(data)]
    strong_candidates = []
    snps_haplotype_mean = []
    snps_haplotype_pos = []
    snps_haplotype_len = []
    start = 0
    snps_haplotype_pos.append(0)
    for index, point in enumerate(change_points):
        sub_list = haplotype_means[start:point]
        if len(sub_list) < 25:
            continue
        else:
            count = 0
            sub_list = [int(i) for i in sub_list]
            mode_hp1 = max(set(sub_list), key=sub_list.count)
            means_hp1 = statistics.mean(sub_list)
            means_hp1 = mode_hp1 if mode_hp1 > means_hp1 else means_hp1
            means_hp1 = statistics.median(sub_list)
            # means_hp1 = statistics.mean(sub_list)
            for i in range(len(sub_list)):
                if means_hp1 - 20 <= sub_list[i] <= means_hp1 + 20:
                    count += 1
            if count > len(sub_list) // 3:
                snps_haplotype_mean.append(means_hp1)
                snps_haplotype_len.append(len(sub_list))

            if count > len(sub_list) / 1.1:
                strong_candidates.append(int(means_hp1))

            for i in range(len(sub_list)):
                if sub_list[i] >= means_hp1 + 10 or sub_list[i] <= means_hp1 - 10:
                    sub_list[i] = randint(int(means_hp1) - 10, int(means_hp1) + 10)
            haplotype_means[start:point] = [sub_list[i] for i in range(len(sub_list))]

        start = point + 1
        snps_haplotype_pos.append(point * bin_size)

    return snps_haplotype_mean, snps_haplotype_len, snps_haplotype_pos