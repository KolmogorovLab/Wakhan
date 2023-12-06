import ruptures as rpt
import numpy as np
import statistics
from random import randint
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
def cluster_coverage_data(cpd_means_input):

    slices = list(slice_lists(lambda x, y: y - x > 2, sorted(cpd_means_input)))
    if len(slices) > 3:
        cpd_centers = []
        cpd_stdev = []
        cpd_score = []
        cpd_means_input = np.array(cpd_means_input, dtype='int')

        MAX = 12
        if len(cpd_means_input) < 12:
            MAX = len(slices)
        cpd_means_input = cpd_means_input.reshape(-1, 1)

        for i in range(3, MAX+1, 1):
            kmeans = KMeans(n_clusters=i, random_state=0, n_init="auto").fit(cpd_means_input)
            center = sorted(list(np.concatenate(kmeans.cluster_centers_)))
            s_score = silhouette_score(cpd_means_input, kmeans.labels_, metric="euclidean")
            sd = stdev_calc(center)
            cpd_centers.append(center)
            cpd_score.append(s_score)
            cpd_stdev.append(sd)
            print(s_score, [round(c, 3) for c in center])

        centers = cpd_centers[np.argmax(cpd_score)]
        stdev = cpd_stdev[np.argmax(cpd_score)]
    else:
        centers = [sum(sub_list) / len(sub_list) for sub_list in slices]
        print(centers)
        stdev = stdev_calc(centers)

    return centers, stdev
def change_point_detection_means(arguments, df_chrom):
    if arguments['without_phasing']:
        haplotype1_means = df_chrom.coverage.values.tolist()
        haplotype2_means = df_chrom.coverage.values.tolist()
    else:
        haplotype1_means = df_chrom.hp1.values.tolist()
        haplotype2_means = df_chrom.hp2.values.tolist()

    data = np.array(haplotype1_means, dtype='uint8')  # numpy.clip(haplotype1_means, a_min=1, a_max=300)
    algo = rpt.Pelt(model="rbf").fit(data)
    result = algo.predict(pen=10)
    change_points = [i for i in result if i <= len(data)]
    strong_candidates = []
    snps_haplotype1_mean = []
    snps_haplotype1_pos = []
    snps_haplotype1_len = []
    start = 0
    snps_haplotype1_pos.append(0)
    for index, point in enumerate(change_points):
        sub_list = haplotype1_means[start:point]
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
                snps_haplotype1_mean.append(means_hp1)
                snps_haplotype1_len.append(len(sub_list))

            if count > len(sub_list) / 1.1:
                strong_candidates.append(int(means_hp1))

            for i in range(len(sub_list)):
                if sub_list[i] >= means_hp1 + 10 or sub_list[i] <= means_hp1 - 10:
                    sub_list[i] = randint(int(means_hp1) - 10, int(means_hp1) + 10)
            haplotype1_means[start:point] = [sub_list[i] for i in range(len(sub_list))]

        start = point + 1
        snps_haplotype1_pos.append(point * 50000)
    ############################################################
    data = np.array(haplotype2_means, dtype='uint8')  # numpy.clip(haplotype2_means, a_min=1, a_max=300)
    algo = rpt.Pelt(model="rbf").fit(data)
    result = algo.predict(pen=10)
    change_points = [i for i in result if i <= len(data)]

    snps_haplotype2_mean = []
    snps_haplotype2_pos = []
    snps_haplotype2_len = []

    start = 0
    snps_haplotype2_pos.append(0)
    for index, point in enumerate(change_points):
        sub_list = haplotype2_means[start:point]
        if len(sub_list) < 25:
            continue
        else:
            count = 0
            sub_list = [int(i) for i in sub_list]
            mode_hp2 = max(set(sub_list), key=sub_list.count)
            means_hp2 = statistics.mean(sub_list)
            means_hp2 = mode_hp2 if mode_hp2 > means_hp2 else means_hp2
            means_hp2 = statistics.median(sub_list)
            for i in range(len(sub_list)):
                if means_hp2 - 20 <= sub_list[i] <= means_hp2 + 20:
                    count += 1
            if count > len(sub_list) // 3:
                snps_haplotype2_mean.append(means_hp2)
                snps_haplotype2_len.append(len(sub_list))

            if count > len(sub_list) / 1.1:
                strong_candidates.append(int(means_hp2))

            for i in range(len(sub_list)):
                if sub_list[i] >= means_hp2 + 10 or sub_list[i] <= means_hp2 - 10:
                    sub_list[i] = randint(int(means_hp2) - 10, int(means_hp2) + 10)
            haplotype2_means[start:point] = [sub_list[i] for i in range(len(sub_list))]

        start = point + 1
        snps_haplotype2_pos.append(point * 50000)
    #return snps_haplotype1_mean, snps_haplotype1_len, snps_haplotype2_mean, snps_haplotype2_len, haplotype1_means, haplotype2_means, strong_candidates
    return snps_haplotype1_mean + snps_haplotype2_mean

def slice_lists(predicate, iterable):
  i, x, size = 0, 0, len(iterable)
  while i < size-1:
    if predicate(iterable[i], iterable[i+1]):
      yield iterable[x:i+1]
      x = i + 1
    i += 1
  yield iterable[x:size]

def stdev_calc(centers):
    if len(centers) == 1:
        centers = [0.25] + centers
    stdev = []
    for i in range(len(centers)):
        if i == 0:
            stdev.append(5)
        elif i == len(centers):
            stdev.append(25)
        else:
           stdev.append(20)
    return stdev
