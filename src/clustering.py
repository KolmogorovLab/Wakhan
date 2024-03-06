import ruptures as rpt
import numpy as np
import statistics
import pandas as pd
from random import randint
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
def cluster_coverage_data(cpd_means_input, collective):

    #slices = list(slice_lists(lambda x, y: y - x >= 2, sorted(cpd_means_input)))
    #cpd_means_input = [sum(sub_list) / len(sub_list) for sub_list in slices]

    if len(cpd_means_input) > 3:
        cpd_centers = []
        cpd_stdev = []
        cpd_score = []
        cpd_means_input = np.array(cpd_means_input, dtype='int')

        MAX = 10
        if len(cpd_means_input) < 10:
            MAX = len(cpd_means_input)
        cpd_means_input = cpd_means_input.reshape(-1, 1)

        for i in range(3, MAX, 1):
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

        if collective:
            centers = slice_list_sums(cpd_means_input)
            stdev = stdev_calc(centers)

        #TODO - removing nearby clusters
        # if collective:
        #     slices = list(slice_lists(lambda x, y: y - x > 2, sorted(centers[1:])))
        #     if len(slices) < len(centers)-1:
        #         zero_center = [centers[0]]
        #         centers = [sum(sub_list) / len(sub_list) for sub_list in slices]
        #         centers = zero_center + centers
        #         stdev = stdev_calc(center)

    else:
        slices = list(slice_lists(lambda x, y: y - x > 2, sorted(cpd_means_input)))
        centers = [sum(sub_list) / len(sub_list) for sub_list in slices]
        #centers = slice_list_sums(sorted(cpd_means_input))
        print(centers)
        stdev = stdev_calc(centers)

    return centers, stdev

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
            stdev.append(1)
        elif i == len(centers):
            stdev.append(5)
        else:
           stdev.append(5)
    return stdev

def slice_list_sums(input):
    res, last = [[]], None
    for x in sorted(input):
        if last is None or abs(last - x) <= 2.5:
            res[-1].append(x)
        else:
            res.append([x])
        last = x
    first = [res[0][0]]
    output = first + [sum(sub_list) / len(sub_list) for sub_list in res[1:]]
    output = [x for xs in [i.tolist() for i in output] for x in xs]

    return output