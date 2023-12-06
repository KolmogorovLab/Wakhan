"""Segmentation by Hidden Markov Model."""
import collections
import logging
import statistics

import numpy as np
import pandas as pd
import pomegranate as pom
import scipy.special

from ..cnary import CopyNumArray as CNA
from ..descriptives import biweight_midvariance
from ..segfilters import squash_by_groups
from ..cluster import kmeans_clustering, hmm_model_select_hatchet, hmm_pome_test, cpd_mean_hps, hmm_validation_custom, add_scatter_trace_coverage



from pomegranate import State, NormalDistribution, IndependentComponentsDistribution, DiscreteDistribution
from pomegranate import HiddenMarkovModel as HMM

def segment_hmm(depth_values_hp1, depth_values_hp2, snps_cpd_means_input, arguments, cnarr, method, window=None, variants=None, processes=1):
    """Segment bins by Hidden Markov Model.

    Use Viterbi method to infer copy number segments from sequential data.

    With b-allele frequencies ('baf' column in `cnarr`), jointly segment
    log-ratios and b-allele frequencies across a chromosome.

    Parameters
    ----------
    cnarr : CopyNumArray
        The bin-level data to segment.
    method : string
        One of 'hmm' (3 states, flexible means), 'hmm-tumor' (5 states, flexible
        means), 'hmm-germline' (3 states, fixed means).

    Results
    -------
    segarr : CopyNumArray
        The segmented data.
    """
    # NB: Incorporate weights into smoothed log2 estimates
    # (Useful kludge until weighted HMM is in place)
    orig_log2 = cnarr["log2"].values.copy()
    #cnarr["log2"] = cnarr.smooth_log2()  # window)
    cnarr.data = cnarr.data.reset_index(drop=True)
    observations = observations_matrix(cnarr, snps_cpd_means_input)#as_observation_matrix(cnarr)

    logging.info("Building model from observations")
    model, centers, stdev, cnarr, snps_cpd_means = hmm_get_model(depth_values_hp1, depth_values_hp2, snps_cpd_means_input, arguments, cnarr, method, processes)

    #mus = np.ravel(model.means_)
    #sigmas = np.ravel(np.sqrt([np.diag(c) for c in model.covars_]))
    #P = model.transmat_
    #print(mus, sigmas, P)

    logging.info("Predicting states from model")
    states = np.concatenate(
        [np.array(model.predict(obs, algorithm="map")) for obs in observations]
    )
    #########################################################
    #from sklearn.metrics import silhouette_score

    #depth_values_hp1 = np.array(depth_values_hp1, dtype='int')
    #depth_values_hp2 = np.array(depth_values_hp2, dtype='int')
    #depth_values_hp1 = depth_values_hp1.reshape(-1, 1)
    #depth_values_hp2 = depth_values_hp2.reshape(-1, 1)
    #score = silhouette_score(np.concatenate([depth_values_hp1, depth_values_hp2]), states, metric="euclidean")
    #print(score)
    #########################################################
    logging.info("Done, now finalizing")
    logging.debug("Model states: %s", model.states)
    logging.debug("Predicted states: %s", states[:100])
    logging.debug(str(collections.Counter(states)))
    logging.debug("Observations: %s", observations[0][:100])
    logging.debug("Edges: %s", model.edges)

    values = cnarr.as_dataframe(cnarr.data)
    values.data.reset_index(drop=True, inplace=True)
    if snps_cpd_means_input:
        half_values = values.data[(values.data['chromosome'] == "chr1") & (values.data['start'] == 0)].index[1]
    else:
        half_values = values.data[(values.data['start'] == 0)].index[1]

    segs = []
    start = 0
    haplotype_dfs = pd.DataFrame()
    for i in range(2): #range(0, len(values), 49592):#len(values) // 2):
        state = states[start:half_values]
        cnarray = values.data.iloc[start:half_values, ]
        start = half_values
        half_values += (len(values) - half_values)#half_values

        meta = {"sample_id": 'sample'}
        cnarray = CNA(cnarray, meta)

        # Merge adjacent bins with the same state to create segments
        #cnarr["log2"] = orig_log2
        cnarray["probes"] = 1
        segarr = squash_by_groups(
            cnarray, pd.Series(state, index=cnarray.data.index), by_arm=True
        )
        if not (segarr.start < segarr.end).all():
            bad_segs = segarr[segarr.start >= segarr.end]
            logging.warning("Bad segments:\n%s", bad_segs.data)

        haplotype_df = segarr.data
        haplotype_df['haplotype'] = i
        haplotype_dfs = haplotype_dfs.append(haplotype_df)

        mean = []
        df = cnarray.data.assign(group=state)
        for i, row in segarr.data.iterrows():
            size = len(segarr.data.index)
            #if row.chromosome == 'chr15':
            #    print("here")
            if (row.end - row.start < 1000000) and (i > 0 and i < size-1):
                if segarr.data.at[i - 1, 'state'] == segarr.data.at[i + 1, 'state'] and segarr.data.at[i, 'chromosome'] == segarr.data.at[i - 1, 'chromosome']:
                    segarr.data.at[i, 'state'] = segarr.data.at[i+1, 'state']
                elif (not segarr.data.at[i - 1, 'state'] == segarr.data.at[i + 1, 'state']) and segarr.data.at[i, 'chromosome'] == segarr.data.at[i - 1, 'chromosome']:
                        segarr.data.at[i, 'state'] = segarr.data.at[i - 1, 'state']

        for i, val in enumerate(np.unique(df['group'])):
            #mean.append(np.mean(df[df['group'] == i]['depth']))
            #segarr['state'] = np.where(segarr['state'] == i, mean[i], segarr['state'])
            segarr['state'] = np.where(segarr['state'] == val, centers[val], segarr['state'])  #TODO Bug, when same center and state values
        segs.append(segarr)

    header = ['chromosome', 'start', 'end', 'state', 'haplotype']
    haplotype_dfs.to_csv('data/copynumbers_segments.csv', sep='\t', columns=header, index=False)

    return segs, states, centers, stdev, cnarr, snps_cpd_means

def squash_regions(df):

    return pd.DataFrame(df)

def slice_lists(predicate, iterable):
  i, x, size = 0, 0, len(iterable)
  while i < size-1:
    if predicate(iterable[i], iterable[i+1]):
      yield iterable[x:i+1]
      x = i + 1
    i += 1
  yield iterable[x:size]

def closest_value(input_list, input_value):
  arr = np.asarray(input_list)
  i = (np.abs(arr - input_value)).argmin()
  return arr[i].tolist()

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

def stdev_updated(centers):
    if len(centers) > 2:
        stdev=[]
        for i in range(len(centers)):
            if not i == 0 and not i == len(centers) - 1:
                stdev.append(min(centers[i] - centers[i - 1], centers[i + 1] - centers[i]) / 2)
        stdev = [stdev[0]] + stdev + [stdev[-1]]
    else:
        stdev = stdev_calc(centers)
    return stdev

def clusters_manual(input, iter):
    for i in range(iter):
        slices = list(slice_lists(lambda x, y: y - x > (i + 1) * 2, sorted(input)))
        centers = [sum(sub_list) / len(sub_list) for sub_list in slices]
    return centers

def optimize_manual_clusters(input, cnarr, X, snps_cpd_means_input):
    manual_centers = clusters_manual(input, 2)
    data = np.array(manual_centers, dtype='int')
    data = data.reshape(-1, 1)
    KMeans_centers = []
    from sklearn.cluster import KMeans
    for i in range(3, 13, 1):
        kmeans = KMeans(n_clusters=i, random_state=0, n_init="auto").fit(data)
        KMeans_centers.append([round(j, 3) for j in sorted(list(np.concatenate(kmeans.cluster_centers_)))])

    cpd_centers_test = []
    cpd_score_test = []
    cpd_stdev_test = []
    for i in range(len(KMeans_centers)):
        ct = sorted(KMeans_centers[i])
        sd = stdev_updated(KMeans_centers[i])
        centers_test, stdev_test, score_test = hmm_pome_test(cnarr, X, ct, sd, snps_cpd_means_input)
        cpd_centers_test.append(centers_test)
        cpd_score_test.append(score_test)
        cpd_stdev_test.append(stdev_test)
        print(score_test, [round(c, 3) for c in centers_test])

    for i in range(len(KMeans_centers)):
        print(cpd_score_test[i], [round(c, 3) for c in cpd_centers_test[i]])

    centers = cpd_centers_test[np.argmax(cpd_score_test)]

    return sorted(centers)

def hmm_get_model(depth_values_hp1, depth_values_hp2, snps_cpd_means_input, arguments, cnarr, method, processes):
    """

    Parameters
    ----------
    cnarr : CopyNumArray
        The normalized bin-level values to be segmented.
    method : string
        One of 'hmm', 'hmm-tumor', 'hmm-germline'.
    processes : int
        Number of parallel jobs to run.

    Returns
    -------
    model :
        A pomegranate HiddenMarkovModel trained on the given dataset.
    """


    #Standard Normal Distribution
    #Mean (μ, mu) -> 0, variance -> σ2 -> 1 (std deviation (σ, sigma))

    # Estimate standard deviation from the full distribution, robustly
    #stdev = biweight_midvariance(np.concatenate(observations), initial=0)

    #depth_values_hp1 = np.clip(depth_values_hp1, a_min=1, a_max=250)
    #depth_values_hp2 = np.clip(depth_values_hp2, a_min=1, a_max=250)
    #depth_values_hp1 = depth_values_hp1.astype(int)
    #depth_values_hp2 = depth_values_hp2.astype(int)

    hp1 = [i for i in depth_values_hp1 if i >= 0]
    hp2 = [i for i in depth_values_hp2 if i >= 0]
    #depth_values_hp1 = [i for i in depth_values_hp1 if i < 600]
    #depth_values_hp2 = [i for i in depth_values_hp2 if i < 600]

    depth_values_hp1 = np.array(hp1, dtype='int')
    depth_values_hp2 = np.array(hp2, dtype='int')

    depth_values_hp1 = depth_values_hp1.reshape(-1, 1)
    depth_values_hp2 = depth_values_hp2.reshape(-1, 1)
    X = np.concatenate([depth_values_hp1, depth_values_hp2])

    centerss = []
    stdevs = []
    scores = []
    if snps_cpd_means_input:
        print(snps_cpd_means_input)
        input = snps_cpd_means_input[0]
        weights = snps_cpd_means_input[1]
        snps_haplotype1_len = []
        snps_haplotype2_len = []
    else:
        snps_haplotype1_mean, snps_haplotype1_len, snps_haplotype2_mean, snps_haplotype2_len, hp1, hp2, strong_candidates  = cpd_mean_hps(hp1, hp2)
        input = snps_haplotype1_mean + snps_haplotype2_mean
        weights = snps_haplotype1_len + snps_haplotype2_len
        print(input)
        print(weights)
        print(strong_candidates)

    slices = list(slice_lists(lambda x, y: y - x > 2, sorted(input)))
    if len(slices) > 2:

        centers_kmeans, sd_kmeans, score_kmeans = kmeans_clustering(hp1 + hp2, None)
        print("centers_kmeans:", centers_kmeans)


        #centers_kmeans, sd_kmeans, score_kmeans = kmeans_clustering(input, None)
        #print("input:", score_kmeans, centers_kmeans, sd_kmeans)
        from sklearn.cluster import KMeans
        from sklearn.metrics import silhouette_score

        cpd_centers = []
        cpd_stdev = []
        cpd_score = []
        input_updated = np.clip(sorted(input), a_min=0, a_max=600)
        input_weight = np.clip(sorted(weights), a_min=0, a_max=12000)
        MAX = 12
        if len(input_updated) < 12:
            MAX = len(slices)
        input_updated = input_updated.reshape(-1, 1)
        #input_weight = input_weight.reshape(-1, 1)
        for i in range(3, MAX+1, 1):
            kmeans = KMeans(n_clusters=i, random_state=0, n_init="auto").fit(input_updated, sample_weight=input_weight)
            cen = sorted(list(np.concatenate(kmeans.cluster_centers_)))
            sr = silhouette_score(input_updated, kmeans.labels_, metric="euclidean")
            st = stdev_calc(cen)
            cpd_centers.append(cen)
            cpd_score.append(sr)
            cpd_stdev.append(st)
            print(sr, [round(c, 3) for c in cen])

        if cpd_score:
           centers = cpd_centers[np.argmax(cpd_score)]
           #sd = cpd_stdev[np.argmax(cpd_score)]

        manual_centers = clusters_manual(input, 2)
        print(manual_centers)

        # print('cpd_centers:', cpd_centers)
        # if snps_cpd_means_input:
        #     for i in range(2):
        #         for i, center in enumerate(cpd_centers):
        #             len_ = len(center)
        #             print(center)
        #             for j in range(len_-1):
        #                 if center[j+1] - center[j] < 9:
        #                     del cpd_centers[i][j]
        #                     del cpd_centers[i][j]
        #                     cpd_centers[i].insert(j+1, (statistics.mean([center[j+1],center[j]])))
        #                     break
        #             print(cpd_centers[i])
        #     print('cpd_centers:', cpd_centers)
        #
        # for i, center in enumerate(cpd_centers): #TODO penalty if diff
        #     final_centers = center
        #     for j in range(len(center)-1):
        #         if center[j+1] - center[j] < 20:
        #             final_centers = cpd_centers[i-1]
        #             break
        #     if final_centers == cpd_centers[i-1]:
        #         break
        #
        # centers = final_centers
    else:
        centers = [sum(sub_list) / len(sub_list) for sub_list in slices]

    if snps_cpd_means_input:
        centers = [2.503, 11.291, 22.451, 44.206]

    #if snps_cpd_means_input:
    #    centers = [1.0607034632034598, 21.010000000000005, 45.210520833333334, 58.67666666666666, 86.22594202898551, 120.01, 138.51, 164.01, 207.51]

    #if snps_cpd_means_input:
    #    centers = optimize_manual_clusters(input, cnarr, X, snps_cpd_means_input)

    # for i in range(2):
    #     slices = list(slice_lists(lambda x, y: y - x > (i + 1) * 5, sorted(input)))
    #     cpd_centers = [sum(sub_list) / len(sub_list) for sub_list in slices]
    #     stdev = stdev_calc(cpd_centers)
    #     print(cpd_centers)
    #     depth_values_hp1 = np.array(hp1, dtype='int')
    #     depth_values_hp2 = np.array(hp2, dtype='int')
    #
    #     depth_values_hp1 = depth_values_hp1.reshape(-1, 1)
    #     depth_values_hp2 = depth_values_hp2.reshape(-1, 1)
    #     X = np.concatenate([depth_values_hp1, depth_values_hp2])
    #     lengths = [len(depth_values_hp1), len(depth_values_hp2)]
    #     cnarr_copy = cnarr.copy()
    #     cnarr_copy['depth'] = hp1 + hp2
    #     cnarr_copy['log2'] = hp1 + hp2
    #     centers, sd, score = hmm_pome_test(cnarr_copy, X, cpd_centers, stdev)
    #     print(score, centers, sd)
    #     centerss.append(centers)
    #     stdevs.append(sd)
    #     scores.append(score)
    #
    # #print(centerss, stdevs, scores)
    # if scores:
    #     centers = centerss[np.argmax(scores)]
    #     stdev = stdevs[np.argmax(scores)]

    #centers = [1.1022376543209877, 19.764044943820224, 62.720703914377665, 93.91878722252301, 120.2457627118644]
    #stdev = [2.087982469231899, 7.678219597001331, 18.48541231592829, 12.177732432494645, 10.185079869112837]




    #centers = [1.1846644169478815, 48.505, 82.5127427184466, 114.8906439854192, 142, 190.78205128205127]
    #stdev = [1.7935393331214273, 23.146584758121712, 10.34797518399463, 12.003251686298688, 20, 15.05770774115513]

    #centers = [0.9621301775147929, 21.55, 51.21556886227545, 73.13020489094514, 94.68897637795276, 116.9951923076923]
    #stdev = [5, 10, 10.317790345277219, 13.74241853604731, 12.261716355645715, 12.01086887228625]

    # for i in range(len(centers) - 1):
    #     for j in range(len(hp1)-2):
    #         if j>0 and (centers[i] + 10 <= hp1[j] <= centers[i] - 10) and (centers[i] + 10 <= hp1[j-1] <= centers[i] - 10) and (centers[i] + 10 <= hp1[j+1] <= centers[i] - 10):
    # b3reeeeeeeeeeeeee             hp1[j] = centers[i]
    #         if j>0 and (centers[i] + 10 <= hp2[j] <= centers[i] - 10) and (centers[i] + 10 <= hp2[j-1] <= centers[i] - 10) and (centers[i] + 10 <= hp2[j+1] <= centers[i] - 10):
    #             hp2[j] = centers[i]

    from random import randint
    # hp1 = [randint(35, 55) for p in range(0, 805)] + [randint(100, 150) for p in range(0, 500)] + [randint(32, 45) for p in range(0, 250)] + [randint(10, 30) for p in range(0, 500)] + [randint(35, 55) for p in range(0, 250)] + [randint(80, 100) for p in range(0, 500)] + [randint(80, 100) for p in range(0, 1000)]
    # hp2 = [randint(35, 55) for p in range(0, 250)] + [randint(10, 30) for p in range(0, 500)] + [randint(45, 79) for p in range(0, 805)] + [randint(80, 100) for p in range(0, 500)] + [randint(80, 100) for p in range(0, 1000)] + [randint(100, 150) for p in range(0, 500)] + [randint(35, 55) for p in range(0, 250)]
    # centers = [20, 45, 90, 125]
    # stdev = [10, 10, 10, 25]
    #cnarr['depth'] = hp1 + hp2
    #cnarr['log2'] = hp1 + hp2

    # INITIAL = 2
    # for i in range(len(centers) - 1):
    #     for j in range(len(hp1)-1):
    #         if centers[i] + ((centers[i + 1] - centers[i]) // INITIAL) - INITIAL <= hp1[j] <= centers[i] + ((centers[i + 1] - centers[i]) // INITIAL) - INITIAL + 2.5:
    #             hp1[j] = centers[i]
    #         elif centers[i] + ((centers[i + 1] - centers[i]) // INITIAL) - INITIAL + 2.5 <= hp1[j] <= centers[i] + ((centers[i + 1] - centers[i]) // INITIAL) + INITIAL:
    #             hp1[j] = centers[i+1]
    #
    #     for j in range(len(hp2) - 1):
    #         if centers[i] + ((centers[i + 1] - centers[i]) // INITIAL) - INITIAL <= hp2[j] <= centers[i] + ((centers[i + 1] - centers[i]) // INITIAL) - INITIAL + 2.5:
    #             hp2[j] = centers[i]
    #         elif centers[i] + ((centers[i + 1] - centers[i]) // INITIAL) - INITIAL + 2.5 <= hp2[j] <= centers[i] + ((centers[i + 1] - centers[i]) // INITIAL) + INITIAL:
    #             hp2[j] = centers[i + 1]
    #
    # cnarr['depth'] = hp1 + hp2
    # cnarr['log2'] = hp1 + hp2

    # mids = []
    # for i in range(len(centers) - 1):
    #     mids.append(centers[i] + (centers[i + 1] - centers[i]) / 2)
    # for i in range(len(centers)):
    #     for j in range(len(hp1) - 1):
    #         if i == 0:
    #             if centers[0] + 5 <= hp1[j] <= mids[0]:
    #                 hp1[j] = randint(int(centers[0]), int(centers[0]) + 5)
    #             if centers[0] + 5 <= hp2[j] <= mids[0]:
    #                 hp2[j] = randint(int(centers[0]), int(centers[0]) + 5)
    #         elif i == len(centers)-1:
    #             if hp1[j] > centers[-1] +15:
    #                 hp1[j] = randint(int(centers[-1]), int(centers[-1]) + 15)
    #             if hp2[j] > centers[-1] +15:
    #                 hp2[j] = randint(int(centers[-1]), int(centers[-1]) + 15)


    #cnarr['depth'] = hp1 + hp2
    #cnarr['log2'] = hp1 + hp2

    #
    #
    #
    #
    # #centers, stdev  = hmm_model_select_hatchet(X, lengths, minK=2, maxK=10, tau=10e-6, tmat='free', decode_alg='viterbi', covar='diag', restarts=5, )
    # #print(centers, stdev)
    #
    # centers1, stdev1, score1 = kmeans_clustering(hp1 + hp2, None)
    # print(centers1, score1)
    # centers = centers1
    #
    # centers2, stdev2, score2 = kmeans_clustering(hp1 + hp2, len(centers1)+1)
    # print(centers2, score2)
    # if score2 > score1:
    #     centers = centers2
    #
    # centers3, stdev3, score3 = kmeans_clustering(hp1 + hp2, len(centers1)+2)
    # print(centers3, score3)
    # if score3 > score2 > score1:
    #     centers = centers3
    #
    # centers = sorted(centers)
    #
    # print("Kmeans: ", centers)
    # for i, val in enumerate(cpd_centers):
    #     closest  = closest_value(centers, val)
    #     if abs(closest - val) > 15:
    #         centers = centers + [val]
    # print("Kmeans: ", centers)
    #
    # centers = cpd_centers
    # for i in range(2):
    #     centers, stdev = zip(*sorted(zip(centers, stdev)))
    #     centers = list(centers)
    #     stdev = list(stdev)
    #
    #     change = []
    #     change_sd = []
    #     new = []
    #     new_sd = []
    #     for i in range(len(centers)-1):
    #         if centers[i+1] - centers[i] < 20:
    #             new.append(statistics.mean([centers[i+1],centers[i]]))
    #             new_sd.append(stdev[i + 1]+stdev[i])
    #             change.append(centers[i+1])
    #             change.append(centers[i])
    #             change_sd.append(stdev[i+1])
    #             change_sd.append(stdev[i])
    #     if new:
    #         centers = centers + new
    #         stdev = stdev + new_sd
    #         centers = [e for e in centers if e not in change]
    #         stdev = [e for e in stdev if e not in change_sd]
    #         centers, stdev = zip(*sorted(zip(centers, stdev)))
    #         centers = list(centers)
    #         stdev = list(stdev)
    #
    #

    #stdev = stdev_calc(centers)

    # print("Kmeans: ", centers)

    #centers = cpd_centers#[int(i) for i in centers]

    ####################################
    #centers = [3, 12, 27, 54, 81, 113]#1437
    #stdev = [1.5, 6, 13.5, 27.5, 40, 56]#1437

    #centers = [1.5, 30, 55, 103, 160, 220]#1937
    #stdev = [.75, 15, 27.5, 52.5, 80, 110]#1937

    #diff not > 10 between each pair
    #when more than one pairs start < 10 - stop


    #centers = [1.06164962, 14.74069362, 27.21298619, 55.05099708, 134.37462712, 78.00189026] #1437
    #centers = [1.37287895, 149.12153768, 54.32031415, 102.23449339, 23.6867503, 80.60718194, 126.76419107] #1937

    #centers = [7.212171659656513, 47.937524853177365, 88.1970515518751, 124.70888217695517, 168]
    #

    # centers = [int(i) for i in centers]
    # for i, strong in enumerate(strong_candidates):
    #     if not any(x in centers for x in range(strong-5, strong+5)):
    #         centers.append(strong)
    # centers = sorted(centers)

    stdev = stdev_updated(centers)
    if stdev[0] < 1:
        stdev[0] = 5
    #centers = [0] + centers
    #stdev = [5] + stdev

    #stdev = [5, 15, 15, 20, 20]

    #centers = [1, 22, 44, 84, 120, 165, 212]
    #stdev = [5, 5, 15, 20, 20, 20, 20]
    #centers = [2, 44, 84, 120, 180]
    #stdev = [20, 20, 20, 20, 20]

    #stdev, centers  = hmm_pome_test(cnarr, X, centers, stdev)
    #print(stdev, centers)

    # centers_change = []
    # stdev_change = []
    # for i in range(len(centers)-1):
    #     if centers[i+1] - centers[i] < stdev[i]:
    #         centers.append(statistics.mean([centers[i+1],centers[i]]))
    #         stdev.append(stdev[i + 1]+ stdev[i])
    #         centers_change.append(centers[i+1])
    #         centers_change.append(centers[i])
    #         stdev_change.append(stdev[i+1])
    #         stdev_change.append(stdev[i])
    #     if stdev[i] > 50:
    #         stdev[i] = 50
    # centers = sorted(centers)
    # centers = [e for e in centers if e not in centers_change]
    # stdev = [e for e in stdev if e not in stdev_change]

    #centers = cpd_centers
    #stdev = stdev_calc(centers)


    #centers = [round(cen, 3)+0.01 for cen in centers]

    #print([round(c, 3) for c in centers])
    print(centers , stdev)
    chroms = cnarr.as_dataframe(cnarr.data)
    chrom = chroms.chromosome.values.tolist()[0]
    # if chrom == 'chr5':
    #     centers = [2.1, 94, 156]
    #     stdev = [60, 30, 30]
    #
    # if chrom == 'chr16':
    #     centers = [2.1, 88, 122]
    #     stdev = [25, 20, 26]

    #centers = [1, 17, 50, 75, 106]
    #stdev = [3, 8, 12, 12, 16]

    #centers  = [1.0276372129233164, 25, 50.102932719954, 80.96294441662494, 110.03218390804598]
    #stdev = [5, 10, 20, 20, 20]


    ####################################
    state_names = []#["copy_1", "copy_2", "copy_3", "copy_4", "copy_5"]
    distributions = []
    for i in range(len(centers)):
        state_names.append("copy_"+str(i))
        distributions.append(pom.NormalDistribution(centers[i], stdev[i], frozen=False))

    # distributions = [
    #     pom.NormalDistribution(centers[0], stdev[0], frozen=False),
    #     pom.NormalDistribution(centers[1], stdev[1], frozen=False),
    #     pom.NormalDistribution(centers[2], stdev[2], frozen=False),
    #     pom.NormalDistribution(centers[3], stdev[3], frozen=False),
    #     pom.NormalDistribution(centers[4], stdev[4], frozen=False),
    # ]

    assert method in ("hmm-tumor", "hmm-germline", "hmm")
    observations = observations_matrix(cnarr, snps_cpd_means_input)#as_observation_matrix(cnarr.autosomes())

    n_states = len(distributions)
    # Starts -- prefer neutral
    binom_coefs = scipy.special.binom(n_states - 1, range(n_states))
    start_probabilities = binom_coefs / binom_coefs.sum()

    # Prefer to keep the current state in each transition
    # All other transitions are equally likely, to start
    transition_matrix = (
        np.identity(n_states) * 100 + np.ones((n_states, n_states)) / n_states
    )
    print(transition_matrix)
    print(start_probabilities)

    model = pom.HiddenMarkovModel.from_matrix(
        transition_matrix,
        distributions,
        start_probabilities,
        state_names=state_names,
        name=method,
    )
    model.fit(
        sequences=observations,
        weights=[len(obs) for obs in observations],
        distribution_inertia=0.9,  # Allow updating dists, but slowly
        edge_inertia=0.05,
        # lr_decay=.75,
        pseudocount=5,
        use_pseudocount=True,
        max_iterations=100000,
        n_jobs=processes,
        verbose=False,
    )

    return model, centers, stdev, cnarr, [input, snps_haplotype1_len + snps_haplotype2_len]


def as_observation_matrix(cnarr, variants=None):
    """Extract HMM fitting values from `cnarr`.

    For each chromosome arm, extract log2 ratios as a numpy array.

    Future: If VCF of variants is given, or 'baf' column has already been
    added to `cnarr` from the same, then the BAF values are a second row/column
    in each numpy array.

    Returns: List of numpy.ndarray, one per chromosome arm.
    """
    # TODO incorporate weights -- currently handled by smoothing
    # TODO incorporate inter-bin distances
    observations = [arm.log2.values for _c, arm in cnarr.by_arm()]
    return observations


def variants_in_segment(varr, segment, min_variants=50):
    if len(varr) > min_variants:
        observations = varr.mirrored_baf(above_half=True)
        state_names = ["neutral", "alt"]
        distributions = [
            pom.NormalDistribution(0.5, 0.1, frozen=True),
            pom.NormalDistribution(0.67, 0.1, frozen=True),
        ]
        n_states = len(distributions)
        # Starts -- prefer neutral
        start_probabilities = [0.95, 0.05]
        # Prefer to keep the current state in each transition
        # All other transitions are equally likely, to start
        transition_matrix = (
            np.identity(n_states) * 100 + np.ones((n_states, n_states)) / n_states
        )
        model = pom.HiddenMarkovModel.from_matrix(
            transition_matrix,
            distributions,
            start_probabilities,
            state_names=state_names,
            name="loh",
        )

        model.fit(
            sequences=[observations],
            edge_inertia=0.1,
            lr_decay=0.75,
            pseudocount=5,
            use_pseudocount=True,
            max_iterations=100000,
            # n_jobs=1,  # processes,
            verbose=False,
        )
        states = np.array(model.predict(observations, algorithm="map"))

        logging.info("Done, now finalizing")
        logging.debug("Model states: %s", model.states)
        logging.debug("Predicted states: %s", states[:100])
        logging.debug(str(collections.Counter(states)))
        # logging.debug("Observations: %s", observations[0][:100])
        logging.debug("Edges: %s", model.edges)

        # Merge adjacent bins with the same state to create segments
        fake_cnarr = CNA(varr.add_columns(weight=1, log2=0, gene=".").data)
        results = squash_by_groups(fake_cnarr, varr.as_series(states), by_arm=True)
        assert (results.start < results.end).all()

    else:
        results = None

    if results is not None and len(results) > 1:
        logging.info(
            "Segment %s:%d-%d on allele freqs for %d additional breakpoints",
            segment.chromosome,
            segment.start,
            segment.end,
            len(results) - 1,
        )
        # Place breakpoints midway between SNVs
        # XXX TODO use original cnarr bin boundaries to select/adjust breakpoint
        mid_breakpoints = (results.start.values[1:] + results.end.values[:-1]) // 2
        starts = np.concatenate([[segment.start], mid_breakpoints])
        ends = np.concatenate([mid_breakpoints, [segment.end]])
        dframe = pd.DataFrame(
            {
                "chromosome": segment.chromosome,
                "start": starts,
                "end": ends,
                # 'baf': results['mean'],
                "gene": segment.gene,  # '-'
                "log2": segment.log2,
                "probes": results["probes"],
                # 'weight': (segment.weight * results['probes']
                #            / (segment.end - segment.start)),
            }
        )
        bad_segs_idx = dframe.start >= dframe.end
        if bad_segs_idx.any():
            raise RuntimeError(
                f"Improper post-processing of segment {segment} -- "
                f"{bad_segs_idx.sum()} bins start >= end:\n{dframe[bad_segs_idx]}\n"
            )

    else:
        dframe = pd.DataFrame(
            {
                "chromosome": segment.chromosome,
                "start": segment.start,
                "end": segment.end,
                "gene": segment.gene,  # '-',
                "log2": segment.log2,
                "probes": segment.probes,
                # 'weight': segment.weight,
            },
            index=[0],
        )

    return dframe

def observations_matrix(cnarray, snps_cpd_means):
    obs = []
    values = cnarray.as_dataframe(cnarray.data)
    values.data.reset_index(drop=True, inplace=True)
    if snps_cpd_means:
        half_values = values.data[(values.data['chromosome'] == "chr1") & (values.data['start'] == 0)].index[1]
    else:
        half_values = values.data[(values.data['start'] == 0)].index[1]

    start = 0
    for i in range(2): #len(values) // 2):
        cnarr = values.data.iloc[start:half_values, ]
        start = half_values
        half_values += (len(values) - half_values)

        meta = {"sample_id": 'sample'}
        cnarr = CNA(cnarr, meta)
        obs.extend(as_observation_matrix(cnarr.autosomes()))

    return obs