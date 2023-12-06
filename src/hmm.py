import collections
import logging
import statistics

import numpy as np
import pandas as pd
import pomegranate as pom

from sklearn.metrics import silhouette_score
import scipy

from cnvlib.cnary import CopyNumArray as CNA
from cnvlib.segfilters import squash_by_groups

from clustering import cluster_coverage_data, change_point_detection_means
def hmm_model_creation(data, observations, means, stdev):
    state_names = []
    distributions = []
    for i in range(len(means)):
        state_names.append("copy_" + str(i))
        distributions.append(pom.NormalDistribution(means[i], stdev[i], frozen=False))

    n_states = len(distributions)
    # Starts -- prefer neutral
    binom_coefs = scipy.special.binom(n_states - 1, range(n_states))
    start_probabilities = binom_coefs / binom_coefs.sum()

    transition_matrix = (
            np.identity(n_states) * 100 + np.ones((n_states, n_states)) / n_states
    )
    model = pom.HiddenMarkovModel.from_matrix(
        transition_matrix,
        distributions,
        start_probabilities,
        state_names=state_names,
        name='hmm',
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
        n_jobs=1,
        verbose=False,
    )
    return model

def call_copynumbers(arguments, cnarr, df_chrom, cpd_means_collective):
    observations = observations_matrix(cnarr, cpd_means_collective)

    if cpd_means_collective:
        means, stdev = cluster_coverage_data(cpd_means_collective)
        individual_cpd_means = None
        print(means, stdev)
        #means = [0.0, 9.0, 22.105, 33.3, 46.0]
        #stdev = [2, 10, 20, 20, 10]
        #means = [7, 26.105, 46.3, 68]
        #stdev = [7, 10, 15, 15]
    else:
        individual_cpd_means = change_point_detection_means(arguments, df_chrom)
        means, stdev = cluster_coverage_data(individual_cpd_means)
        print(means, stdev)
    model = hmm_model_creation(cnarr, observations, means, stdev)
    logging.info("Predicting states from model")
    states = np.concatenate(
        [np.array(model.predict(obs, algorithm="map")) for obs in observations]
    )
    values = cnarr.as_dataframe(cnarr.data)
    values.data.reset_index(drop=True, inplace=True)
    if cpd_means_collective:
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

        # for i, row in segarr.data.iterrows():
        #     size = len(segarr.data.index)
        #     #if row.chromosome == 'chr15':
        #     #    print("here")
        #     if (row.end - row.start < 1000000) and (i > 0 and i < size-1):
        #         if segarr.data.at[i - 1, 'state'] == segarr.data.at[i + 1, 'state'] and segarr.data.at[i, 'chromosome'] == segarr.data.at[i - 1, 'chromosome']:
        #             segarr.data.at[i, 'state'] = segarr.data.at[i+1, 'state']
        #         elif (not segarr.data.at[i - 1, 'state'] == segarr.data.at[i + 1, 'state']) and segarr.data.at[i, 'chromosome'] == segarr.data.at[i - 1, 'chromosome']:
        #                 segarr.data.at[i, 'state'] = segarr.data.at[i - 1, 'state']

        for i, val in enumerate(np.unique(df['group'])):
            #mean.append(np.mean(df[df['group'] == i]['depth']))
            #segarr['state'] = np.where(segarr['state'] == i, mean[i], segarr['state'])
            segarr['state'] = np.where(segarr['state'] == val, means[val], segarr['state'])  #TODO Bug, when same center and state values
        segs.append(segarr)

    header = ['chromosome', 'start', 'end', 'state', 'haplotype']
    haplotype_dfs.to_csv('data/copynumbers_segments.csv', sep='\t', columns=header, index=False)

    return segs, states, means, stdev, individual_cpd_means

def silhouette_score_from_states(X, states):
    y = np.concatenate(list(X)).ravel().tolist()
    lst = [i for i in range(0, (len(y) // 2) * 50000, 50000)]
    x = lst + lst
    stdev = []
    means = []
    for g in np.unique(states):
        ix = [index for index, i in enumerate(x) if states[index] == g]
        xn = [x[i] for i in ix]
        yn = [y[i] for i in ix]
        if len(yn) > 1:
            stdev.append(statistics.stdev(yn))
            means.append(statistics.mean(yn))
    return silhouette_score(X, states, metric="euclidean")

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

def as_observation_matrix(cnarr):
    observations = [arm.log2.values for _c, arm in cnarr.by_arm()]
    return observations