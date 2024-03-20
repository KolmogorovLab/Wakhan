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

from clustering import cluster_coverage_data, stdev_calc
from extras import get_contigs_list
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

def call_copynumbers(arguments, cnarr, df_chrom, cpd_means, cpd_means_collective):
    observations = observations_matrix(cnarr, cpd_means_collective, arguments)

    if cpd_means_collective:
        print(cpd_means_collective)
        means, stdev = cluster_coverage_data(cpd_means_collective, True)
        individual_cpd_means = None
        print(cpd_means_collective)
        print(means, stdev)

        #means = [0.0, 2.0, 6.583333333333333, 12.0, 15.75, 19.0, 23.428571428571427, 31.0]
        #stdev = stdev_calc(means)

        # depth_values_hp1 = np.array(cnarr.data['depth'], dtype='int')
        # depth_values_hp2 = np.array(cnarr.data['depth'], dtype='int')
        # depth_values_hp1 = depth_values_hp1.reshape(-1, 1)
        # depth_values_hp2 = depth_values_hp2.reshape(-1, 1)
        # X = np.concatenate([depth_values_hp1])
        # means_, stdev_, score_ = hmm_pome_test(cnarr, X, means, stdev, cpd_means_collective, arguments)
        # print(score_, means_)
        # for i in range(len(means)):
        #     means_, stdev_, score_ = hmm_pome_test(cnarr, X, means[i], stdev[i], cpd_means_collective, arguments)
        #     print(score_, means_)
        #Debug
        #means = [0, 8, 19, 35, 46, 58, 120]
        #stdev = [1, 5, 5, 5, 5, 5, 5]

        #C21
        #means = [0, 11, 22, 26.5, 32, 37, 41]
        #stdev = [1, 5, 5, 5, 5, 5, 5]

        #C18
        #means = [0, 6, 12, 18, 25]
        #stdev = [1, 5, 5, 5, 5]

        means = [mean + 0.025 for mean in means]

        #df_means_chr = df_means_chr_all
    else:
        print(cpd_means)
        #individual_cpd_means = [0.0,0.2,0.3, 33,34,35, 66,67,68, 99,98,97]
        means, stdev = cluster_coverage_data(cpd_means, False)
        print(means, stdev)

    model = hmm_model_creation(cnarr, observations, means, stdev)
    logging.info("Predicting states from model")
    states = np.concatenate(
        [np.array(model.predict(obs, algorithm="map")) for obs in observations]
    )
    values = cnarr.as_dataframe(cnarr.data)
    values.data.reset_index(drop=True, inplace=True)
    if cpd_means_collective:
        half_values = values.data[(values.data['chromosome'] == arguments['contigs'].split('-')[0]) & (values.data['start'] == 0)].index[1]
    else:
        half_values = values.data[(values.data['start'] == 0)].index[1]

    segs = []
    start = 0

    write_csv = True
    for k in range(2): #range(0, len(values), 49592):#len(values) // 2):
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

        mean = []
        df = cnarray.data.assign(group=state)

        if cpd_means_collective:
            for i, row in segarr.data.iterrows():
                size = len(segarr.data.index)
                if (row.end - row.start < 100000) and (i > 0 and i < size-1):
                    if segarr.data.at[i - 1, 'state'] == segarr.data.at[i + 1, 'state'] and segarr.data.at[i, 'chromosome'] == segarr.data.at[i - 1, 'chromosome']:
                        segarr.data.at[i, 'state'] = segarr.data.at[i+1, 'state']
                    elif (not segarr.data.at[i - 1, 'state'] == segarr.data.at[i + 1, 'state']) and segarr.data.at[i, 'chromosome'] == segarr.data.at[i - 1, 'chromosome']:
                            segarr.data.at[i, 'state'] = segarr.data.at[i - 1, 'state']

        for i, val in enumerate(np.unique(df['group'])):
            segarr['depth'] = np.where(segarr['state'] == val, i, segarr['state'])  #TODO Bug, when same center and state values

        for i, val in enumerate(np.unique(df['group'])):
            #mean.append(np.mean(df[df['group'] == i]['depth']))
            #segarr['state'] = np.where(segarr['state'] == i, mean[i], segarr['state'])
            segarr['state'] = np.where(segarr['state'] == val, means[val], segarr['state'])  #TODO Bug, when same center and state values

        # if cpd_means_collective:
        #     steps = np.array(means)
        #     df_means_chr_all["state"] = df_means_chr_all["mean"].apply(lambda x: steps[np.argmin(np.abs(x - steps))])

        if cpd_means_collective:
            segarr = merge_adjacent_regions(segarr, arguments)
            #TODO: merge_adjacent_regions prior?
            # haplotype_df = segarr
            # if not arguments['without_phasing']:
            #     haplotype_df['haplotype'] = k + 1
            #     header = ['chromosome', 'start', 'end', 'depth', 'state', 'haplotype']
            # else:
            #     header = ['chromosome', 'start', 'end', 'depth', 'state']
            # if write_csv:
            #     haplotype_df.to_csv('data/' + arguments['genome_name'] + '_copynumbers_segments.csv', sep='\t', columns=header,
            #                      index=False, mode='a', header=False)
            #     if arguments['without_phasing']:
            #         write_csv = False
        else:
            segarr = merge_adjacent_regions_1(segarr)

        segs.append(segarr)

    return segs, states, means, stdev

def merge_adjacent_regions(segarr, arguments):
    chroms = get_contigs_list(arguments['contigs'])
    dfs = []
    for index, chrom in enumerate(chroms):
        seg = segarr.data[segarr.data['chromosome'] == chrom]
        label_groups = seg['state'].ne(seg['state'].shift()).cumsum()
        df = (seg.groupby(label_groups).agg({'chromosome': 'first', 'start': 'min', 'end': 'max', 'depth': 'first', 'state': 'first'}).reset_index(drop=True))
        dfs.append(df)
    out = pd.concat(dfs)
    return out

def merge_adjacent_regions_1(segarr):
    #label_groups = segarr.data['state'].ne(segarr.data['state'].shift()).cumsum()
    #out = (segarr.data.groupby(label_groups).agg({'chromosome': 'first', 'start': 'min', 'end': 'max', 'depth': 'first', 'state': 'first'}).reset_index(drop=True))
    out = segarr.data
    return out
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

def observations_matrix(cnarray, snps_cpd_means, arguments):
    obs = []
    values = cnarray.as_dataframe(cnarray.data)
    values.data.reset_index(drop=True, inplace=True)
    if snps_cpd_means:
        half_values = values.data[(values.data['chromosome'] == arguments['contigs'].split('-')[0]) & (values.data['start'] == 0)].index[1]
    else:
        half_values = values.data[(values.data['start'] == 0)].index[1]

    start = 0
    for i in range(2): #len(values) // 2):
        cnarr = values.data.iloc[start:half_values, ]
        start = half_values
        half_values += (len(values) - half_values)

        meta = {"sample_id": 'sample'}
        cnarr = CNA(cnarr, meta)
        obs.extend(as_observation_matrix(cnarr))

    return obs

def as_observation_matrix(cnarr):
    observations = [arm.log2.values for _c, arm in cnarr.by_arm()]
    return observations

def hmm_pome_test(cnarr, X, means, stdev, snps_cpd_means_input, arguments):
    import scipy
    import pomegranate as pom
    import scipy.special

    observations = observations_matrix(cnarr, snps_cpd_means_input, arguments)

    state_names = []#["copy_1", "copy_2", "copy_3", "copy_4", "copy_5"]
    distributions = []
    for i in range(len(means)):
        state_names.append("copy_"+str(i))
        distributions.append(pom.NormalDistribution(means[i], stdev[i], frozen=False))

    n_states = len(distributions)
    # Starts -- prefer neutral
    binom_coefs = scipy.special.binom(n_states - 1, range(n_states))
    start_probabilities = binom_coefs / binom_coefs.sum()

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
    states = np.concatenate(
        [np.array(model.predict(obs, algorithm="map")) for obs in observations]
    )

    # y = np.concatenate(list(X)).ravel().tolist()
    # lst = [i for i in range(0, (len(y) // 2) * 50000, 50000)]
    # x = lst + lst
    # stdev = []
    # means = []
    # for g in np.unique(states):
    #     ix = [index for index, i in enumerate(x) if states[index] == g]
    #     xn = [x[i] for i in ix]
    #     yn = [y[i] for i in ix]
    #     if len(yn) > 1:
    #         stdev.append(statistics.stdev(yn))
    #         means.append(statistics.mean(yn))

    score = silhouette_score(X, states, metric="euclidean")

    return means, stdev, score