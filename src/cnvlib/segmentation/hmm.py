"""Segmentation by Hidden Markov Model."""
import collections
import logging

import numpy as np
import pandas as pd
import pomegranate as pom
import scipy.special

from ..cnary import CopyNumArray as CNA
from ..descriptives import biweight_midvariance
from ..segfilters import squash_by_groups
from ..cluster import kmeans_clustering



from pomegranate import State, NormalDistribution, IndependentComponentsDistribution, DiscreteDistribution
from pomegranate import HiddenMarkovModel as HMM

def segment_hmm(depth, arguments, cnarr, method, window=None, variants=None, processes=1):
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

    logging.info("Building model from observations")
    model = hmm_get_model(depth, arguments, cnarr, method, processes)

    logging.info("Predicting states from model")
    observations = observations_matrix(cnarr)#as_observation_matrix(cnarr)
    states = np.concatenate(
        [np.array(model.predict(obs, algorithm="map")) for obs in observations]
    )

    logging.info("Done, now finalizing")
    logging.debug("Model states: %s", model.states)
    logging.debug("Predicted states: %s", states[:100])
    logging.debug(str(collections.Counter(states)))
    logging.debug("Observations: %s", observations[0][:100])
    logging.debug("Edges: %s", model.edges)

    values = cnarr.as_dataframe(cnarr.data)
    values.data.reset_index(drop=True, inplace=True)
    half_values = values.data[(values.data['chromosome'] == "chr1") & (values.data['start'] == 0)].index[1]

    segs = []
    start = 0
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

        mean = []
        df = cnarray.data.assign(group=state)
        for i, row in segarr.data.iterrows():
            size = len(segarr.data.index)
            if row.chromosome == 'chr15':
                print("here")
            if (row.end - row.start < 1000000) and (i > 0 and i < size-1):
                if segarr.data.at[i - 1, 'state'] == segarr.data.at[i + 1, 'state'] and segarr.data.at[i, 'chromosome'] == segarr.data.at[i - 1, 'chromosome']:
                    segarr.data.at[i, 'state'] = segarr.data.at[i+1, 'state']
                elif (not segarr.data.at[i - 1, 'state'] == segarr.data.at[i + 1, 'state']) and segarr.data.at[i, 'chromosome'] == segarr.data.at[i - 1, 'chromosome']:
                        segarr.data.at[i, 'state'] = segarr.data.at[i - 1, 'state']

        for i in range(len(np.unique(df['group']))):
            mean.append(np.mean(df[df['group'] == i]['depth']))
            segarr['state'] = np.where(segarr['state'] == i, mean[i], segarr['state'])
        segs.append(segarr)

    return segs

def squash_regions(df):

    return pd.DataFrame(df)
def hmm_get_model(depth_values, arguments, cnarr, method, processes):
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
    assert method in ("hmm-tumor", "hmm-germline", "hmm")
    observations = observations_matrix(cnarr)#as_observation_matrix(cnarr.autosomes())

    # Estimate standard deviation from the full distribution, robustly
    stdev = biweight_midvariance(np.concatenate(observations), initial=0)

    u_labels, labels, centers, stdev, clusters = kmeans_clustering(depth_values, arguments['no_of_clusters'])

    state_names = []#["copy_1", "copy_2", "copy_3", "copy_4", "copy_5"]
    distributions = []
    for i in range(len(u_labels)):
        state_names.append("copy_"+str(i))
        distributions.append(pom.NormalDistribution(centers[i], stdev[i], frozen=False))

    # distributions = [
    #     pom.NormalDistribution(centers[0], stdev[0], frozen=False),
    #     pom.NormalDistribution(centers[1], stdev[1], frozen=False),
    #     pom.NormalDistribution(centers[2], stdev[2], frozen=False),
    #     pom.NormalDistribution(centers[3], stdev[3], frozen=False),
    #     pom.NormalDistribution(centers[4], stdev[4], frozen=False),
    # ]

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
        name=method,
    )

    model.fit(
        sequences=observations,
        weights=[len(obs) for obs in observations],
        distribution_inertia=0.8,  # Allow updating dists, but slowly
        edge_inertia=0.1,
        # lr_decay=.75,
        pseudocount=5,
        use_pseudocount=True,
        max_iterations=100000,
        n_jobs=processes,
        verbose=False,
    )

    return model


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

def observations_matrix(cnarray):
    obs = []
    values = cnarray.as_dataframe(cnarray.data)
    values.data.reset_index(drop=True, inplace=True)

    half_values = values.data[(values.data['chromosome'] == "chr1") & (values.data['start'] == 0)].index[1]
    start = 0
    for i in range(2): #len(values) // 2):
        cnarr = values.data.iloc[start:half_values, ]
        start = half_values
        half_values += (len(values) - half_values)

        meta = {"sample_id": 'sample'}
        cnarr = CNA(cnarr, meta)
        obs.extend(as_observation_matrix(cnarr.autosomes()))

    return obs