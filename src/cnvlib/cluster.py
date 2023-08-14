#!/usr/bin/env python
"""Markov-cluster control samples' correlation matrix.

See:

    Stijn van Dongen, Graph Clustering by Flow Simulation,
    PhD thesis, University of Utrecht, May 2000.
    https://micans.org/mcl/

"""
import logging

from sklearn.cluster import KMeans
from kneed import KneeLocator
import plotly.figure_factory as ff
import numpy as np
from hmmlearn import hmm
from sklearn.metrics import silhouette_score
from scipy.special import logsumexp


def kmeans(samples, k=None):
    from scipy.cluster import vq

    if k is None:
        from math import log

        k = max(1, int(round(log(len(samples), 3))))
        # E.g. n=66 -> k=2, 16 -> 3, 47 -> 4, 141 -> 5, 421 -> 6, 1263 -> 7

    print("Clustering", len(samples), "samples by k-means, where k =", k)
    obs = pca_sk(samples, 3)
    obs = vq.whiten(obs)  # Needed?
    _centroids, labels = vq.kmeans2(obs, k, minit="++")
    # XXX shorter way to do this grouping?
    from collections import defaultdict

    clusters = defaultdict(list)
    for idx, label in enumerate(labels):
        clusters[label].append(idx)
    clusters = clusters.values()
    # plot_clusters(obs, clusters)
    return clusters


def markov(samples, inflation=5, max_iterations=100, by_pca=True):
    """Markov-cluster control samples by their read depths' correlation.

    Each of the matrices in the resulting iterable (list) can be processed the
    same as the input to calculate average log2 and spread values for that
    cluster.

    Parameters
    ----------
    samples : array
        Matrix of samples' read depths or normalized log2 values, as columns.
    inflation : float
        Inflation parameter for MCL. Must be >1; higher more granular clusters.
    by_pca : bool
        If true, similarity is by PCA; otherwise, by Pearson correlation.

    Return
    ------
    results : list
        A list of matrices representing non-overlapping column-subsets of the
        input, where each set of samples represents a cluster.
    """
    if inflation <= 1:
        raise ValueError("inflation must be > 1")

    if by_pca:
        pca_matrix = pca_sk(samples, 2)  # pca_plain
        # Convert to similarity matrix
        from scipy.spatial import distance

        dists = distance.squareform(distance.pdist(pca_matrix))
        M = 1 - (dists / dists.max())
    else:
        M = np.corrcoef(samples)

    M, clusters = mcl(M, max_iterations, inflation)
    # plot_clusters(M, clusters)
    return clusters


# https://github.com/koteth/python_mcl/blob/master/mcl/mcl_clustering.py
# https://github.com/GuyAllard/markov_clustering/blob/master/markov_clustering/mcl.py
# https://stackoverflow.com/questions/44243525/mcl-clustering-implementation-in-python-deal-with-overlap
def mcl(M, max_iterations, inflation, expansion=2):
    """Markov cluster algorithm."""
    print("M_init:\n", M)
    M = normalize(M)
    # print("M_norm:\n", M)
    for i in range(max_iterations):
        M_prev = M
        M = inflate(expand(M, expansion), inflation)
        # print("M_inflate_%d:\n" % i, M)
        if converged(M, M_prev):
            logging.debug("Converged at iteration %d", i)
            break
        M = prune(M)
        # print("M_prune_%d:\n" % i, M)

    # print("M_final:\n", M)
    clusters = get_clusters(M)
    return M, clusters


def normalize(A):
    """Normalize matrix columns."""
    return A / A.sum(axis=0)


def inflate(A, inflation):
    """Apply cluster inflation with the given element-wise exponent.

    From the mcl manual:

    This value is the main handle for affecting cluster granularity.
    This parameter is the usually only one that may require tuning.

    By default it is set to 2.0 and this is a good way to start. If you want to
    explore cluster structure in graphs with MCL, vary this parameter to obtain
    clusterings at different levels of granularity.  It is usually chosen
    somewhere in the range [1.2-5.0]. -I 5.0 will tend to result in fine-grained
    clusterings, and -I 1.2 will tend to result in very coarse grained
    clusterings. A good set of starting values is 1.4, 2, 4, and 6.
    Your mileage will vary depending on the characteristics of your data.

    Low values for -I, like -I 1.2, will use more CPU/RAM resources.

    Use mcl's cluster validation tools 'clm dist' and 'clm info' to test the
    quality and coherency of your clusterings.
    """
    return normalize(np.power(A, inflation))


def expand(A, expansion):
    """Apply cluster expansion with the given matrix power."""
    return np.linalg.matrix_power(A, expansion)


def converged(M, M_prev):
    """Test convergence.

    Criterion: homogeneity(??) or no change from previous round.
    """
    return np.allclose(M, M_prev)


# https://stackoverflow.com/questions/17772506/markov-clustering
def get_clusters(M):
    """Extract clusters from the matrix.

    Interpretation: "Attractors" are the non-zero elements of the matrix
    diagonal. The nodes in the same row as each attractor form a cluster.

    Overlapping clusterings produced by MCL are extremely rare, and always a
    result of symmetry in the input graph.

    Returns
    -------
    result : list
        A list of arrays of sample indices. The indices in each list item
        indicate the elements of that cluster; the length of the list is the
        number of clusters.
    """
    attractors_idx = M.diagonal().nonzero()[0]
    clusters_idx = [M[idx].nonzero()[0] for idx in attractors_idx]
    return clusters_idx


# https://github.com/GuyAllard/markov_clustering/blob/master/markov_clustering/mcl.py
def prune(M, threshold=0.001):
    """Remove many small entries while retaining most of M's stochastic mass.

    After pruning, vectors are rescaled to be stochastic again.
    (stochastic: values are all non-negative and sum to 1.)

    This step is purely to keep computation tractable in mcl by making the
    matrix more sparse (i.e. full of zeros), enabling sparse-matrix tricks to
    work.

    ----

    mcl:
        The default setting is something like -P 4000 -S 500 -R 600, where:

      -P <int> (1/cutoff)
      -S <int> (selection number)
      -R <int> (recover number)
      ---
      -pct <pct> (recover percentage)
      -p <num> (cutoff)

    After computing a new (column stochastic) matrix vector during expansion
    (which  is  matrix  multiplication c.q.  squaring), the vector is
    successively exposed to different pruning strategies. Pruning effectively
    perturbs the MCL process a little in order to obtain matrices that are
    genuinely sparse, thus keeping the computation tractable.

    mcl proceeds as follows:

    First, entries that are smaller than cutoff are
    removed, resulting in a vector with  at most 1/cutoff entries.

        * The cutoff can be supplied either by -p, or as the inverse value by
        -P.  The latter is more intuitive, if your intuition is like mine (P
        stands for precision or pruning).

    Second, if the remaining stochastic mass (i.e. the sum of all remaining
    entries) is less than <pct>/100 and the number of remaining entries is
    less than <r> (as specified by the -R flag), mcl will try to regain ground
    by recovering the largest discarded entries. If recovery was not necessary,
    mcl tries to prune the vector further down to at most s entries (if
    applicable), as specified by the -S flag. If this results in a vector that
    satisfies the recovery condition then recovery is attempted, exactly as
    described above. The latter will not occur of course if <r> <= <s>.

    """
    pruned = M.copy()
    pruned[pruned < threshold] = 0
    return pruned


def pca_sk(data, n_components=None):
    """Principal component analysis using scikit-learn.

    Parameters
    ----------
    data : 2D NumPy array
    n_components : int

    Returns: PCA-transformed data with `n_components` columns.
    """
    from sklearn.decomposition import PCA

    return PCA(n_components=n_components).fit_transform(data)


def pca_plain(data, n_components=None):
    """Principal component analysis using numpy eigenvalues.

    Source:
    https://stackoverflow.com/questions/13224362/principal-component-analysis-pca-in-python

    Parameters
    ----------
    data : 2D NumPy array
    n_components : int

    Returns: PCA-transformed data with `n_components` columns.
    """
    # Standarize the data
    data = data.copy()
    data -= data.mean(axis=0)
    data /= data.std(axis=0)
    # Covariance matrix
    C = np.cov(data)
    # Eigenvectors & eigenvalues of covariance matrix
    # ('eigh' rather than 'eig' since C is symmetric, for performance)
    E, V = np.linalg.eigh(C)
    # Sort eigenvalues in decreasing order, & eigenvectors to match,
    # keeping only the first `n_components` components
    key = np.argsort(E)[::-1][:n_components]
    E, V = E[key], V[:, key]
    # Transformation the data using eigenvectors
    U = np.dot(data, V)  # or: np.dot(V.T, data.T).T
    # U = np.dot(V.T, data.T).T
    # Return the re-scaled data, eigenvalues, and eigenvectors
    return U  # , E, V


def plot_clusters(M, cluster_indices):
    """Scatter plot first 2 components, colorized by cluster.

    Parameters
    ----------
    M : np.array
        PCA'd matrix. Rows are samples.
    cluster_indices : iterable of np.array
        Indices of samples in each cluster, as present in M.
    """
    from matplotlib import pyplot as plt

    # colors = list("krgbo")[:len(cluster_indices)]
    _fig, ax = plt.subplots(1, 1)
    for cl_idx in cluster_indices:
        ax.scatter(M[cl_idx, 0], M[cl_idx, 1])  # c=color
    plt.show()

######################################################################
#
#
#
######################################################################

def kmeans_clustering_seed(depth_values):
    sum_of_squared_distances = []
    K = range(1, 15)
    for k in K:
        km = KMeans(n_clusters=k)
        km = km.fit(depth_values)
        sum_of_squared_distances.append(km.inertia_)

    kn = optimal_clusters_size_knee(K, sum_of_squared_distances)
    return kn.knee
def kmeans_clustering(depth_values, no_of_clusters=None):
    depth_values = np.clip(depth_values, a_min=1, a_max=150)
    depth_values[depth_values != 0]

    hp = list(depth_values)
    depth_values = depth_values.reshape(-1, 1)

    if no_of_clusters == None:
        no_of_clusters = kmeans_clustering_seed(depth_values)

    km = KMeans(n_clusters=no_of_clusters, n_init=25, max_iter=600, random_state=0)
    km = km.fit(depth_values, len(depth_values))
    labels = list(km.labels_)
    u_labels = np.unique(labels)
    centers = list(np.concatenate(km.cluster_centers_))

    stdev = []
    clusters = []
    for i in range(len(u_labels)):
        clusters.append([int(a) for a, b in zip(hp, labels) if b == i])
        stdev.append(np.std([(a, b) for a, b in zip(hp, labels) if b == i]))

    return u_labels, labels, centers, stdev, clusters

def optimal_clusters_size_knee(K, sum_of_squared_distances):
    return KneeLocator(x=K, y=sum_of_squared_distances, curve='convex', direction='decreasing')

def optimal_clusters_size_silhouette_score(depth_values, km):
    score = silhouette_score(depth_values, km.labels_, sample_size=200)
def plot_optimal_clusters(hp1, hp2, unphased, arguments, chrom=None):

    hp1=np.clip(hp1, a_min=1, a_max=300)
    hp2=np.clip(hp2, a_min=1, a_max=300)

    hp1[hp1 != 0]
    hp1=hp1.reshape(-1, 1)

    hp2[hp2 != 0]
    hp2=hp2.reshape(-1, 1)
    depth_values = np.concatenate((hp1, hp2))

    depth_values = [int(x) for x in depth_values]

    u_labels, labels, centers, stdev, clusters = kmeans_clustering(depth_values, arguments['no_of_clusters'])
    plot_clsters_settings(clusters, arguments, chrom)

    if arguments['no_of_clusters'] == None:
        #nearest_clusters = []
        for no_of_clusters in range(len(clusters)-1, len(clusters)+2):
            u_labels, labels, centers, stdev, clusters = kmeans_clustering(depth_values, no_of_clusters)
            plot_clsters_settings(clusters, arguments, chrom)


def plot_clusters(plt, X, y=None):
    plt.scatter(X[:, 0], X[:, 1], c=y, s=500, cmap='autumn')
    plt.xlabel("$x_1$", fontsize=14)
    plt.ylabel("$x_2$", fontsize=14, rotation=0)



def plot_clsters_settings(clusters, arguments, chrom):
    group_labels = ['Cluster '+str(i) for i in range(len(clusters))]
    fig = ff.create_distplot(clusters, group_labels, curve_type='normal', # override default 'kde'
    )

    fig.update_layout(
        title={
            'text': 'Genome - ' + arguments['genome_name'],
            'y': 0.96,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'},

        font_family="Courier New",
        font_color="dimgray",
        title_font_family="Times New Roman",
        title_font_color="red",
        legend_title_font_color="green",
    )
    if not chrom == None:
        fig.write_html(arguments['out_dir_plots'] +'/'+ arguments['genome_name'] + "_genome_bins_clusters["+str(len(clusters))+"]_histogram_"+chrom+".html")
    else:
        fig.write_html(arguments['out_dir_plots'] +'/'+ arguments['genome_name'] + "_genome_bins_clusters["+str(len(clusters))+"]_histogram.html")

def bayesian_gaussian_mixture_clustering(X, length, clouds=None, concentration_prior=None, K=10, restarts=5, seed=0):

    from sklearn.mixture import BayesianGaussianMixture
    from collections import Counter

    #total = list(points)
    # if clouds is not None:
    #     total.extend(list(clouds))
    # npArray = np.array(total)

    # bgm = BayesianGaussianMixture(n_components=10, n_init=10, random_state=42)
    # bgm.fit(npArray)
    # np.round(bgm.weights_, 2)

    gmm = BayesianGaussianMixture(
        n_components=K,
        n_init=restarts,
        weight_concentration_prior=concentration_prior,
        max_iter=int(1e6),
        random_state=seed,
    )
    targetAssignments = gmm.fit_predict(X, length)
    #targetAssignments = targetAssignments[: length]
    mus = gmm.means_
    sigmas = gmm.covariances_
    cntr = Counter(targetAssignments)
    numPoints = [cntr[i] if i in cntr else 0 for i in range(K)]
    numClusters = len(cntr)

    return mus, sigmas, targetAssignments, numPoints, numClusters

def gaussian_hmm_clustering(X, length):
    from hmmlearn import hmm
    model = hmm.GaussianHMM(n_components=10, n_iter=2_000).fit(X, length)

    mus = np.ravel(model.means_)
    sigmas = np.ravel(np.sqrt([np.diag(c) for c in model.covars_]))
    P = model.transmat_

    return mus, sigmas, P

#def optimal_clusters_size_silhouette_score():
    #score = silhouette_score(hp1, km.labels_, sample_size=200)

    #mus, sigmas, targetAssignments, numPoints, numClusters = clusters_gmm(hp1, )
    #mus, sigmas, targetAssignments, numPoints, numClusters = clusters_gmm(hp2, )

    #For generating some data
    # from sklearn.datasets import make_blobs
    #
    # from matplotlib import pyplot as plt
    # from sklearn.mixture import BayesianGaussianMixture
    #
    # X, y = make_blobs(n_samples=350, centers=10, cluster_std=0.60)
    # plt.scatter(X[:, 0], X[:, 1], cmap='viridis')
    # bay_gmm = BayesianGaussianMixture(n_components=10, n_init=10)
    #
    # bay_gmm.fit(X)
    # bay_gmm_weights = bay_gmm.weights_
    # np.round(bay_gmm_weights, 2)
    # n_clusters_ = (np.round(bay_gmm_weights, 2) > 0).sum()
    # print('Estimated number of clusters: ' + str(n_clusters_))
    # y_pred = bay_gmm.predict(X)
    # plt.scatter(X[:, 0], X[:, 1], c=y_pred, cmap="viridis")
    # plt.xlabel("Feature 1")
    # plt.ylabel("Feature 2")
    # props = bay_gmm.predict_proba(X)
    # props = props.round(3)
    # props
    # size = 50 * props.max(1) ** 2  # square emphasizes differences
    # plt.scatter(X[:, 0], X[:, 1], c=y_pred, cmap='viridis', s=size)
    # plt.show()
    #
    # fig = px.histogram(df, x="coverage", color="Haplotypes", barmode="overlay", marginal="violin")#, log_y=True)
    #
    # fig.update_yaxes(range=[1, 1000])
    # fig.update_xaxes(range=[1, 36])

def clusters_tests(depth_values, depth_values_1):

    #################################################
    import time
    depth_values = np.clip(depth_values, a_min=1, a_max=300)
    depth_values[depth_values != 0]
    depth_values = depth_values.reshape(-1, 1)

    depth_values_1 = np.clip(depth_values_1, a_min=1, a_max=300)
    depth_values_1[depth_values_1 != 0]
    depth_values_1 = depth_values_1.reshape(-1, 1)

    X = np.concatenate([depth_values, depth_values_1])
    lengths = [len(depth_values), len(depth_values_1)]

    # start = time.time()
    # mus, sigmas, P = gaussian_hmm_clustering(X, lengths)
    # print("gaussian hmm clustering")
    # print(mus, sigmas, P)
    # end = time.time()
    # print(end - start)

    start = time.time()
    best_score, best_model, best_labels, best_K, rs = hmm_model_select(X, lengths)
    print("ghmm hatchet clustering")
    print(best_score, best_model, best_labels, best_K, rs)
    end = time.time()
    print(end - start)

    # start = time.time()
    # mus, sigmas, targetAssignments, numPoints, clusters_m = bayesian_gaussian_mixture_clustering(X, lengths)
    # print("bayesian gaussian mixture clustering")
    # print(mus, sigmas, targetAssignments, numPoints, clusters_m)
    # end = time.time()
    # print(end - start)
    #################################################

def hmm_model_select(X, lengths, minK=2, maxK=10, tau=10e-6, tmat='diag', decode_alg='viterbi', covar='diag', restarts=10):
    assert tmat in ['fixed', 'diag', 'free']
    assert decode_alg in ['map', 'viterbi']

    # format input
    # tracks = [a for a in tracks if a.shape[0] > 0 and a.shape[1] > 0]
    # if len(tracks) > 1:
    #     X = np.concatenate(tracks, axis=1).T
    #     lengths = [a.shape[1] for a in tracks]
    # else:
    #     X = tracks[0].T
    #     lengths = [tracks[0].shape[1]]

    best_K = 0
    best_score = -1.01   # below minimum silhouette score value
    best_model = None
    best_labels = None

    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import silhouette_score
    from scipy.special import logsumexp
    from scipy.spatial.distance import pdist, squareform
    from hmmlearn import hmm

    #scaler = StandardScaler()
    #X_scaled = scaler.fit_transform(X, lengths)
    #C = squareform(pdist(X_scaled))

    rs = {}
    for K in range(minK, maxK + 1):
        # print(K, datetime.now())

        my_best_ll = -1 * np.inf
        my_best_labels = None
        my_best_model = None
        for s in range(restarts):
            # construct initial transition matrix
            A = make_transmat(1 - tau, K)
            assert np.all(A > 0), (
                'Found 0 or negative elements in transition matrix.'
                'This is likely a numerical precision issue -- try increasing tau.',
                A,
            )
            assert np.allclose(np.sum(A, axis=1), 1), ('Not all rows in transition matrix sum to 1.', A)

            model = hmm.GaussianHMM(
                n_components=K,
                init_params='mc',
                params='smc',
                covariance_type=covar,
                random_state=s,
            )

            model.startprob_ = np.ones(K) / K
            model.transmat_ = A
            model.fit(X, lengths)

            prob, labels = model.decode(X, lengths, algorithm=decode_alg)
            if prob > my_best_ll:
                my_best_labels = labels
                my_best_ll = prob
                my_best_model = model

        score = silhouette_score(X, my_best_labels, metric="euclidean", sample_size=3000)
        print(score, K)

        rs[K] = my_best_ll, score, my_best_labels
        if score > best_score:
            best_score = score
            best_model = my_best_model
            best_labels = my_best_labels
            best_K = K

    return best_score, best_model, best_labels, best_K, rs

def make_transmat(diag, K):
    offdiag = (1 - diag) / (K - 1)
    transmat_ = np.diag([diag - offdiag] * K)
    transmat_ += offdiag
    return transmat_


############################################################################
def hmm_validation():
    import numpy as np
    import matplotlib.pyplot as plt

    from hmmlearn import hmm

    # Prepare parameters for a 4-components HMM
    # Initial population probability
    startprob = np.array([0.6, 0.3, 0.1, 0.0])
    # The transition matrix, note that there are no transitions possible
    # between component 1 and 3
    transmat = np.array([[0.7, 0.2, 0.0, 0.1],
                         [0.3, 0.5, 0.2, 0.0],
                         [0.0, 0.3, 0.5, 0.2],
                         [0.2, 0.0, 0.2, 0.6]])
    # The means of each component
    means = np.array([[0.0, 0.0],
                      [0.0, 11.0],
                      [9.0, 10.0],
                      [11.0, -1.0]])
    # The covariance of each component
    covars = .5 * np.tile(np.identity(2), (4, 1, 1))

    # Build an HMM instance and set parameters
    gen_model = hmm.GaussianHMM(n_components=4, covariance_type="full")

    # Instead of fitting it from the data, we directly set the estimated
    # parameters, the means and covariance of the components
    gen_model.startprob_ = startprob
    gen_model.transmat_ = transmat
    gen_model.means_ = means
    gen_model.covars_ = covars

    # Generate samples
    X, Z = gen_model.sample(500)
    print(X, Z)

    # Plot the sampled data
    fig, ax = plt.subplots()
    ax.plot(X[:, 0], X[:, 1], ".-", label="observations", ms=6,
            mfc="orange", alpha=0.7)

    # Indicate the component numbers
    for i, m in enumerate(means):
        ax.text(m[0], m[1], 'Component %i' % (i + 1),
                size=17, horizontalalignment='center',
                bbox=dict(alpha=.7, facecolor='w'))
    ax.legend(loc='best')
    fig.show()

    ########################################################

    scores = list()
    models = list()
    for n_components in (3, 4, 5):
        for idx in range(10):
            # define our hidden Markov model
            model = hmm.GaussianHMM(n_components=n_components,
                                    covariance_type='full',
                                    random_state=idx)
            model.fit(X[:X.shape[0] // 2])  # 50/50 train/validate
            models.append(model)
            scores.append(model.score(X[X.shape[0] // 2:]))
            print(f'Converged: {model.monitor_.converged}'
                  f'\tScore: {scores[-1]}')

    # get the best model
    model = models[np.argmax(scores)]
    n_states = model.n_components
    print(f'The best model had a score of {max(scores)} and {n_states} '
          'states')

    # use the Viterbi algorithm to predict the most likely sequence of states
    # given the model
    states = model.predict(X)
    ##############################################################

    # plot model states over time
    fig, ax = plt.subplots()
    ax.plot(Z, states)
    ax.set_title('States compared to generated')
    ax.set_xlabel('Generated State')
    ax.set_ylabel('Recovered State')
    fig.show()

    # plot the transition matrix
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5))
    ax1.imshow(gen_model.transmat_, aspect='auto', cmap='spring')
    ax1.set_title('Generated Transition Matrix')
    ax2.imshow(model.transmat_, aspect='auto', cmap='spring')
    ax2.set_title('Recovered Transition Matrix')
    for ax in (ax1, ax2):
        ax.set_xlabel('State To')
        ax.set_ylabel('State From')

    fig.tight_layout()
    fig.show()

##############################################################################

def hmm_validation_depth(hp1, hp2):
    import numpy as np
    import matplotlib.pyplot as plt
    from hmmlearn import hmm

    hp1 = [int(x) for x in hp1]
    hp2 = [int(x) for x in hp2]
    X = np.concatenate([hp1, hp2])
    lengths = [len(hp1), len(hp2)]

    depth = np.clip(X, a_min=1, a_max=150)
    depth[depth != 0]
    depth = depth.reshape(-1, 1)

    import numpy as np
    from scipy.stats import norm

    mu, std = norm.fit(depth)
    print(mu, std)


    n_states = 7
    # Prepare parameters for a 5-components HMM

    # Initial population probability
    startprob = 1. / n_states * np.ones(n_states)

    startprob = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    # The transition matrix
    transmat = np.array([[0.93, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01],
                         [0.02, 0.93, 0.01, 0.01, 0.01, 0.01, 0.01],
                         [0.02, 0.01, 0.93, 0.01, 0.01, 0.01, 0.01],
                         [0.02, 0.01, 0.01, 0.93, 0.01, 0.01, 0.01],
                         [0.02, 0.01, 0.01, 0.01, 0.93, 0.01, 0.01],
                         [0.02, 0.01, 0.01, 0.01, 0.01, 0.93, 0.01],
                         [0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.93]])

    # The means of each component
    means = np.array([[2.44980],
                      [38.5740],
                      [55.5904],
                      [78.9980],
                      [94.5037],
                      [116.50378],
                      [144.987]])

    # The covariance of each component
    covars = np.array([[[1.44980]],
                      [[19.57403]],
                      [[27.59040]],
                      [[38.99801]],
                      [[47.50378]],
                      [[58.50378]],
                      [[72.50378]]])

    # Build an HMM instance and set parameters
    gen_model = hmm.GaussianHMM(n_components=n_states, covariance_type="full")

    gen_model.startprob_ = startprob
    gen_model.transmat_ = transmat
    gen_model.means_ = means
    gen_model.covars_ = covars

    gen_model.fit(depth, lengths)
    states = gen_model.predict(depth)

    # Instead of fitting it from the data, we directly set the estimated
    # parameters, the means and covariance of the components
    # gen_model.startprob_ = startprob
    # gen_model.transmat_ = transmat
    # gen_model.means_ = means
    # gen_model.covars_ = covars

    # Generate samples
    #X, Z = gen_model.sample(500)
    #print(X, Z)

    #import scipy.stats as stats
    #stats.zscore(data) #https://www.kaggle.com/discussions/general/210726

    # Plot the sampled data
    # fig, ax = plt.subplots()
    # ax.plot(X[:, 0], X[:, 1], ".-", label="observations", ms=6,
    #         mfc="orange", alpha=0.7)

    # Indicate the component numbers
    # for i, m in enumerate(means):
    #     ax.text(m[0], m[1], 'Component %i' % (i + 1),
    #             size=17, horizontalalignment='center',
    #             bbox=dict(alpha=.7, facecolor='w'))
    # ax.legend(loc='best')
    # fig.show()

    ########################################################
    model = hmm.GaussianHMM(n_components=n_states, covariance_type="full")
    model.startprob_ = startprob
    model.transmat_ = transmat
    model.means_ = means

    aic = list()
    bic = list()
    scores = list()
    models = list()
    for n_components in (n_states - 1, n_states, n_states + 1):
        for idx in range(25):
            # define our hidden Markov model
            #model = hmm.GaussianHMM(n_components=n_components,
            #                        covariance_type='full')
            model.fit(depth, lengths)
            #Initial guess
            print(model.means_)
            models.append(model)
            scores.append(model.score(depth, lengths))

            aic.append(model.aic(depth, lengths))
            bic.append(model.bic(depth, lengths))

            print(f'Converged: {model.monitor_.converged}'
                  f'\tScore: {scores[-1]}')

    # get the best model
    model = models[np.argmax(scores)]
    n_states = model.n_components
    print(f'The best model had a score of {max(scores)} and {n_states} '
          'states')

    startprob_best = model.startprob_
    transmat_best = model.transmat_
    means_best = model.means_
    covars_best = model.covars_

    model_s = models[np.argmin(aic)]
    print(model_s.means_, model_s.n_components)
    model_s = models[np.argmin(bic)]
    print(model_s.means_, model_s.n_components)
    # use the Viterbi algorithm to predict the most likely sequence of states
    # given the model
    states = model.predict(depth, lengths)

    ##############################################################

    # plot model states over time
    # fig, ax = plt.subplots()
    # ax.plot(Z, states)
    # ax.set_title('States compared to generated')
    # ax.set_xlabel('Generated State')
    # ax.set_ylabel('Recovered State')
    # fig.show()
    #
    # plot the transition matrix
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5))
    ax1.imshow(gen_model.transmat_, aspect='auto', cmap='spring')
    ax1.set_title('Generated Transition Matrix')
    ax2.imshow(model.transmat_, aspect='auto', cmap='spring')
    ax2.set_title('Recovered Transition Matrix')
    for ax in (ax1, ax2):
        ax.set_xlabel('State To')
        ax.set_ylabel('State From')

    fig.tight_layout()
    fig.show()

def poission_hmm():
    import numpy as np
    import matplotlib.pyplot as plt

    from scipy.stats import poisson
    from hmmlearn import hmm

    # earthquake data from http://earthquake.usgs.gov/
    earthquakes = np.array([
        13, 14, 8, 10, 16, 26, 32, 27, 18, 32, 36, 24, 22, 23, 22, 18,
        25, 21, 21, 14, 8, 11, 14, 23, 18, 17, 19, 20, 22, 19, 13, 26,
        13, 14, 22, 24, 21, 22, 26, 21, 23, 24, 27, 41, 31, 27, 35, 26,
        28, 36, 39, 21, 17, 22, 17, 19, 15, 34, 10, 15, 22, 18, 15, 20,
        15, 22, 19, 16, 30, 27, 29, 23, 20, 16, 21, 21, 25, 16, 18, 15,
        18, 14, 10, 15, 8, 15, 6, 11, 8, 7, 18, 16, 13, 12, 13, 20,
        15, 16, 12, 18, 15, 16, 13, 15, 16, 11, 11])

    # Plot the sampled data
    fig, ax = plt.subplots()
    ax.plot(earthquakes, ".-", ms=6, mfc="orange", alpha=0.7)
    ax.set_xticks(range(0, earthquakes.size, 10))
    ax.set_xticklabels(range(1906, 2007, 10))
    ax.set_xlabel('Year')
    ax.set_ylabel('Count')
    fig.show()

    scores = list()
    models = list()
    for n_components in range(1, 4):
        for idx in range(10):  # ten different random starting states
            # define our hidden Markov model
            model = hmm.PoissonHMM(n_components=n_components, random_state=idx,
                                   n_iter=10)
            model.fit(earthquakes[:, None])
            models.append(model)
            scores.append(model.score(earthquakes[:, None]))
            print(f'Converged: {model.monitor_.converged}\t\t'
                  f'Score: {scores[-1]}')

    # get the best model
    model = models[np.argmax(scores)]
    print(f'The best model had a score of {max(scores)} and '
          f'{model.n_components} components')

    # use the Viterbi algorithm to predict the most likely sequence of states
    # given the model
    states = model.predict(earthquakes[:, None])

    # plot model states over time
    fig, ax = plt.subplots()
    ax.plot(model.lambdas_[states], ".-", ms=6, mfc="orange")
    ax.plot(earthquakes)
    ax.set_title('States compared to generated')
    ax.set_xlabel('State')
    plt.show()

def custom_hmmlearn(depth_values):
    import numpy as np
    import matplotlib.pyplot as plt
    from hmmlearn import hmm

    depth_values = np.clip(depth_values, a_min=1, a_max=350)
    #depth_values[depth_values != 0]

    hp = list(depth_values)
    X = depth_values.reshape(-1, 1)

    n_states = 7
    startprob = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    # The transition matrix
    transmat = np.array([[0.93, 0.02, 0.01, 0.01, 0.01, 0.01, 0.01],
                         [0.02, 0.93, 0.01, 0.01, 0.01, 0.01, 0.01],
                         [0.02, 0.01, 0.93, 0.01, 0.01, 0.01, 0.01],
                         [0.02, 0.01, 0.01, 0.93, 0.01, 0.01, 0.01],
                         [0.02, 0.01, 0.01, 0.01, 0.93, 0.01, 0.01],
                         [0.02, 0.01, 0.01, 0.01, 0.01, 0.93, 0.01],
                         [0.02, 0.01, 0.01, 0.01, 0.01, 0.01, 0.93]])

    # The means of each component
    means = np.array([[2.44980],
                      [38.5740],
                      [55.5904],
                      [78.9980],
                      [94.5037],
                      [116.50378],
                      [144.987]])

    # The covariance of each component
    covars = np.array([[[1.44980]],
                      [[19.57403]],
                      [[27.59040]],
                      [[38.99801]],
                      [[47.50378]],
                      [[58.50378]],
                      [[72.50378]]])

    # Build an HMM instance and set parameters
    gen_model = hmm.GaussianHMM(n_components=n_states, covariance_type="full")

    gen_model.startprob_ = startprob
    gen_model.transmat_ = transmat
    gen_model.means_ = means
    gen_model.covars_ = covars

    gen_model.fit(X)
    states = gen_model.predict(X)

    return states


########################################################
def hmm_model_select_hatchet(X, lengths, minK=20, maxK=50, tau=10e-6, tmat='diag', decode_alg='viterbi', covar='diag', restarts=10):
    assert tmat in ['fixed', 'diag', 'free']
    assert decode_alg in ['map', 'viterbi']

    # format input
    # tracks = [a for a in tracks if a.shape[0] > 0 and a.shape[1] > 0]
    # if len(tracks) > 1:
    #     X = np.concatenate(tracks, axis=1).T
    #     lengths = [a.shape[1] for a in tracks]
    # else:
    #     X = tracks[0].T
    #     lengths = [tracks[0].shape[1]]

    best_K = 0
    best_score = -1.01   # below minimum silhouette score value
    best_model = None
    best_labels = None

    # scaler = StandardScaler()
    # X_scaled = scaler.fit_transform(X)
    # C = squareform(pdist(X_scaled))

    rs = {}
    for K in range(minK, maxK + 1):
        # print(K, datetime.now())

        my_best_ll = -1 * np.inf
        my_best_labels = None
        my_best_model = None
        for s in range(restarts):
            # construct initial transition matrix
            A = make_transmat(1 - tau, K)
            assert np.all(A > 0), (
                'Found 0 or negative elements in transition matrix.'
                'This is likely a numerical precision issue -- try increasing tau.',
                A,
            )
            assert np.allclose(np.sum(A, axis=1), 1), ('Not all rows in transition matrix sum to 1.', A)

            if tmat == 'fixed':
                model = hmm.GaussianHMM(
                    n_components=K,
                    init_params='mc',
                    params='smc',
                    covariance_type=covar,
                    random_state=s,
                )
            elif tmat == 'free':
                model = hmm.GaussianHMM(
                    n_components=K,
                    init_params='mc',
                    params='smct',
                    covariance_type=covar,
                    random_state=s,
                )
            else:
                model = DiagGHMM(
                    n_components=K,
                    init_params='mc',
                    params='smct',
                    covariance_type=covar,
                    random_state=s,
                )

            model.startprob_ = np.ones(K) / K
            model.transmat_ = A
            model.fit(X, lengths)

            print(K, s)

            prob, labels = model.decode(X, lengths, algorithm=decode_alg)
            if prob > my_best_ll:
                my_best_labels = labels
                my_best_ll = prob
                my_best_model = model

        score = silhouette_score(X, my_best_labels, metric="euclidean", sample_size=3000)

        rs[K] = my_best_ll, score, my_best_labels
        if score > best_score:
            best_score = score
            best_model = my_best_model
            best_labels = my_best_labels
            best_K = K

    print(best_model.means_, best_model.n_components)

    import matplotlib.pyplot as plt
    # plot the transition matrix
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5))
    ax1.imshow(A, aspect='auto', cmap='spring')
    ax1.set_title('Generated Transition Matrix')
    ax2.imshow(best_model.transmat_, aspect='auto', cmap='spring')
    ax2.set_title('Recovered Transition Matrix')
    for ax in (ax1, ax2):
        ax.set_xlabel('State To')
        ax.set_ylabel('State From')

    fig.tight_layout()
    plt.savefig("tmat.pdf", format="pdf", bbox_inches="tight")

    return best_model.n_components, best_model.means_, best_model.covars_


class DiagGHMM(hmm.GaussianHMM):
    def _accumulate_sufficient_statistics(self, stats, obs, framelogprob, posteriors, fwdlattice, bwdlattice):
        super()._accumulate_sufficient_statistics(stats, obs, framelogprob, posteriors, fwdlattice, bwdlattice)

        if 't' in self.params:
            # for each ij, recover sum_t xi_ij from the inferred transition matrix
            bothlattice = fwdlattice + bwdlattice
            loggamma = (bothlattice.T - logsumexp(bothlattice, axis=1)).T

            # denominator for each ij is the sum of gammas over i
            denoms = np.sum(np.exp(loggamma), axis=0)
            # transpose to perform row-wise multiplication
            stats['denoms'] = denoms

    def _do_mstep(self, stats):
        super()._do_mstep(stats)
        if 't' in self.params:

            denoms = stats['denoms']
            x = (self.transmat_.T * denoms).T

            # numerator is the sum of ii elements
            num = np.sum(np.diag(x))
            # denominator is the sum of all elements
            denom = np.sum(x)

            # (this is the same as sum_i gamma_i)
            # assert np.isclose(denom, np.sum(denoms))

            stats['diag'] = num / denom
            # print(num.shape)
            # print(denom.shape)

            self.transmat_ = self.form_transition_matrix(stats['diag'])

    def form_transition_matrix(self, diag):
        tol = 1e-10
        diag = np.clip(diag, tol, 1 - tol)

        offdiag = (1 - diag) / (self.n_components - 1)
        transmat_ = np.diag([diag - offdiag] * self.n_components)
        transmat_ += offdiag
        # assert np.all(transmat_ > 0), (diag, offdiag, transmat_)
        return transmat_


def make_transmat(diag, K):
    offdiag = (1 - diag) / (K - 1)
    transmat_ = np.diag([diag - offdiag] * K)
    transmat_ += offdiag
    return transmat_
########################################################