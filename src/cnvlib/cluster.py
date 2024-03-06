#!/usr/bin/env python
"""Markov-cluster control samples' correlation matrix.

See:

    Stijn van Dongen, Graph Clustering by Flow Simulation,
    PhD thesis, University of Utrecht, May 2000.
    https://micans.org/mcl/

"""
import logging
import statistics

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from kneed import KneeLocator
import plotly.figure_factory as ff
import numpy as np
from hmmlearn import hmm
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
    K = range(2, 15)
    for k in K:
        km = KMeans(n_clusters=k)
        km = km.fit(depth_values)
        sum_of_squared_distances.append(km.inertia_)

    kn = optimal_clusters_size_knee(K, sum_of_squared_distances)
    return kn.knee
def kmeans_clustering(depth_values, no_of_clusters=None):
    depth_values = np.clip(depth_values, a_min=0, a_max=600)
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

    score = silhouette_score(depth_values, labels, metric="euclidean")

    stdev = []
    clusters = []
    for i in range(len(u_labels)):
        clusters.append([int(a) for a, b in zip(hp, labels) if b == i])
        stdev.append(np.std([(a, b) for a, b in zip(hp, labels) if b == i]))

    return centers, stdev, score

def optimal_clusters_size_knee(K, sum_of_squared_distances):
    return KneeLocator(x=K, y=sum_of_squared_distances, curve='convex', direction='decreasing')

def optimal_clusters_size_silhouette_score(depth_values, km):
    score = silhouette_score(depth_values, km.labels_, sample_size=2000)
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
def hmm_validation_custom(depth, lengths,  means):
    startprob = 1. / len(means) * np.ones(len(means))
    transmat = np.diag(np.full(len(means),1))
    means = [1, 46, 88, 132, 200]
    means_new = np.array(list([[i] for i in means]))
    covars = []
    for i in range(len(means)):
        if i == 0:
            covars.append(5)
        else:
           covars.append(20)
    covars_new = np.array(list([[i] for i in covars]))

    #covars_new = np.array([l[0] for l in x])
    #covars_new = 400 * np.tile(np.identity(1), (len(means), 1))
    #covars_new = 400 * np.tile(np.identity(2), (len(means), 1, 1))
    #means_new = np.array([[[1]], [[46]], [[88]], [[132]], [[180]]])
    covars_new = np.array([[[5]], [[20]], [[20]], [[20]], [[50]]])

    gen_model = hmm.GaussianHMM(n_components=len(means), covariance_type="full")
    gen_model.startprob_ = startprob
    gen_model.transmat_ = transmat
    gen_model.means_ = means_new
    gen_model.covars_ = covars_new

    gen_model.fit(depth, lengths)
    states = gen_model.predict(depth)
    score = silhouette_score(depth, states, metric="euclidean")
    means = gen_model.means_
    covars = gen_model.covars_
    return means, covars, score




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
    means = []
    stdev = []
    models = []
    scores = []

    rs = {}
    repeats = False
    for K in range(minK, maxK + 1):
        # print(K, datetime.now())
        if repeats == True:
            break
        my_best_ll = -1 * np.inf
        my_best_labels = None
        my_best_model = None
        for s in range(restarts):
            if repeats == True:
                break
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
                    init_params='m',
                    params='smc',
                    covariance_type=covar,
                    random_state=s,
                )
            elif tmat == 'free':
                model = hmm.GaussianHMM(
                    n_components=K,
                    init_params='m',
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

            x = [[[50]] for i in range(K)]
            #np.array(x).shape
            #np.array([l[0] for l in x]).shape

            model.startprob_ = np.ones(K) / K
            model.transmat_ = A
            #model.covars_ = np.array([l[0] for l in x])
            #model.covars_ = 400 * np.tile(np.identity(1), (K, 1))

            model.fit(X, lengths)

            #print(K, s)

            prob, labels = model.decode(X, lengths, algorithm=decode_alg)
            #print(sorted(list(model.means_)), model.n_components)

            #res = np.concatenate(sorted(list(model.means_))).ravel().tolist()
            # for x in range(len(res) - 1):
            #     if (res[x + 1] - res[x] < 15):
            #         repeats = True

            #print(res, model.n_components)

            if prob > my_best_ll:
                my_best_labels = labels
                my_best_ll = prob
                my_best_model = model

        means.append(np.concatenate(sorted(list(my_best_model.means_))).ravel().tolist())
        stdev.append(np.concatenate(sorted(list(my_best_model.covars_))).ravel().tolist())
        models.append(my_best_model)

        score = silhouette_score(X, my_best_labels, metric="euclidean")
        scores.append(score)
        print(score, K, np.concatenate(sorted(list(my_best_model.means_))).ravel().tolist())
        #print(score, K, np.concatenate(sorted(list(my_best_model.covars_))).ravel().tolist())

        add_scatter_trace_coverage(None, K, X, my_best_labels, np.concatenate(list(my_best_model.means_)).ravel().tolist(), np.concatenate(list(my_best_model.covars_)).ravel().tolist())

        rs[K] = my_best_ll, score, my_best_labels
        if score > best_score:
            best_score = score
            best_model = my_best_model
            best_labels = my_best_labels
            best_K = K

    # for i in range(len(scores)-1):
    #     best_score_final = scores[i]
    #     means_diff = [j-i for i, j in zip(means[i][:-1], means[i][1:])]
    #     print(means_diff)
    #     best_model = models[i]
    #     if abs(scores[i] - scores[i+1]) > 0.08 and all(i >= 12 for i in means_diff):
    #         break

    print(sorted(list(best_model.means_)), best_model.n_components)

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
    #plt.savefig("data/tmat.pdf", format="pdf", bbox_inches="tight")

    #plot_clusters_means(X, lengths)

    #print(means)
    #return means[-1], stdev[-1]
    return best_model.means_, best_model.covars_


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

def plot_clusters_means(X, lengths):
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from scipy.optimize import curve_fit

    # data generation
    np.random.seed(123)
    data = np.concatenate((np.random.normal(1, .2, 5000), np.random.normal(1.6, .3, 2500)))
    y, x, _ = plt.hist(X, alpha=.3, label='data')
    #x = (x[1:] + x[:-1]) / 2  # for len(x)==len(y)

    # x, y inputs can be lists or 1D numpy arrays

    def gauss(x, mu, sigma, A):
        res = A * np.exp(-(x - mu) ** 2 / 2 / sigma ** 2)
        return res

    def bimodal(x, mu1, sigma1, A1, mu2, sigma2, A2, mu3, sigma3, A3, mu4, sigma4, A4, mu5, sigma5, A5, mu6, sigma6, A6):
        res = gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2) + gauss(x, mu3, sigma3, A3) + gauss(x, mu4, sigma4, A4) + \
              gauss(x, mu5, sigma5, A5) + gauss(x, mu6, sigma6, A6)
        return res

    expected = (1, .2, 250, 2, .2, 125)
    params, cov = curve_fit(bimodal, x, y, expected)
    sigma = np.sqrt(np.diag(cov))
    x_fit = np.linspace(x.min(), x.max(), 500)
    # plot combined...
    plt.plot(x_fit, bimodal(x_fit, *params), color='red', lw=3, label='model')
    # ...and individual Gauss curves
    plt.plot(x_fit, gauss(x_fit, *params[:3]), color='red', lw=1, ls="--", label='distribution 1')
    plt.plot(x_fit, gauss(x_fit, *params[3:]), color='red', lw=1, ls=":", label='distribution 2')
    # and the original data points if no histogram has been created before
    # plt.scatter(x, y, marker="X", color="black", label="original data")
    plt.legend()
    print(pd.DataFrame(data={'params': params, 'sigma': sigma}, index=bimodal.__code__.co_varnames[1:]))
    plt.show()


def hmm_pome_test(cnarr, X, means, stdev, snps_cpd_means_input):
    import scipy
    import pomegranate as pom
    import scipy.special

    #X = np.concatenate(list(X)).ravel().tolist()

    from src.cnvlib.segmentation.hmm import observations_matrix
    observations = observations_matrix(cnarr, snps_cpd_means_input)

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

    score = silhouette_score(X, states, metric="euclidean")

    return means, stdev, score

def add_scatter_trace_coverage(chrom, depth_values_hp1, depth_values_hp2, labels, means, covar):
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    import math
    #fig = go.Figure()

    depth_values_hp1 = np.array(depth_values_hp1, dtype='int')
    depth_values_hp2 = np.array(depth_values_hp2, dtype='int')

    depth_values_hp1 = depth_values_hp1.reshape(-1, 1)
    depth_values_hp2 = depth_values_hp2.reshape(-1, 1)
    Y = np.concatenate([depth_values_hp1, depth_values_hp2])

    fig = make_subplots(rows=1, cols=2, shared_yaxes=True, column_widths=[0.8, 0.2], vertical_spacing=0.01, horizontal_spacing=0.01)

    y = np.concatenate(list(Y)).ravel().tolist()
    lst = [i for i in range(0, (len(y) // 2) * 50000, 50000)]
    x = lst + lst

    cdict = {0: '#1f77b4',  # muted blue
    1: '#ff7f0e',  # safety orange
    2: '#2ca02c',  # cooked asparagus green
    3: '#d62728',  # brick red
    4: '#9467bd',  # muted purple
    5: '#8c564b',  # chestnut brown
    6: '#e377c2',  # raspberry yogurt pink
    7: '#7f7f7f',  # middle gray
    8: '#bcbd22',  # curry yellow-green
    9: '#17becf'   # blue-teal
                 }
    ldict = {0: 'Cluster_1', 1: 'Cluster_2', 2: 'Cluster_3', 3: 'Cluster_4', 4: 'Cluster_5',\
             5: 'Cluster_6', 6: 'Cluster_7', 7: 'Cluster_8', 8: 'Cluster_9', 9: 'Cluster_10'}

    for g in np.unique(labels):
        ix = [index for index, i in enumerate(x) if labels[index] == g]
        xn = [x[i] for i in ix]
        yn = [y[i] for i in ix]
        fig.add_trace(go.Scatter(x=xn, y=yn, mode='markers', marker=dict(color=cdict[g]), name=ldict[g], opacity=0.7,),
                      row=1, col=1)

    for i in range(len(means)):
        fig.add_hline(y=means[i], line_width=2,
                  line=dict(dash='solid'), line_color=cdict[i], annotation_text="mean_"+str(i+1))

        fig.add_hline(y=means[i] + covar[i], line_width=2,
                      line=dict(dash='dash'), line_color=cdict[i], annotation_text="stdev_"+str(i+1))
        fig.add_hline(y=means[i] - covar[i], line_width=2,
                      line=dict(dash='dash'), line_color=cdict[i], annotation_text="stdev_" + str(i + 1))
    fig.update_yaxes(range=[0, 150])

    #fig.add_trace(go.Histogram(y=y, orientation='h', nbinsy=5000,), row=1, col=2)
    fig.add_trace(go.Histogram(y=np.concatenate(list(Y)).ravel().tolist(), orientation='h', marker_color='gray', name='Histogram', nbinsy=8000), row=1, col=2)
    fig.update_layout(xaxis2=dict(range=[0, 200]))
    #fig.add_trace(go.Histogram(y=np.concatenate(list(Y[0:len(Y)//2-1])).ravel().tolist(), name='HP-1', orientation='h', nbinsy=8000, marker_color='#6A5ACD'), row=1, col=2)
    #fig.add_trace(go.Histogram(y=np.concatenate(list(Y[len(Y)//2:len(Y)])).ravel().tolist(), name='HP-2', orientation='h', nbinsy=8000, marker_color='#2E8B57'), row=1, col=2)

    # Overlay both histograms
    #fig.update_layout(barmode='overlay')
    # Reduce opacity to see both histograms
    #fig.update_traces(opacity=0.75)

    score = silhouette_score(Y, labels, metric="euclidean")
    print(score)

    fig.write_html("coverage_plots/"+chrom+"_cluster_"+str(K)+".html")

def cpd_mean_hps(haplotype1_means, haplotype2_means):
    import ruptures as rpt
    import numpy as np
    import statistics
    from random import randint

    data = np.array(haplotype1_means, dtype='uint8') #numpy.clip(haplotype1_means, a_min=1, a_max=300)
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
        sub_list=haplotype1_means[start:point]
        if len(sub_list) < 25:
            continue
        else:
            count = 0
            sub_list = [int(i) for i in sub_list]
            mode_hp1 = max(set(sub_list), key=sub_list.count)
            means_hp1 = statistics.mean(sub_list)
            means_hp1 = mode_hp1 if mode_hp1 > means_hp1  else means_hp1
            means_hp1 = statistics.median(sub_list)
            #means_hp1 = statistics.mean(sub_list)
            for i in range(len(sub_list)):
                if means_hp1 - 20 <= sub_list[i] <= means_hp1 + 20:
                    count += 1
            if count > len(sub_list)//3:
                snps_haplotype1_mean.append(means_hp1)
                snps_haplotype1_len.append(len(sub_list))

            if count > len(sub_list)/1.1:
                strong_candidates.append(int(means_hp1))

            for i in range(len(sub_list)):
                if sub_list[i] >= means_hp1 + 10 or sub_list[i] <= means_hp1 - 10:
                    sub_list[i] = randint(int(means_hp1)-10, int(means_hp1)+10)
            haplotype1_means[start:point] = [sub_list[i] for i in range(len(sub_list))]

        start = point + 1
        snps_haplotype1_pos.append(point * 50000)
    ############################################################
    data = np.array(haplotype2_means, dtype='uint8') #numpy.clip(haplotype2_means, a_min=1, a_max=300)
    algo = rpt.Pelt(model="rbf").fit(data)
    result = algo.predict(pen=10)
    change_points = [i for i in result if i <= len(data)]

    snps_haplotype2_mean = []
    snps_haplotype2_pos = []
    snps_haplotype2_len = []

    start = 0
    snps_haplotype2_pos.append(0)
    for index, point in enumerate(change_points):
        sub_list=haplotype2_means[start:point]
        if len(sub_list) < 25:
            continue
        else:
            count = 0
            sub_list = [int(i) for i in sub_list]
            mode_hp2 = max(set(sub_list), key=sub_list.count)
            means_hp2 = statistics.mean(sub_list)
            means_hp2 = mode_hp2 if mode_hp2 > means_hp2  else means_hp2
            means_hp2 = statistics.median(sub_list)
            for i in range(len(sub_list)):
                if means_hp2 - 20 <= sub_list[i] <= means_hp2 + 20:
                    count += 1
            if count > len(sub_list)//3:
                snps_haplotype2_mean.append(means_hp2)
                snps_haplotype2_len.append(len(sub_list))

            if count > len(sub_list)/1.1:
                strong_candidates.append(int(means_hp2))

            for i in range(len(sub_list)):
                if sub_list[i] >= means_hp2 + 10 or sub_list[i] <= means_hp2 - 10:
                    sub_list[i] = randint(int(means_hp2)-10, int(means_hp2)+10)
            haplotype2_means[start:point] = [sub_list[i] for i in range(len(sub_list))]

        start = point + 1
        snps_haplotype2_pos.append(point * 50000)

    # print(snps_haplotype1_mean, snps_haplotype1_len)
    # print(snps_haplotype2_mean, snps_haplotype2_len)

    # return [i for i in final if i >= 1]

    return snps_haplotype1_mean, snps_haplotype1_len, snps_haplotype2_mean, snps_haplotype2_len, haplotype1_means, haplotype2_means, strong_candidates

