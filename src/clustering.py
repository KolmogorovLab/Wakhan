from sklearn.cluster import KMeans
from kneed import KneeLocator
from sklearn.metrics import silhouette_score
import plotly.figure_factory as ff
import numpy as np
import pandas as pd
from utils import flatten
from smoothing import smoothing

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
    depth_values = np.clip(depth_values, a_min=1, a_max=300)
    depth_values[depth_values != 0]
    hp = list(depth_values)
    depth_values = depth_values.reshape(-1, 1)

    if no_of_clusters == None:
        no_of_clusters = kmeans_clustering_seed(depth_values)
    #_, _, _, _, no_clusters = bayesian_gaussian_mixture_clustering(depth_values)
    #print(no_clusters)
    gmm(depth_values)

    km = KMeans(n_clusters=no_of_clusters, n_init=25, max_iter=600, random_state=0)
    km = km.fit(depth_values)
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
def plot_optimal_clusters(hp1, hp2, unphased, arguments):
    hp1 = flatten(hp1)
    hp2 = flatten(hp2)
    unphased = flatten(unphased)
    unphased, hp1, hp2 = smoothing(unphased, hp1, hp2, conv_window_size=15)

    hp1=np.clip(hp1, a_min=1, a_max=300)
    hp2=np.clip(hp2, a_min=1, a_max=300)

    hp1[hp1 != 0]
    hp1=hp1.reshape(-1, 1)

    hp2[hp2 != 0]
    hp2=hp2.reshape(-1, 1)
    depth_values = np.concatenate((hp1, hp2))

    u_labels, labels, centers, stdev, clusters = kmeans_clustering(depth_values, arguments['no_of_clusters'])
    plot_clsters_settings(clusters, arguments)

    if arguments['no_of_clusters'] == None:
        #nearest_clusters = []
        for no_of_clusters in range(len(clusters)-1, len(clusters)+2):
            u_labels, labels, centers, stdev, clusters = kmeans_clustering(depth_values, no_of_clusters)
            plot_clsters_settings(clusters, arguments)
def plot_clsters_settings(clusters, arguments):
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

    fig.write_html(arguments['out_dir_plots'] +'/'+ arguments['genome_name'] + "_genome_bins_clusters["+str(len(clusters))+"]_histogram.html")

def bayesian_gaussian_mixture_clustering(points, clouds=None, concentration_prior=None, K=10, restarts=5, seed=0):

    from sklearn.mixture import BayesianGaussianMixture
    from collections import Counter

    total = list(points)
    if clouds is not None:
        total.extend(list(clouds))
    npArray = np.array(total)

    gmm = BayesianGaussianMixture(
        n_components=K,
        n_init=restarts,
        weight_concentration_prior=concentration_prior,
        max_iter=int(1e6),
        random_state=seed,
    )
    targetAssignments = gmm.fit_predict(npArray)
    targetAssignments = targetAssignments[: len(points)]
    mus = gmm.means_
    sigmas = gmm.covariances_
    cntr = Counter(targetAssignments)
    numPoints = [cntr[i] if i in cntr else 0 for i in range(K)]
    numClusters = len(cntr)

    return mus, sigmas, targetAssignments, numPoints, numClusters

def gmm(X):
    from sklearn.mixture import GaussianMixture
    gm = GaussianMixture(n_components=5, random_state=0).fit(X)
    labels = gm.predict(X)
    probs = gm.predict_proba(X)

    labels = list(labels)
    u_labels = np.unique(labels)

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