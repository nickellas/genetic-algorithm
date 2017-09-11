""" Utilities and object for Agglomerative Hierarchical Clustering of atomic
    structures.
"""
from copy import copy
from scipy.cluster.hierarchy import linkage as get_linkage_matrix
from scipy.cluster.hierarchy import fcluster, leaders
from scipy.spatial import distance
import numpy as np

def assign_to_cluster(a, centroids, comp, a_d=None, t=None, outliers=True):
    """ Assigns a structure to an existing cluster based on the
        minimum distance to the centroids of all clusters.

        --Input--

        a: atoms object
            Structure to be assigned.

        centroids: list
            Cluster centroids (calculated with the cut_tree function).

        comp: comparator object
            Object that calculates features and similarity with
            get_features and get_similarity.

        a_d: dictionary
            Dictionary with features. If not given, it will be calculated.

        t: float
            Maximum distance to a centroid before assigning structure as
            outlier.

        outliers: boolean
            Include possibility of outlier assignment using threshold t.

        --Output--

        Cluster label (int).
    """
    if centroids is None:
        return None, None
    elif len(centroids) == 0:
        return 0, None

    if a_d is None:
        a_d = comp.get_features(a)

    # Delete feature types from structure
    # if they were excluded in the centroids
    for key in sorted(a_d.keys()):
        if key not in centroids[0].keys():
            print 'Features, {}, are not in centroids (deleted)'.format(key)
            del a_d[key]

    # If only one centroid, structure belongs to that cluster
    if len(centroids) == 1 and t is None:
        cluster_label = 1
        traj_centroid = comp.get_similarity(a_d, centroids[0])
        return cluster_label, [traj_centroid]

    if t is None:
    # Calculate distances between centroids of clusters
        nn_dist = []
        for c1 in centroids:
            c1_c2 = []
            for c2 in centroids:
                if all(np.hstack([i == j for i, j
                                  in zip(c1.values(), c2.values())])):
                    continue
                c1_c2.append(comp.get_similarity(c1, c2))
            nn_dist.append(min(c1_c2))

        max_dist = [max(nn_dist)]*len(centroids)
    else:
        max_dist = t

    # Return cluster label corresponding to the minimum distance
    a_c = []
    for i,c in enumerate(centroids):
        traj_centroid = comp.get_similarity(a_d, c)
        a_c.append(traj_centroid/t[i])
    if min(a_c) < 1. or outliers is False:
        cluster_label = a_c.index(min(a_c)) + 1
    else:
        cluster_label = 0  # Outlier (not belonging to any cluster)

    return cluster_label

class AgglomerativeClustering(object):
    """ Class object used for clustering of atomic structures.
        Clustering algorithm: Agglomerative Hierarchical Clustering
        (SciPy implementation).

        --Input--

        comp: comparator object
            Object that calculates features and similarity with
            get_features and get_similarity (and determines structures
            uniqueness with looks_like if unique_data=True).

        linakge: str
            Method for calculating distances between clusters.

        unique_data: boolean
            Eliminate duplicates in provided data with the comparator function
            looks_like.

        For details on the structure of the linkage matrix and the linkage
        methods, see the documentation on the scipy.cluster.hierarchy.linkage
        module.
    """

    def __init__(self, comp=None, linkage='average', unique_data=True):
        self.comp = comp
        self.linkage = linkage
        self.unique_data = unique_data
        self.data = []
        self.n_data = len(self.data)
        self.feature_matrix = None
        self.similarity_matrix = None
        self.linkage_matrix = None

    def grow_tree(self, data=None, feature_matrix=None,
                  similarity_matrix=None, linkage_matrix=None):
        """ Performs the clustering based on either a
            feature matrix (n_data,n_features) or a
            similarity matrix (n_data,n_data).
            The output is a linkage matrix which can be used for plotting
            a dendrogram with scipy.cluster.hierarchy.dendrogram.

            --Input--

            data: list
                List of n_data data to be clustered.

            feature_matrix: array, optional
                A matrix with size (n_data,n_features). If not given, it is
                computed with the comparator.

            similarity_matrix: array, optional
                A symmetric matrix with size (n_data,n_data) containing the
                similarities between all data.

            linkage_matrix: array, optional
                Matrix containing the result of the clustering.

            If neither a feature matrix nor a similarity matrix has been
            provided, they will both be computed when calling this function.
        """
        if data is not None and self.unique_data:
            self._get_unique_data_()
        elif len(self.data) == 0:
            self.data = data
        self.feature_matrix = feature_matrix
        self.similarity_matrix = similarity_matrix
        self.linkage_matrix = linkage_matrix

        if linkage_matrix is not None:
            self.n_data = int(linkage_matrix[-1, -1])
            self.linkage_matrix = linkage_matrix
            return feature_matrix, similarity_matrix, linkage_matrix

        if self.feature_matrix is None and self.similarity_matrix is None:
            # Compute feature matrix
            self.get_feature_matrix()

        if self.similarity_matrix is not None:
            Y = distance.squareform(self.similarity_matrix)
        elif self.feature_matrix is not None:
            # Compute similarity matrix
            self.get_similarity_matrix()
            Y = distance.squareform(self.similarity_matrix)

        print 'Computing linkage matrix...'
        # Compute linkage matrix
        linkage_matrix = get_linkage_matrix(Y, self.linkage)
        self.linkage_matrix = linkage_matrix
        self.n_data = int(linkage_matrix[-1, -1])

        return feature_matrix, similarity_matrix, linkage_matrix


    def cut_tree(self, t=None, criterion='inconsistent', depth=None,
                 cluster_min=3):
        """ Groups data into clusters based on the linkage matrix.

            --Input--

            t: float
                Threshold to be used by the criterion.

            criterion: str
                Criterion for grouping (cutting) the dendrogram.
                Can be either 'distance' or 'inconsistency'.

            depth: int
                Depth used when calculating inconsistency coefficient of a
                branch. Uses all banches by default.

            cluster_min: int
                Minimum cluster size. Data in clusters below this size will be
                assigned as outliers.

            --Output--

            labels: list
                List of cluster labels (integers) for each data.

            clusters: list
                Atoms objects grouped after cluster labels.
                The first group (list) is all the outliers.

            branches: list
                The branch index of each cluster.

            centroids: list
                The centroids of all clusters calculated as averaged features.

            cluster_energies: list
                The average energy of outliers (first element - empty if no
                outliers) and structures in the clusters. Is an empty list if
                data has no .get_potential_energy() attribute.

                The first energy is the average energy of the outliers. This
                is an empty list if there are no outliers.

            avg_width: float
                Average cluster width. Can be used as an outlier threshold in
                the assign_to_cluster function.
        """
        if t is None:
            if criterion == 'distance':
                t = 0.7 * max(self.linkage_matrix[:, 2])
            elif criterion == 'inconsistent':
                t = 4.0
        if depth is None:
            depth = self.n_data

        labels = fcluster(self.linkage_matrix, t, criterion, depth)
        branches = leaders(self.linkage_matrix, labels)

        # Data in clusters below cluster_size are outliers with label = 0
        for label in sorted(set(labels)):
            cluster_size = sum([i == label for i in labels])
            if cluster_size < cluster_min:
                labels = np.where(labels == label, 0, labels)

                x = np.delete(branches[0], np.where(branches[1] == label))
                y = np.delete(branches[1], np.where(branches[1] == label))
                branches = (x, y)

        # Rearrange labels to fill gaps from assignment of outliers
        for label in sorted(set(labels)):
            if label < 2:
                continue
            while len(np.where(labels == label-1)[0]) == 0:
                y = np.where(branches[1] == label, label-1, branches[1])
                branches = (branches[0], y)
                labels = np.where(labels == label, label-1, labels)
                label -= 1

        n_clusters = max(labels)
        print 'Number of clusters: {}'.format(n_clusters)

        # Group data of equal cluster label and calculate centroid
        clusters = [[] for i in range(n_clusters+1)]
        centroids = [[] for i in range(n_clusters)]
        i = 0
        for label, sample in zip(labels, self.data):
            clusters[label].append(sample)
            if label > 0:
                centroids[label-1].append(self.feature_matrix[i])
            i += 1

        print 'Number of outliers: {}'.format(len(clusters[0]))

        for i, c in enumerate(centroids):
            c_mean = dict()
            for key in sorted(self.feature_matrix[0].keys()):
                c_mean[key] = np.mean([x[key] for x in c], axis=0)
            centroids[i] = c_mean

        # Determine average width of clusters
        cluster_widths = []

        for branch in branches[0]:
            width = self.linkage_matrix[branch-self.n_data][-2]
            cluster_widths.append(width)
        print 'Cluster widths: {}'.format(cluster_widths)

        # Calculate average cluster energies
        cluster_energies = []
        for cluster in clusters:
            try:  # Check if data has a potential energy
                self.data[0].get_potential_energy()
            except: #IndexError, SyntaxError:
                mean_energy = 0
                cluster_energies.append(mean_energy)
                continue
            if len(cluster) == 0:
                mean_energy = []
            else:
                mean_energy = np.mean([x.get_potential_energy()
                                       for x in cluster])
            cluster_energies.append(mean_energy)

        return labels, clusters, branches, centroids, cluster_energies, \
            cluster_widths

    def get_leaves(self, branch_id):
        """ Returns indices of all the sub-branches (twigs) and data (leaves)
            of a branch with index branch_id. Can be used for coloring
            clusters in a dendrogram using the returned branch indices from
            cut_tree.
        """
        twigs = []
        # tree: matrix with all branch indices
        tree = np.array(np.delete(self.linkage_matrix, [2, 3], 1),dtype=int)
        branch = [branch_id]
        while True in [x >= self.n_data for x in branch]:
            leaves = copy(branch)
            for x in leaves:
                if x >= self.n_data:
                    branch = np.delete(branch, np.where(branch == x)[0])
                    branch = np.append(branch, tree[x-self.n_data])
                    twigs.append(x)

        leaves = [int(leaf) for leaf in branch]
        return leaves, twigs

    def get_similarity_matrix(self, feature_matrix=None):
        """ Computes similarity matrix from self.feature_matrix or an input
            feature matrix.
        """
        print 'Computing similarity matrix...'
        if feature_matrix is None:
            feature_matrix = self.feature_matrix
        n_data = len(feature_matrix)

        similarity_matrix = np.zeros((n_data, n_data))
        for i, a1 in enumerate(feature_matrix):
            for j, a2 in enumerate(feature_matrix[i+1:]):
                j = j+i+1
                similarity = self.comp.get_similarity(a1, a2)
                similarity_matrix[i][j] = similarity
                similarity_matrix[j][i] = similarity

        self.similarity_matrix = similarity_matrix
        return self.similarity_matrix

    def get_feature_matrix(self, data=None):
        """ Computes feature matrix from self.data or input data.
        """
        if data != None:
            self.data = [a.copy() for a in data]
            self._get_unique_data_()

        print 'Computing feature matrix...'
        feature_matrix = []
        for a in self.data:
            f = self.comp.get_features(a)
            feature_matrix.append(f)

        self.feature_matrix = feature_matrix
        return self.feature_matrix

    def add_data(self, data):
        """ Add data and update similarity matrix.
        """
        assert self.similarity_matrix is not None

        for a in data:
            feature = self.comp.get_features(a)

            idx = self.n_data

            s_row = np.zeros(self.n_data)
            for i, f2 in enumerate(self.feature_matrix):
                s_row[i] = self.comp.get_similarity(feature, f2)

            self.similarity_matrix = np.insert(self.similarity_matrix, idx,
                                               s_row, axis=0)
            s_col = np.insert(s_row, idx, 0.)
            self.similarity_matrix = np.insert(self.similarity_matrix, idx,
                                               s_col, axis=1)

            if self.feature_matrix is not None:
                self.feature_matrix = np.insert(self.feature_matrix, idx, feature)


            self.data.append(a)
            self.n_data += 1

    def _get_unique_data_(self):
        """ Remove duplicates in data with the comparator.
        """
        unique_data = []
        for a in self.data:
            dublicate = False
            for b in unique_data:
                if self.comp.looks_like(a, b):
                    dublicate = True
                    break
            if not dublicate:
                unique_data.append(a)

        self.data = unique_data
