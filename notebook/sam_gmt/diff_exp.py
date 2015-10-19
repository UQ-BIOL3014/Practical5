"""
    diff_exp
    ========

    This submodule implements the differential expression methods
    of the sam_gmt project.
"""

from __future__ import division

import numpy as np

from matplotlib import pyplot as plt
# from scipy.cluster import hierarchy
import scipy
from scipy import stats as scipy_stats


def welch_ttest(rpkm):
    """ Performs a one-tailed Welchs T-statistic test
        Returns a dictionary mapping identifiers to their original value
        lists, with the calculated metric appended to the end.
    """
    # Generate T scores.
    t_scores = []
    for i, ii, s, ss in rpkm.values():
        sd1 = np.std([i, ii]) / 2
        sd2 = np.std([s, ss]) / 2

        if sd2 == 0 and sd1 == 0:
            t = np.nan
        else:
            t = (np.mean([s, ss]) - np.mean([i, ii]) / np.sqrt(sd2 + sd1))
        t_scores.append(t)

    # Generate P values.
    pvals = []
    for i, t in enumerate(t_scores):
        pval = scipy_stats.t.sf(abs(t), 1)
        if not pval:
            pvals.append(1)
        else:
            pvals.append(pval)

    corr_pvals = correct_pvalues(pvals)

    # Assemble result.
    result = {}
    for (k, v), corr in zip(rpkm.items(), corr_pvals):
        result[k] = v + [corr]

    return result


def correct_pvalues(p_values, correction_type="Benjamini-Hochberg"):
    """ Consistent with R printcorrect_pvalues_for_multiple_testing(
            [0.0, 0.01, 0.029, 0.03, 0.031, 0.05,
            0.069, 0.07, 0.071, 0.09, 0.1])
    """
    pvals = np.array(p_values)
    n = float(pvals.shape[0])
    new_pvals = np.empty(n)

    if correction_type == "Bonferroni":
        #
        new_pvals = n * pvals

    elif correction_type == "Bonferroni-Holm":
        values = sorted((pval, i) for i, pval in enumerate(pvals))
        for rank, (pval, i) in enumerate(values):
            new_pvals[i] = (n - rank) * pval

    elif correction_type == "Benjamini-Hochberg":
        values = list(reversed(sorted((pvalue, i) for i, pvalue
                                      in enumerate(pvals))))

        new_vals = []

        for i, (pval, index) in enumerate(values):
            rank = n - i
            val = (n / rank) * pval
            new_vals.append(val)

        for i in range(int(n) - 1):
            if new_vals[i] < new_vals[i + 1]:
                new_vals[i + 1] = new_vals[i]

        for i, (pval, index) in enumerate(values):
            new_pvals[index] = new_vals[i]

    return new_pvals


def cluster_data(data_matrix, genenames, timepoint):
    """ One replicates a specific time. (?)
    """
    dist_matrix = np.zeros([np.shape(data_matrix)[0], 1])

    # Generate a distance matrix
    for i in range(np.shape(data_matrix)[0]):
        for j in (0, 1):
            dist_matrix[i, j] = abs(data_matrix[i] - data_matrix[j]) ** 2

    # Define labels.
    labels = ["{}, {}".format(*item) for item in genenames]

    # Create clustered objects.
    linked = scipy.hierarchy.linkage(dist_matrix, method="centroid")
    dend = scipy.hierarchy.dendrogram(linked,
                                      orientation="right",
                                      labels=labels)

    # Save figure.
    fig = plt.figure(1, figsize=(17, 8))
    plt.title(timepoint)
    fig.savefig("{}-dendrogram.png")

    return dend["ivl"]
