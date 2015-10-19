"""
    plotting
    ========

    This module implements the plotting and visualising used in
    the sam_gmt package.

"""

from __future__ import division
import itertools
import operator

import numpy as np
from matplotlib import pyplot as plt

import stats


def aligned_reads(groups, numfastq):
    """ Return a list of mapped_reads and the number of reads in the
        fastq file.

        :param groups: An iterable of 3(?) groups.
        :param numfastq: The number of (?)
    """
    mapped = sum(map(len, groups)) / numfastq

    labels = ("p<1e-3", "1e-3<=p<1e-2", "1e-2<=p<1", "Unmapped")
    x = [len(group) / numfastq for group in groups] + [1 - mapped]

    plt.figure(1, figsize=(8, 8))
    plt.axes([0.1, 0.1, 0.8, 0.8])

    plt.pie(x, labels=labels, autopct="%1.1f%%", shadow=True)
    plt.title("Mapping Statistics")
    plt.show()

    return mapped


def base_composition(reads, base):
    """ Plot the compositon of a base.
    """
    assert base.upper() in set("ACGT")

    """ Reports nucelotide frequencies at each position in the
        sam sequences
    """
    # DNA_Alphabet=["A","C","T","G","N"]
    all_nucs = []
    for read in reads:
        nucs = {}  # Dictionary to store nucleotide data.
        seq = read[9]
        for i in range(0, len(seq)):
            nucs[str(i + 1)] = seq[i]
        all_nucs.append(nucs)
    all_items = []
    counts = []
    for dicts in all_nucs:
        for item in dicts.items():
            all_items.append(item)
    all_items.sort(key=operator.itemgetter(0))
    groups = [map(operator.itemgetter(1), list(group))
              for key, group in itertools.groupby(
                  all_items, operator.itemgetter(0))]
    for group in groups:
        counts.append(group.count(base))

    pos = range(1, len(seq) + 1)

    # Create plot.
    plt.figure(1, figsize=(8, 8))
    plt.axes([0.1, 0.1, 0.8, 0.8])
    plt.bar(pos, counts, facecolor='g')
    plt.xlabel("Position")
    plt.ylabel("number of mapped reads")
    plt.title(base)
    plt.show()


def length_distrib(group, name):
    """ Distribution of lengths of all the sam reads.
    """
    lengths = []
    for read in group:
        if int(read[8]) < 0:
            length = -1 * int(read[8])
        else:
            length = int(read[8])
        lengths.append(length)
    # Visualize length distribution
    plt.figure(1, figsize=(8, 8))
    plt.axes([0.1, 0.1, 0.8, 0.8])
    n, bins, patches = plt.hist(lengths, 100, normed=0, facecolor='g')
    plt.xlabel("lengths")
    plt.ylabel("number of mapped reads")
    plt.title(name)
    plt.show()


def rpkm_repr(rpkm_data, timepoint):
    """ Plot showing level of agreement between technical replicates
        for RPKM between replicates and plots coefficient of
        determination.
    """
    one = []
    two = []
    if timepoint == "t1":
        for i in range(0, len(rpkm_data.values())):
            one.append(int(rpkm_data.values()[i][0]))
            two.append(int(rpkm_data.values()[i][1]))
    else:
        for i in range(0, len(rpkm_data.values())):
            one.append(int(rpkm_data.values()[i][2]))
            two.append(int(rpkm_data.values()[i][3]))
    plt.plot(one, two, 'o')
    pcc = stats.pearson_correlation(one, two)
    R2 = pcc ** 2
    name = """Technical Replicates\nR2=""" + str(R2)
    m, b = np.polyfit(one, two, 1)
    plt.figure(1, figsize=(8, 8))
    plt.plot(one, np.array(one) * m + b, 'r-')
    plt.text(3000, max(two) - 1000, name, fontsize=12)
    plt.xlabel("RPKM replicate 1")
    plt.ylabel("RPKM replicate 2")
    plt.title(timepoint)
    plt.show()


def box_plots(norm, original):
    """ Distribution of the coeficient of variation across samples
        (replicates) normalised using the methods provided
    """
    bp = plt.boxplot([norm, original], notch=False, patch_artist=True)
    for box in bp['boxes']:
        box.set(color="red")
        box.set(color="blue")
    plt.ylabel("coefficient of variation")
    plt.xlabel("Methods")
    my_xticks = ['RPKM', 'raw counts']
    x = [1, 2]
    plt.xticks(x, my_xticks)
    plt.ylim(0, 400)
    plt.show()


def avgerage_cv(norm, original):
    """ Distribution of the coeficient of variation across samples
        (replicates) normalised using the methods provided
    """
    x = [1, 2]
    y = [np.mean(norm), np.mean(original)]
    plt.figure(1, figsize=(8, 8))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    plt.bar(x[0], y[0], color="red", label="RPKM")
    plt.bar(x[1], y[1], color="blue", label="Raw counts")
    plt.ylabel("Average coefficient of variation")
    plt.xlabel("Methods")
    ax.xaxis.set_ticklabels([])
    plt.legend(loc="upper right")
    plt.show()


def ma(rpkm_data, cutoff=[-1.5, 1.5]):
    logfc = []
    avg_rpkm = []
    sig_logfc = []
    sig_avg_rpkm = []
    logfc2 = []
    avg_rpkm2 = []
    sig_logfc2 = []
    sig_avg_rpkm2 = []
    for i, ii, s, ss in rpkm_data.values():
        fc = np.log2(float(s + 1) / (i + 1))
        if fc < cutoff[0] or fc > cutoff[1]:
            sig_logfc.append(fc)
            sig_avg_rpkm.append(np.log2(s + 1) + np.log2(i + 1) / 2)
        else:
            logfc.append(fc)
            avg_rpkm.append(np.log2(s + 1) + np.log2(i + 1) / 2)
    for i, ii, s, ss in rpkm_data.values():
        fc2 = np.log2(float(ss + 1) / (ii + 1))
        if fc2 < cutoff[0] or fc2 > cutoff[1]:
            sig_logfc2.append(fc2)
            sig_avg_rpkm2.append(np.log2(ss + 1) + np.log2(ii + 1) / 2)
        else:
            logfc2.append(fc2)
            avg_rpkm2.append(np.log2(ss + 1) + np.log2(ii + 1) / 2)
    plt.figure(1, figsize=(8, 8))
    plt.axes([0.1, 0.1, 0.8, 0.8])
    plt.plot(avg_rpkm, logfc, 'o', color="blue", label="rep1")
    plt.plot(avg_rpkm2, logfc2, 'x', color="blue", label="rep2")
    plt.plot(sig_avg_rpkm, sig_logfc, 'o', color="red", label="sig rep1")
    plt.plot(sig_avg_rpkm2, sig_logfc2, 'x', color="red", label="sig rep2")
    plt.axhline(cutoff[0], color="orange")
    plt.axhline(cutoff[1], color="orange")
    plt.ylabel("Fold Change (log2)")
    plt.xlabel("Average RPKM (log2)")
    plt.title("MA plot")
    plt.legend(loc="upper left")
    plt.show()


def ma_pval(rpkm_data, cutoff=0.05):
    logfc = []
    avg_rpkm = []
    sig_logfc = []
    sig_avg_rpkm = []
    for i, ii, s, ss, pval in rpkm_data.values():
        fc = np.log2(float(s + 1) / (i + 1))
        if float(pval) < cutoff:
            sig_logfc.append(fc)
            sig_avg_rpkm.append(np.log2(s + 1) + np.log2(i + 1) / 2)
        else:
            logfc.append(fc)
            avg_rpkm.append(np.log2(s + 1) + np.log2(i + 1) / 2)
    plt.figure(1, figsize=(8, 8))
    ax = plt.axes([0.1, 0.1, 0.8, 0.8])
    plt.plot(avg_rpkm, logfc, 'o', color="blue", label="rep1")
    plt.plot(sig_avg_rpkm, sig_logfc, 'o', color="red", label="sig rep1")
    plt.ylabel("Fold Change (log2)")
    plt.xlabel("Average RPKM (log2)")
    plt.title("MA plot")
    plt.legend(loc="upper left")
    plt.show()
