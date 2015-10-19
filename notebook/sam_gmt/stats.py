"""
    stats
    =====

    This submodule implements the statistical and numerical methods
    of the sam_gmt project.
"""

from __future__ import division
from collections import Counter

import numpy as np


def nucleotide_freqs(reads):
    """ Returns a dictionary of nucleotide percentages ("ACGT").
    """
    nucs = "".join(read[9] for read in reads)
    counts = Counter(nucs)

    freqs = {k: v / sum(counts.values()) for k, v in counts.items()}
    return freqs


def dinucleotide_freqs(reads):
    """ Returns a dictionary of dinucleotide frequencies using
        p(Rho) statistics for over representation.
    """
    single_freqs = nucleotide_freqs(reads)

    # Create list of consequtive dinucleotide pairs.
    di_nucleotides = [seq[i:i + 2] for seq in (read[9] for read in reads)
                      for i in range(len(seq) - 1)]

    counts = Counter(di_nucleotides)
    freqs = {k: (v / sum(counts.values())) / (single_freqs[k[0]]
                                              * single_freqs[k[1]])
             for k, v in counts.items()}

    return freqs


def mapped_reads(reads, paired_end=True):
    """ If duplicate tracking was enabled, then this attempts
        to recapitulate the number of unique, mapped, probe-id's in
        the original sam file. It is multiplied by 2 for paired-end
        data with duplicate read id's.
        The idea is that if you divide this by the number of reads in
        the fastq you aligned (possibly from the output of
        fastq-stats), you will get an accurate "percentage of reads
        aligned" statistic.
        "mapped" is something with a non-negative position, and a
        "non-asterisk" cigar string.
    """
    s_reads = [read for read in reads if read[3] > 0 and read[5] != '*']
    return s_reads

    # m_reads = [read[0] for read in reads if read[3] > 0 and read[5] != '*']
    ## Woo crazy boolean coercion!
    # mapped = len(set(m_reads)) * (paired_end + 1)


def mapped_bases(mapped_reads):
    """ Return the number of mapped bases in SAM file.
    """
    seq = "".join(read[9] for read in mapped_reads)
    return len(seq)


def length_stats(*args):
    """ Returns a dictionary of basic stats relating to the length of
        the reads.
        Calculations are based on the length of the (possibly
        hard-clipped) sequence in the same file.

        Three groups should be passed in as arguments:
            >>> length_stats(group1, group2, group3)
    """
    data = {}

    for i, group in enumerate(args):
        lengths = [abs(int(read[8])) for read in group]

        mean_len = np.mean(lengths)
        max_len = np.max(lengths)
        min_len = np.min(lengths)

        data[i] = (mean_len, max_len, min_len)

    return data


def pearson_correlation(x, y):
    """ Return a float between -1 and 1, representing the pearson
        correlation.
    """
    assert len(x) == len(y), "Input lists must be the same length."
    n = len(x)
    assert n > 0, "Input lists must have one or more elements."

    x_avg, y_avg = np.mean(x), np.mean(y)

    diff_product, x_diff_total, y_diff_total = 0, 0, 0

    for x_elem, y_elem in zip(x, y):
        x_diff = x_elem - x_avg
        y_diff = y_elem - y_avg

        diff_product += x_diff * y_diff

        x_diff_total += x_diff ** 2
        y_diff_total += y_diff ** 2

    return diff_product / np.sqrt(x_diff_total * y_diff_total)


def inv_logit(p):
    """ TODO: is this necessary as a function?
    """
    return 10 ** (p / -10)
