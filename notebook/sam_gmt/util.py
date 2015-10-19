"""
    util
    ====

    This submodule implements additional functions that
    are useful but not categorized withing the other submodules.
"""

from __future__ import division
from collections import defaultdict

import numpy as np


def get_cv(data, condition):
    """ Returns a list of cv data.
    """
    slice_index = 0 if condition == "t1" else 2
    slc = slice(slice_index, slice_index + 2)

    cvs = []

    for val in data.values():
        mean = np.mean(val[slc])
        std = np.std(val[slc])

        if mean != 0 or std != 0:
            cvs.append((mean + 1) / (std + 1))

    return cvs


def get_rpkm(data, num_maps):
    """ Return a dictionary of Reads per Kilobase per Million.
        Provide number of mapped reads for the two groups of interest
        and raw count data. This method provides length normalisation
        to prevent length and total count bias.
    """
    count_lists, lengths = zip(*data.values())[:-1], zip(*data.values())[-1]

    # Use a defaultdict from collections, unknow keys return a list.
    rpkms = defaultdict(list)

    for i, count_list in enumerate(count_lists):
        for count, length in zip(count_list, lengths):
            if count == 0:
                rpkm = 0
            else:
                rpkm = count / (length * (num_maps[i] / 1e6))
            rpkms[i].append(rpkm)

    final = {}

    for i, key in enumerate(data.keys()):
        result = [lst[i] for lst in rpkms.values()]
        final[key] = result

    return final


def strand_filter(reads, forward=True):
    """ Return a list of lines from a sam file read
        that were aligned to a particular strand.
        No accounting is done for duplicate lines.
    """
    if forward:
        return [read for read in reads if read[9] > 0]
    else:
        return [read for read in reads if read[9] < 0]


def sub_groups(mapped_reads):
    """ Return a tuple of three "groups",
        seperated by p<1e-3, 1e-3<=p<1e-2 and 1e-2<=p<1
    """
    groups = [], [], []

    for read in mapped_reads:
        read_val = int(read[4])
        if read_val > 29:
            group_index = 0
        elif 17 < read_val <= 29:
            group_index = 1
        elif read_val <= 17:
            group_index = 2

        groups[group_index].append(read)

    return groups
