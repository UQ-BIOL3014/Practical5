#!/usr/bin/env python

"""
    exp_redux
    =========

    This file executes the logic originally run globally at the
    bottom of the original sam.py file that the sam_gmt module is
    based off.

    Usage
    -----
        python perform_comparison.py <raw_count_file> [p-value cutoff (0-1)]

"""


from __future__ import division, print_function

import os
import sys

import numpy as np

import diff_exp, fileio, plotting, util


def run_differential(count_file, pval_cutoff=0.05, display_plots=False):
    """ Performs a differential expression analysis.
        Modeled off global code found in sam.py
    """
    # Scope variables.
    replicates = ["1", "10", "1_2", "10_2"]
    num_map = (118898, 121634, 136286, 135102)

    # Load in data.
    raw_data = fileio.read_rawcounts(count_file)

    # Normalize and Visualize.
    rpkm = util.get_rpkm(raw_data, num_map)

    # Get cv
    meth = util.get_cv(rpkm, "t1")
    orig = util.get_cv(raw_data, "t1")

    # Perform T test.
    result_ttest = diff_exp.welch_ttest(rpkm)

    if display_plots:
        for rep in replicates:
            plotting.rpkm_repr(rpkm, rep)

        # Visualize variation.
        plotting.box_plots(meth, orig)
        plotting.avgerage_cv(meth, orig)

        # Create MA plots
        plotting.ma(rpkm)
        plotting.ma_pval(result_ttest, 0.01)

    # Map genes to pvalues if significant.
    diff = {k: v[-1] for k, v in result_ttest.items()
            if np.logical_not(np.isnan(v[-1])) and v[-1] <= pval_cutoff}

    print("{} Differentially Expressed Genes Found".format(len(diff)))
    if diff:
        for gene, pval in sorted(diff.items(),
                                 key=lambda item: item[-1]):
            print("\t{:>10}: {:0.2e}".format(gene, pval))
    else:
        print("\t*NONE*")


def main():
    """ Entry point for the program.
        Simulates the global code from the original sam.py
    """

    if not 1 < len(sys.argv) < 4 or not os.path.exists(sys.argv[1]):
        print("Please specify a raw count file.")
        return 1

    pval_threshold = 0.05
    if len(sys.argv) == 3:
        try:
            pval_threshold = float(sys.argv[2])
            if not (0 <= pval_threshold <= 1):
                raise ValueError("pvalue '{}' not between 0 and 1".format(
                    pval_threshold))
        except Exception:
            print("p value threshold is not valid.",
                  "Must be a float between 0 and 1.")
            return 2

    count_file = sys.argv[1]
    run_differential(count_file, pval_threshold)

    return 0


if __name__ == '__main__':
    exit(main())
