#!/usr/bin/env python3


########################################################################
# File: pairwise_fisher.py
#  executable:
# Purpose:
#
#
# Author: Cameron M. Soulette
# History:      cms 01/08/2020 Created
#
########################################################################

########################################################################
# Hot Imports & Global Variable
########################################################################

import argparse
import sys
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests


########################################################################
# Helper Functions
#
#
########################################################################


def loadNPZ(x):
    """
    takes in npz formatted matrix.
    """

    try:
        data = np.load(x)
    except:
        print("ERR ** Cannot load matrix %s. Check path or format." % x)
        sys.exit(1)
    return data


def getColIndexFromArray(x, y):
    """
    takes in list of strings = x
    and finds list index in array = y
    """

    return np.nonzero(np.isin(y, x))


def returnSamplesFromManifest(x):
    """
    reads in mesa formatted manifest
    returns list of samples
    """
    s = list()
    with open(x) as fin:
        for i in fin:
            s.append(i.split()[0])

    return s


def getClust(fname):
    data = dict()
    with open(fname) as fin:
        for i in fin:
            left, right = i.rstrip().split()
            mxes = right.split(",")
            data[left] = mxes
    return data


########################################################################
# MAINE
#
#
########################################################################

def add_parser(parser):
    parser.add_argument(
        "--inclusionMESA",
        type=str,
        required=True,
        help="Compressed NPZ formatted Inclusion count matrix from quantMESA.",
    )
    parser.add_argument(
        "-c",
        "--clusters",
        type=str,
        required=True,
        help="Clusters table.",
    )
    parser.add_argument(
        "--chi2",
        action="store_true",
        default=False,
        help="Use X^2 instead of fishers. Quicker, not as sensitive.",
    )
    parser.add_argument(
        "--no-correction",
        action="store_true",
        default=False,
        help="Output raw p-values instead of corrected ones. Correction is "
        "done via Benjamini-Hochberg",
    )


def run_with(args):
    pmesa = args.inclusionMESA
    cmesa = args.clusters
    if args.chi2:
        test_method = chi2_contingency
    else:
        test_method = fisher_exact

    # load psi
    data = loadNPZ(pmesa)
    clusters = getClust(cmesa)

    # table has 3 arrays, cols, rows and data
    cols, rows, matrix = data["cols"], data["rows"], data["data"]

    comparisons = set()
    for i, v in enumerate(cols):
        for j, v in enumerate(cols):
            if i == j:
                continue
            comparisons.add(tuple(sorted([i, j])))

    # do the math
    comps = list(comparisons)

    print(
        "clusterID",
        "\t".join("%s_%s" % (cols[j[0]], cols[j[1]]) for j in comps),
        sep="\t",
    )

    for n, vals in enumerate(matrix):
        event_id = rows[n]
        mxes = matrix[np.isin(rows, clusters[event_id])]

        inc = vals
        exc = np.sum(mxes, axis=0)

        pvalues = list()

        for i in comps:
            left, right = i
            table = [[inc[left], inc[right]], [exc[left], exc[right]]]

            # chi^2 test will fail on contingency tables where any column or
            # rows are zero, but not if diagonal is zero
            # [[0, 0], [2, 3]] fail
            # [[0, 2], [0, 3]] fail
            # [[0, 2], [2, 0]] still valid
            if args.chi2:
                col_total = [inc[left] + inc[right], exc[left] + exc[right]]
                row_total = [inc[left] + exc[left], inc[right] + exc[right]]
                if any(count == 0 for count in row_total + col_total):
                    continue

            data = test_method(table)[1]
            pvalues.append(data)
        if not args.no_correction and len(pvalues) > 0:
            pvalues = multipletests(pvalues, method="fdr_bh")[1]
        print(event_id, "\t".join(str(x) for x in pvalues), sep="\t")


def main():
    parser = argparse.ArgumentParser()
    add_parser(parser)
    args = parser.parse_args()
    run_with(args)


if __name__ == "__main__":
    main()
