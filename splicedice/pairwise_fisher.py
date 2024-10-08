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
# Helper Functions
#
#
########################################################################


import numpy as np


def getClusters(filename, filter_list=None):
    clusters = {}
    with open(filename) as allClusters:
        for line in allClusters:
            try:
                event, mxes = line.rstrip().split()
                mxes = mxes.split(",")
            except ValueError:
                event = line.strip()
                mxes = []

            clusters[event] = mxes
            if filter_list is not None:
                if event in filter_list:
                    clusters[event] = mxes
            else:
                clusters[event] = mxes
    return clusters


def getEventCounts(filename, filter_list=None):
    events = []
    counts = []
    with open(filename) as tsv:
        samples = tsv.readline().rstrip().split("\t")[1:]
        for line in tsv:
            row = line.rstrip().split("\t")
            if filter_list is not None:
                if row[0] in filter_list:
                    events.append(row[0])
                    counts.append(row[1:])
            else:
                events.append(row[0])
                counts.append(row[1:])
    counts = np.array(counts, dtype=float)
    return samples, events, counts


########################################################################
# MAINE
#
#
########################################################################


def add_parser(parser):
    parser.add_argument(
        "--inclusionSPLICEDICE",
        type=str,
        required=True,
        help="Compressed NPZ formatted Inclusion count matrix from quantSPLICEDICE.",
    )

    parser.add_argument(
        "-c", "--clusters", type=str, required=True, help="Clusters table."
    )

    parser.add_argument(
        "--chi2",
        action="store_true",
        default=False,
        help="Use X^2 instead of fishers. Quicker, not as sensitive.",
    )

    parser.add_argument(
        "--multiple_test_correction",
        default="pairwise",
        choices=["pairwise", "all", "none"],
        help="Correction via Benjamini-Hochberg. Options are "
        "[pairwise (default): all p-values per sample pair; "
        "all: all p-values for every sample pair together;\n"
        "none: no multiple test correction]",
    )

    parser.add_argument(
        "-f",
        "--filter_list",
        help="txt file where each line is a event to analyze, can speed up fisher test analysis",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="pairwise.tsv",
        help="tab-separated output filename (default uses input prefix)",
    )


def run_with(args):
    """ """
    import sys
    import numpy as np
    from scipy.stats import fisher_exact
    from scipy.stats import chi2_contingency
    from statsmodels.stats.multitest import multipletests

    # Input file arguments
    inclusionCounts = args.inclusionSPLICEDICE
    allClusters = args.clusters
    outfilename = args.output

    if args.filter_list is not None:
        with open(args.filter_list, "r") as filter_list_file:
            filter_list = set([line.rstrip() for line in filter_list_file])
    else:
        filter_list = None

    if args.chi2:
        test_method = chi2_contingency
    else:
        test_method = fisher_exact

    samples, events, counts = getEventCounts(inclusionCounts, filter_list)
    print("Counts loaded from", inclusionCounts, "...")
    clusters = getClusters(allClusters, filter_list)
    print("Clusters loaded from", allClusters, "...")
    pairs = []
    for i in range(len(samples) - 1):
        for j in range(i + 1, len(samples)):
            pairs.append((i, j))

    columns = [f"{samples[a]}_{samples[b]}" for a, b in pairs]
    print("Analyzing pairs:")
    print(",".join(columns))

    totaln = len(events)
    parray = []

    for n, inclusions in enumerate(counts):

        if n % 50 == 0:
            print(f"[{n} / {totaln}] events analyzed...")
        mxes = counts[np.isin(events, clusters[events[n]])]

        exclusions = np.sum(mxes, axis=0)

        pvalues = []

        for a, b in pairs:
            table = [[inclusions[a], inclusions[b]], [exclusions[a], exclusions[b]]]

            # chi^2 test will fail on contingency tables where any column or
            # rows are zero, but not if diagonal is zero
            # [[0, 0], [2, 3]] fail
            # [[0, 2], [0, 3]] fail
            # [[0, 2], [2, 0]] still valid

            # if args.chi2:
            #    col_total = [inclusions[a] + inclusions[b], exclusions[a] + exclusions[b]]
            #    row_total = [inclusions[a] + exclusions[a], inclusions[b] + exclusions[b]]
            #    if any(count == 0 for count in row_total + col_total):
            #        continue

            pvalues.append(test_method(table)[1])
        parray.append(pvalues)

    if args.multiple_test_correction == "all":
        parray = np.array(parray)
        shape = parray.shape
        corrected = multipletests(parray.flatten(), method="fdr_bh")[1]
        parray = np.array(corrected).reshape(shape)
    elif args.multiple_test_correction == "pairwise":
        parray = np.array(parray)
        for i in range(parray.shape[1]):
            corrected = multipletests(parray[:, i], method="fdr_bh")[1]
            parray[:, i] = corrected
    elif args.multiple_test_correction == "none":
        pass

    with open(outfilename, "w") as outfile:
        tab = "\t"
        columns = tab.join(f"{samples[a]}_{samples[b]}" for a, b in pairs)
        outfile.write(f"clusterID\t{columns}\n")
        for n, pvalues in enumerate(parray):
            outfile.write(f"{events[n]}\t{tab.join(str(p) for p in pvalues)}\n")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    add_parser(parser)
    args = parser.parse_args()
    run_with(args)
