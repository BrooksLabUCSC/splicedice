#!/usr/bin/env python3


########################################################################
# File: compareSampleSets.py
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

import sys
import numpy as np
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests

########################################################################
# CommandLine
########################################################################


class CommandLine(object):
    """
    Handle the command line, usage and help requests.  CommandLine uses
    argparse, now standard in 2.7 and beyond.  it implements a standard command
    line argument parser with various argument options, and a standard usage
    and help,
    attributes:
    myCommandLine.args is a dictionary which includes each of the available
    command line arguments as myCommandLine.args['option']

    methods:

    """

    def __init__(self, inOpts=None):
        """
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using
        argparse.
        """
        import argparse

        self.parser = argparse.ArgumentParser(
            description="TBD",
            epilog="Please feel free to forward any usage questions or concerns",
            add_help=True,  # default is True
            prefix_chars="-",
            usage="%(prog)s -m1 manifest1.txt -m2 manifest2.txt",
        )
        # Add args
        self.parser.add_argument(
            "--psiMESA",
            type=str,
            action="store",
            required=True,
            help="Compressed NPZ formatted PSI matrix from quantMESA.",
        )
        self.parser.add_argument(
            "-m1",
            "--manifest1",
            type=str,
            action="store",
            required=True,
            help="Manifest containing samples for sample set group1",
        )
        self.parser.add_argument(
            "-m2",
            "--manifest2",
            type=str,
            action="store",
            required=True,
            help="Manifest containing samples for sample set group2",
        )
        self.parser.add_argument(
            "-o",
            "--out_prefix",
            type=str,
            action="store",
            required=False,
            help="Prefix for output file.",
        )

        if inOpts is None:
            self.args = vars(self.parser.parse_args())
        else:
            self.args = vars(self.parser.parse_args(inOpts))


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


########################################################################
# MAINE
#
#
########################################################################

def add_parser(subparser):
    css_parser = subparser.add_parser("compare_sample_sets")
    css_parser.add_argument(
        "--psiMESA",
        type=str,
        required=True,
        help="Compressed NPZ formatted PSI matrix from 'mesa quant'.",
    )
    css_parser.add_argument(
        "-m1",
        "--manifest1",
        type=str,
        required=True,
        help="Manifest containing samples for sample set group1",
    )
    css_parser.add_argument(
        "-m2",
        "--manifest2",
        type=str,
        required=True,
        help="Manifest containing samples for sample set group2",
    )
    css_parser.set_defaults(func=run_with)


def run_with(args):
    """
    A workflow to compute the significance difference
    between two distributions of PSI values.
    Values are assumed to not be normall distributed, thus
    we invoke the wilcoxon ranksum test as the statistical analysis.
    """
    pmesa = args.psiMESA
    group1 = args.manifest1
    group2 = args.manifest2

    # get sample lists
    g1 = returnSamplesFromManifest(group1)
    g2 = returnSamplesFromManifest(group2)

    if len(g1) < 3 or len(g2) < 3:
        print(
            "Cannot conduct wilcoxon with less than 3 samples in either group. Exit.",
            file=sys.stderr,
        )
        sys.exit(1)

    # load psi
    data = loadNPZ(pmesa)

    # table has 3 arrays, cols, rows and data
    cols, rows, matrix = data["cols"], data["rows"], data["data"]

    # get sample indices
    g1Indices = getColIndexFromArray(g1, cols)
    g2Indices = getColIndexFromArray(g2, cols)

    # do the math
    pvals = list()
    testedEvents = list()

    for n, event in enumerate(matrix):
        d1, d2 = event[g1Indices], event[g2Indices]
        nonans1 = np.invert(np.isnan(d1))
        nonans2 = np.invert(np.isnan(d2))
        data1 = d1[nonans1]
        data2 = d2[nonans2]

        if len(data1) < 3 or len(data2) < 3:
            continue

        D, pval = ranksums(d1, d2)
        testedEvents.append((rows[n], np.mean(data1) - np.mean(data2)))
        pvals.append(pval)

    # correct pvals
    corrected = multipletests(pvals, method="fdr_bh")[1]
    for n, i in enumerate(testedEvents):
        print(pvals[n], corrected[n], i[0], i[1])
