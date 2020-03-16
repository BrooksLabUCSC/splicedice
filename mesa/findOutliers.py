#!/usr/bin/env python3


########################################################################
# File: findOutliers.py
#  executable:
# Purpose:
#
#
# Author: Cameron M. Soulette
# History:      cms 03/04/2020 Created
#
########################################################################

########################################################################
# Hot Imports & Global Variable
########################################################################


import sys
import numpy as np


########################################################################
# CommandLine
########################################################################

def add_parser(subparser):
    """Add command line argument parser for findOutliers module"""
    outlier_parser = subparser.add_parser("findOutliers")
    outlier_parser.add_argument(
        "--psiMESA",
        type=str,
        action="store",
        required=True,
        help="Compressed NPZ formatted PSI matrix from quantMESA.",
    )

    outlier_parser.add_argument(
        "-m",
        "--manifest",
        action="store",
        required=True,
        help="Sample manifest for samples you want ",
    )
    outlier_parser.add_argument(
        "--nullMan",
        action="store",
        required=False,
        default=None,
        help="List of samples you want to use to use as baseline/background. "
             "Samples must be in table.",
    )
    outlier_parser.add_argument(
        "--outlierCutoff",
        type=int,
        action="store",
        required=False,
        default=3,
        help="Report events where zScore >= N (default 3).",
    )
    outlier_parser.add_argument(
        "--dpsiThrsh",
        type=float,
        action="store",
        required=False,
        default=0.1,
        help="Report events where dPSI >= N (default 0.1).",
    )
    outlier_parser.set_defaults(func=run_with)


########################################################################
# Funktions
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


########################################################################
# Main
########################################################################

def run_with(args):
    pmesa = args.psiMESA
    man = args.manifest
    null = args.nullMan
    zthresh = args.outlierCutoff
    pthresh = args.outlierCutoff

    # get samples
    with open(man) as fin:
        samples = [x.split("\t")[0] for x in fin]

    if null is None:
        nullSamps = samples
    else:
        with open(null) as fin:
            nullSamps = [x.split("\t")[0] for x in fin]

    if len(nullSamps) < 10:
        print(
            "Less than 10 samples detected in null group. Too few. Exit.",
            file=sys.stderr,
        )
        sys.exit(1)
    # load psi
    data = loadNPZ(pmesa)
    mSamps, mEvents, matrix = data["cols"], data["rows"], data["data"]
    # m is for matrix
    nullIndices = np.argwhere(np.isin(mSamps, nullSamps))[:, 0]
    outlierIndices = np.argwhere(np.isin(mSamps, samples))[:, 0]

    # ok. it is more intuitive to go line by line.. should I?....

    for ePos, i in enumerate(matrix):
        # yea.
        nullVals = i[nullIndices]
        if np.sum(np.isnan(nullVals)) / len(nullVals) > 0.2:
            continue
        std = np.nanstd(nullVals)

        if std < 0.001:
            continue

        mean = np.nanmean(nullVals)

        vals = i[outlierIndices]
        zScores = (vals - mean) / std

        [
            print(
                mEvents[ePos],
                mSamps[outlierIndices[pos]],
                x,
                vals[pos],
                mean,
                std,
                sep="\t",
            )
            for pos, x in enumerate(zScores)
            if x >= zthresh and abs(x - mean) >= pthresh
        ]
