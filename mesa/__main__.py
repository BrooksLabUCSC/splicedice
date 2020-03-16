#!/usr/bin/env python3
"""This script is the main entry point for the different MESA commands."""
import argparse

from . import quantMESA
from . import star_junc_to_bed
from . import pairwise_fisher
from . import compareSampleSets
from . import clusterMESA
from . import findOutliers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers()

    star_junc_to_bed.add_parser(subparsers)
    quantMESA.add_parser(subparsers)
    clusterMESA.add_parser(subparsers)
    compareSampleSets.add_parser(subparsers)
    pairwise_fisher.add_parser(subparsers)
    findOutliers.add_parser(subparsers)

    args = parser.parse_args()

    if "func" in args.__dict__:
        args.func(args)
    else:
        parser.print_usage()


if __name__ == "__main__":
    main()
