#!/usr/bin/env python3
"""This script is the main entry point for the different MESA commands."""
import argparse

from . import MESA as quant
from . import bam_to_junc_bed as bjb
from . import pairwise_fisher as pf
from . import compareSampleSets as css
from . import clusterMESA as cm
from . import findOutliers as fo
from . import intron_coverage as ic
from . import ir_table as it
from . import subset as subset


def add_cmd(name, arg_parser_fn, run_with, subparser):
    """Adds a subparser command to the parser.  name is a string that denotes
    name of the command being added. This function helps a) simplify adding new
    commands to the mesa tool b) also allows for using each python file as its
    own script instead of relying on the __main__.py parser solely.
    ie allows for `mesa quant` and `python quantMESA.py` to be easier to use.

    arg_parser_fn is a function that takes a subparser and adds arguments based
    on the command. Modules typically expose an add_parser() function that does
    this so you can pass that directly to this function.

    subparser is the object the add_subparsers() method that adds the subparser
    to the ArgumentParser
    """
    command_parser = subparser.add_parser(name)
    arg_parser_fn(command_parser)
    command_parser.set_defaults(main=run_with)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers()

    add_cmd("star_junc_to_bed", sjb.add_parser, sjb.run_with, subparsers)
    add_cmd("bam_to_junc_bed", bjb.add_parser, bjb.run_with, subparsers)
    add_cmd("quant", quant.add_parser, quant.run_with, subparsers)
    add_cmd("intron_coverage", ic.add_parser, ic.run_with, subparsers)
    add_cmd("ir_table", it.add_parser, it.run_with, subparsers)
    add_cmd("cluster", cm.add_parser, cm.run_with, subparsers)
    add_cmd("compare_sample_sets", css.add_parser, css.run_with, subparsers)
    add_cmd("pairwise", pf.add_parser, pf.run_with, subparsers)
    add_cmd("findOutliers", fo.add_parser, fo.run_with, subparsers)
    add_cmd("subset", subset.add_parser, subset.run_with, subparsers)

    args = parser.parse_args()

    if "main" in args.__dict__:
        args.main(args)
    else:
        parser.print_usage()


if __name__ == "__main__":
    main()
