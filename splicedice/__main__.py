#!/usr/bin/env python3
"""This script is the main entry point for the different SPLICEDICE commands."""
import argparse
from . import SPLICEDICE as quant
from . import intron_coverage as ic
from . import ir_table as it

def add_cmd(name, arg_parser_fn, run_with, subparser):
    command_parser = subparser.add_parser(name)
    arg_parser_fn(command_parser)
    command_parser.set_defaults(main=run_with)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers()
    add_cmd("quant", quant.add_parser, quant.run_with, subparsers)
    add_cmd("intron_coverage", ic.add_parser, ic.run_with, subparsers)
    add_cmd("ir_table", it.add_parser, it.run_with, subparsers)
    args = parser.parse_args()
    if "main" in args.__dict__:
        args.main(args)
    else:
        parser.print_usage()

if __name__ == "__main__":
    main()
