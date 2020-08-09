#!/usr/bin/env python3
"""This script is used to convert SJ.out data files outputted from the STAR
aligner into bed files for further processing."""


def add_parser(parser):
    parser.add_argument("-s", "--star-tab", required=True)
    parser.add_argument("--min-overhang", type=int, default=5)
    parser.add_argument("--min-unique", type=int, default=5)
    parser.add_argument("--no-multimap", default=False, action="store_true")


def run_with(args):
  
    """Parse every line from the star_tab file, and converts to a bed file,
    with the 5th column of the bed file holding either the number of uniquely
    mapping reads to that junction or the number of unique and multimapping
    reads."""

    with open(f"{args.star_tab}.bed", "w") as bed_file:
        with open(args.star_tab, "r") as star_file:
            for line in star_file:
                cols = line.rstrip().split()
                chrom, start, end, strand = cols[:4]
                start, end, strand = int(start), int(end), int(strand)
                unique = int(cols[6])
                multi = int(cols[7])
                overhang = int(cols[-1])
                
                # Convert strand number to symbol
                if strand == 1:
                    strand = "+"
                elif strand == 2:
                    strand = "-"
               
                if (
                    overhang < args.min_overhang
                    or unique < args.min_unique
                    or strand == 0
                ):
                    continue

                if args.no_multimap:
                    print(
                        chrom,
                        start - 1,
                        end,
                        ".",
                        unique,
                        strand,
                        sep="\t",
                        file=bed_file,
                    )
                else:
                    print(
                        chrom,
                        start - 1,
                        end,
                        ".",
                        unique + multi,
                        strand,
                        sep="\t",
                        file=bed_file,
                    )


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    add_parser(parser)
    args = parser.parse_args()
    run_with(args)
