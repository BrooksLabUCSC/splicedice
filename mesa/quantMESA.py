#!/usr/bin/env python3

########################################################################
# File: quantMESA.py
#  executable: quantMESA.py
# Purpose:
#
#
# Author: Cameron M. Soulette
# History:      cms 06/20/2018 Created
#
########################################################################

########################################################################
# Hot Imports & Global Variable
########################################################################

import sys
import numpy as np
from multiprocessing import Pool
from tqdm import tqdm

from . import constructMESA

########################################################################
# CommandLine
########################################################################


def is_valid_filter(parser, arg):
    if arg not in ["gtag_only", "gc_at", "all"]:
        parser.error("Filter argument %s is invalid!" % arg)
        sys.exit(1)
    else:
        return arg


def add_parser(parser):
    """Add parsing options for quant module of MESA"""
    parser.add_argument(
        "-m",
        "--manifest",
        required=True,
        help="List of bed files containing stranded junctions.",
    )
    parser.add_argument("-g", "--genome", required=True)
    parser.add_argument("-o", "--output-prefix", required=True)

    # TODO can use choices parameter to avoid invalidation
    parser.add_argument(
        "--filter",
        default="gtag_only",
        type=lambda x: is_valid_filter(parser, x),
        help="Filter splice sites that\n"
        'filter non GT-AG splice sites : "gtag_only" (default),\n'
        'filter non GT/GC-AG and AT-AC introns: "gc_at",\n'
        'keep all splice sites: "all"',
    )
    parser.add_argument(
        "--drim-table", default=False, action="store_true", help="Generate a table to use with DRIM-Seq."
    )
    parser.add_argument(
        "--resolve-strands",
        action="store_true",
        default=False,
        help="Resolve splice site strand based on dinucleotide.",
    )
    parser.add_argument(
        "--support-threshold",
        default=10,
        type=int,
        help="Filter junctions < non-normalized read counts (10)",
    )
    parser.add_argument(
        "--max-length",
        default=50000,
        type=int,
        help="Filter junctions > N length (50,000)",
    )
    parser.add_argument(
        "--min-length", default=50, type=int, help="Filter junctions < N length (50)"
    )
    parser.add_argument(
        "--min-samp-threshold",
        default=0.01,
        type=float,
        help="Filter junctions found in < N samples (0.1)",
    )

    parser.add_argument(
        "--uncompressed-only",
        action="store_false",
        required=False,
        default=True,
        help="NPZ compression",
    )
    parser.add_argument(
        "-p", "--threads", type=int, default=2, help="Number of threads."
    )
    parser.add_argument(
        "--inclusion-threshold",
        type=int,
        default=10,
        help="Filter events with less than N reads (default 10)",
    )
    # TODO: Maybe better named just as --read-threshold
    parser.add_argument(
        "--psi-threshold",
        type=int,
        default=20,
        help="Do not compute PSI for samples with less than N reads (default 20)",
    )


########################################################################
# Clusters & Junctions
########################################################################


class Junction(object):
    def __init__(self, chrom=None, c1=None, c2=None, strand=None):
        self.chrom = chrom
        self.c1 = c1
        self.c2 = c2
        self.strand = strand
        self.mxes = list()
        self.name = "%s:%s-%s" % (chrom, c1, c2)
        self.quant = list()


class Sample(object):
    def __init__(self, name=None):
        self.sid = name
        self.introns = dict()


########################################################################
# Helper Functions
#
#
########################################################################


def runCMD2(x):
    """
    mutiprocessing helper function to run
    multiple instances of bedToCoordSet
    """
    return clusterQuant(x)


def clusterQuant(data):
    """
    stuff
    """
    root, fname, order = data

    order = np.load(order)["introns"]
    d = {x: np.nan for x in order}
    with open(fname) as fin:
        for i in fin:

            cols = i.rstrip().split()
            intronID = "%s:%s-%s" % (cols[0], cols[1], cols[2])
            if intronID in d:
                d[intronID] = int(cols[4])
    return [root, np.array([d[x] for x in order], dtype=np.float32)]


########################################################################
# Main
# Here is the main program
#
########################################################################
def run_with(args):
    constructMESA.run_with(args)

    bedList = args.manifest
    outPrefix = args.output_prefix
    junctionBed = "{}_junctions.bed".format(outPrefix)
    clusters = "{}_all_clusters2.tsv".format(outPrefix)
    comp = args.uncompressed_only
    threads = args.threads
    incT = args.inclusion_threshold
    psiT = args.psi_threshold
    drim = args.drim_table

    beds = list()
    with open(bedList, "r") as fnames:
        for num, fdata in enumerate(fnames, 0):
            root, fname, group1, group2 = fdata.rstrip().split("\t")
            group1, group2 = group1.replace(" ", ""), group2.replace(" ", "")
            beds.append((root, fname, "order.npz"))

    # Make junction objects
    n = 0
    order = list()
    orderDict = dict()
    with open(junctionBed, "r") as lines:
        for line in tqdm(lines, desc="Initializing junction counts", position=1):
            j = line.rstrip().split()
            c1, c2 = map(int, [j[1], j[2]])
            jid = "%s:%s-%s" % (j[0], j[1], j[2])
            order.append(jid)
            orderDict[jid] = n
            n += 1

    # this compressed order array will be used to re-order junction counts
    # from individual bed files
    order = np.array(order)
    np.savez_compressed("order.npz", introns=order)
    allData = list()
    samps = list()
    with Pool(threads) as p:
        for i in tqdm(
            p.imap_unordered(runCMD2, beds),
            total=len(beds),
            desc="Computing counts",
            position=1,
        ):
            allData.append(i[-1])
            samps.append(i[0])

    # print("current mem %s" % memory_usage_psutil())
    allData = np.nan_to_num(np.transpose(np.array(allData)))
    # print(allData.shape)

    dataQ = list()
    dataP = list()
    events = list()
    incOut = open("%s_inclusionCounts.tsv" % outPrefix, "w")
    psiOut = open("%s_allPSI.tsv" % outPrefix, "w")

    if drim:
        drimOut = open("%s_drimTable.tsv" % outPrefix, "w")

    print("\t".join(samps), "cluster", sep="\t", file=incOut)
    print("\t".join(samps), "cluster", sep="\t", file=psiOut)

    if drim:
        print("gene", "feature_id", "\t".join(samps), sep="\t", file=drimOut)

    cluster_num = 0
    with open(clusters, "r") as lines:
        for line in tqdm(lines, desc="Compute inclusion/exclusion counts", position=1):
            intron, mxes = line.rstrip().split()
            mxes = mxes.split(",")

            intronPos = orderDict[intron]
            mxePos = [orderDict[x] for x in mxes]

            inclusions = allData[intronPos]
            exclusions = allData[mxePos]

            total = exclusions.sum(axis=0) + inclusions

            # lets filter...

            # if np.max(inclusions)<incT:
            #     continue

            events.append(intron)
            psiQuants = inclusions / total
            psiQuants[total < psiT] = np.nan
            dataQ.append(inclusions)
            dataP.append(psiQuants)

            print(
                "\t".join("%.2f" % x for x in psiQuants), intron, sep="\t", file=psiOut
            )
            print("\t".join(str(x) for x in inclusions), intron, sep="\t", file=incOut)

            if drim:
                print(
                    "cl_%s_%s\t%s_%s" % (cluster_num, intron, intron, cluster_num),
                    "\t".join(str(x) for x in inclusions),
                    sep="\t",
                    file=drimOut,
                )
                # print other intron values too
                for mxeNum, x in enumerate(exclusions):
                    print(
                        "cl_%s_%s\t%s_%s"
                        % (cluster_num, intron, mxes[mxeNum], cluster_num),
                        "\t".join(str(y) for y in x),
                        sep="\t",
                        file=drimOut,
                    )

            cluster_num += 1

    cols, rows, data1 = np.array(samps), np.array(events), np.array(dataQ)
    data2 = np.array(dataP, dtype=np.float32)

    if comp:
        np.savez_compressed(
            "%s_inclusionCounts.npz" % outPrefix, cols=cols, rows=rows, data=data1
        )
        np.savez_compressed(
            "%s_allPSI.npz" % outPrefix, cols=cols, rows=rows, data=data2
        )

    incOut.close()
    psiOut.close()
    if drim:
        drimOut.close()
