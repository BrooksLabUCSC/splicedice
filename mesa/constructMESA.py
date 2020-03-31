#!/usr/bin/env python3


########################################################################
# File: constructMESA.py
#  executable: constructMESA.py
# Purpose:
#
#
# Author: Cameron M. Soulette
# History:      cms 06/20/2018 Created
#
########################################################################

# IMPORT
import os
import sys
from tqdm import tqdm
from multiprocessing import Pool
import tempfile
import pybedtools


########################################################################
# Helper Functions
#
#
########################################################################


def bedToCoordSet(f):
    """
    takes bed file and returns a set of junctions
    """
    sampName, fname, tdir = f
    newOut = fname.split("/")[-1]
    outName = os.path.join(tdir, "%s.tmp" % sampName)
    jSet = set()
    with open(outName, "w") as fout:
        with open(fname, "r") as lines:
            for line in lines:
                cols = line.rstrip().split()
                chrom, c1, c2, i, score, strand = cols[:6]
                c1, c2, score = map(int, [c1, c2, score])
                if score < readThresh or c2 - c1 > lengthThresh:
                    continue
                cols[4] = "0"

                print("\t".join(cols[:6]), file=fout)


def runCMD(x):
    """
    mutiprocessing helper function to run
    multiple instances of bedToCoordSet
    """
    return bedToCoordSet(x)


def makeClusters(intersection):
    """
    takes a bedtools intersection with -wa -wb options and returns
    a list of overlaps between -a and all -b

    e.g. of intersection

    GL000008.2  164884  170271  .   0   -   GL000008.2  164884  170271  .   0   -
    """
    mxeDict = dict()
    for i in intersection:
        left = "%s:%s-%s" % (i[0], i[1], i[2])
        right = "%s:%s-%s" % (i[6], i[7], i[8])
        if left == right:
            continue
        else:
            if left not in mxeDict:
                mxeDict[left] = list()
            mxeDict[left].append(right)

    return mxeDict


def checkFileExist(f):
    if os.path.isfile(f):
        pass
    else:
        print("** ERR: %s file does not exist. **" % f, file=sys.stderr)
        sys.exit(1)


def juncFilesFromManifest(manifest, tdir):
    """function takes in 2 column manifest
    and returns list of file names (second column)"""

    beds = list()
    with open(manifest, "r") as lines:
        for line in lines:
            data = line.rstrip().split("\t")
            beds.append((data[0], data[1], tdir))
    return beds


def filterJunctions(bed, filta, correct, genome):
    """
    reads in junction bed file, and filters based on ss motif.
    also resolves if true
    """
    finalBed = str()
    btJuncs = pybedtools.BedTool(bed)
    dinucSeq = btJuncs.sequence(fi=genome, s=True, tab=True)
    if correct:
        print("** Under construction **", file=sys.stderr)
        sys.exit(1)
        with open(dinucSeq.seqfn) as fileObj:
            for i in fileObj:
                head, seq = i.split()
                donor, acceptor = seq[:2].upper(), seq[-2:].upper()
    else:
        dinucSeq = btJuncs.sequence(fi=genome, s=True, tab=True)
        acceptableDinucs = set(["GTAG"])
        if filta == "gc_at":
            acceptableDinucs.add(["GTAG"])
            acceptableDinucs.add(["ATAC"])

        with open(dinucSeq.seqfn) as fileObj:
            for i in fileObj:
                head, seq = i.split()
                l, strand = head.split("(")
                seqID = l.replace("-", ":")
                seqID = seqID.split(":")
                if filta == "all":
                    finalBed += "%s\t%s\t%s\t%s\t%s\t%s\n" % (
                        seqID[0],
                        seqID[1],
                        seqID[2],
                        head,
                        0,
                        strand.rstrip(")"),
                    )
                else:
                    donor, acceptor = seq[:2].upper(), seq[-2:].upper()
                    if donor + acceptor in acceptableDinucs:
                        finalBed += "%s\t%s\t%s\t%s\t%s\t%s\n" % (
                            seqID[0],
                            seqID[1],
                            seqID[2],
                            head,
                            0,
                            strand.rstrip(")"),
                        )
    return finalBed


########################################################################
# Main
# Here is the main program
#
########################################################################


def run_with(args):
    global minLengthThresh
    global lengthThresh
    global readThresh

    # filtering
    minLengthThresh = args.min_length
    lengthThresh = args.max_length
    readThresh = args.support_threshold
    threads = args.threads

    # required args
    fPrefix = args.output_prefix
    bedList = args.manifest
    genome = args.genome
    ssFilter = args.filter
    correct = args.resolve_strands

    # initialize tempdir
    with tempfile.TemporaryDirectory() as tmpdirname:
        print(tmpdirname)
        # Collect junction bed files
        beds = juncFilesFromManifest(bedList, tmpdirname)
        # start worker pool
        with Pool(threads) as p:
            for fname in tqdm(
                p.imap_unordered(runCMD, beds),
                total=len(beds[:]),
                position=1,
                desc="1/2 Filtering sample junctions",
            ):
                continue

        # get list of temp files
        flist = [os.path.join(tmpdirname, x) for x in os.listdir(tmpdirname)]
        junctions = set()
        with open("%s_junctions.bed" % fPrefix, "w") as fout:
            for fname in tqdm(
                flist, total=len(flist), position=1, desc="2/2 Concatenating junctions"
            ):
                with open(fname, "rb") as fin:
                    for i in fin:
                        if i in junctions:
                            continue
                        else:
                            junctions.add(i)
                            print(i.rstrip().decode(), file=fout)
    # Filter junctions
    filteredBed = filterJunctions(
        "%s_junctions.bed" % fPrefix, ssFilter, correct, genome
    )
    bTool = pybedtools.BedTool(filteredBed, from_string=True)
    bTool.sort()
    intersection = bTool.intersect(bTool, s=True, wa=True, wb=True)
    bTool.saveas("%s_junctions.bed" % fPrefix)

    clusters = makeClusters(intersection)

    with open("%s_all_clusters2.tsv" % fPrefix, "w") as out:
        for i, c in clusters.items():
            print(i, ",".join(c), file=out)
    print("done.", file=sys.stderr)
