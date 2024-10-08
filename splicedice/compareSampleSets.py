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
# Helper Functions
#
#
########################################################################


def getAnnotated(annotation):
    gid_coords = {}
    gene_coords = {}
    genes = {}
    transcripts = {}

    with open(annotation) as gtf:
        for line in gtf:

            if line.startswith("#"):
                continue

            row = line.rstrip().split('\t')
            info = [x.split('"') for x in row[8].split(';')]
            chrom = row[0]
            strand = row[6]

            start = int(row[3])
            stop = int(row[4])-1

            if row[2] == "transcript":
                tid =  [x[1] for x in info if 'transcript_id' in x[0]][0]
                gene_name = [x[1] for x in info if 'gene_name' in x[0]][0]
                genes[tid] = gene_name
                transcripts[(tid,chrom,strand)] = []

            elif row[2] == "exon":
                tid =  [x[1] for x in info if 'transcript_id' in x[0]][0]
                transcripts[(tid,chrom,strand)].append((start,stop))

            elif row[2] == "gene":
                gene_name = [x[1] for x in info if 'gene_name' in x[0]][0]
                gid = [x[1] for x in info if 'gene_id' in x[0]][0]
                try:
                    try:
                        gene_coords[(chrom,strand)] [(start,stop)].append(gene_name)
                    except KeyError:    
                        gene_coords[(chrom,strand)] [(start,stop)] = [gene_name]

                    try:
                        gid_coords[(chrom,strand)] [(start,stop) ].append(gid)
                    except KeyError:
                        gid_coords[(chrom,strand)] [(start,stop) ] = [gid]
                except KeyError:
                    gene_coords[(chrom,strand)] = {(start,stop) : [gene_name]}
                    gid_coords[(chrom,strand)] = {(start,stop) : [gid]}              
                
    annotated = {}
    transcript_ids = {}
    for transcript,exons in transcripts.items():
        tid,chromosome,strand = transcript
        for i in range(len(exons)-1):
            junction = (chromosome,exons[i][1],exons[i+1][0],strand)
            if junction in annotated:
                if genes[tid] not in annotated[junction]:
                    annotated[junction].append(genes[tid])
                    transcript_ids[junction].append(tid)
            else:
                annotated[junction] = [genes[tid]]
                transcript_ids[junction] = [tid]

    return annotated,gene_coords,transcript_ids


def getColIndexFromArray(x, y):
    """
    takes in list of strings = x
    and finds list index in array = y
    """

    return np.nonzero(np.isin(y, x))


def returnSamplesFromManifest(x):
    """
    reads in splicedice formatted manifest
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

def add_parser(parser):
    parser.add_argument(
        "--psiSPLICEDICE",
        type=str,
        required=True,
        help="Compressed NPZ formatted PSI matrix from 'splicedice quant'.",
    )
    parser.add_argument(
        "-m1",
        "--manifest1",
        type=str,
        required=True,
        help="Manifest containing samples for sample set group1",
    )
    parser.add_argument(
        "-m2",
        "--manifest2",
        type=str,
        required=True,
        help="Manifest containing samples for sample set group2",
    )
    parser.add_argument(
        "-a",
        "--annotation",
        type=str,
        required=False,
        default="",
        help="Optional GTF file to label known splice junctions and genes",
    )
    parser.add_argument(
        "-o",
        "--outputFile",
        type=str,
        required=True,
        help="Output filename for tab-separated table",
    )



def run_with(args):
    """
    A workflow to compute the significance difference
    between two distributions of PSI values.
    Values are assumed to not be normall distributed, thus
    we invoke the wilcoxon ranksum test as the statistical analysis.
    """
    psplicedice = args.psiSPLICEDICE
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
    #data = loadNPZ(psplicedice)

    # table has 3 arrays, cols, rows and data
    #cols, rows, matrix = data["cols"], data["rows"], data["data"]

    
    ######
    rows = []
    data = []
    with open(psplicedice) as tsv:
        headers = tsv.readline().strip().split('\t')[1:]
        for line in tsv:
            row = line.strip().split('\t')
            rows.append(row[0])
            data.append(row[1:])
            
    matrix = np.array(data,dtype="float32")
    rows = np.array(rows)
    cols = np.array(headers)
    

    
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

        D, pval = ranksums(data1, data2)
        med1 = np.median(data1)
        med2 = np.median(data2)
        mean1 = np.mean(data1)
        mean2 = np.mean(data2)
        testedEvents.append((rows[n], med1, med2, mean1, mean2, med1-med2))
        pvals.append(pval)

    # correct pvals
    corrected = multipletests(pvals, method="fdr_bh")[1]
    
    with open(args.outputFile,"w") as tsv:
        if args.annotation:
            print("event\tmean1\tmean2\tmedian1\tmedian2\tdelta\tp-value\tcorrected\tgene\toverlapping\ttranscript_id",file=tsv)

            annotated,gene_coords,transcript_ids = getAnnotated(args.annotation)
            for n, event in enumerate(testedEvents):
                name,med1,med2,mean1,mean2,delta = event

                chromosome,coords,strand = name.split(":")
                start,stop = [int(x) for x in coords.split("-")]

                start -= 1
                stop += 1

                junction = (chromosome,start,stop,strand)

                overlaps = []
                try:
                    for gene_start, gene_stop in gene_coords[(chromosome,strand)].keys():
                        if (start >= gene_start and start <= gene_stop) or (stop >= gene_start and stop <= gene_stop):
                            overlaps.extend(gene_coords[(chromosome,strand)][gene_start,gene_stop])
                except KeyError:
                    pass

                print(name,mean1,mean2,med1,med2,delta, pvals[n], corrected[n], 
                      ','.join(annotated.get(junction,["nan"])), ','.join(overlaps),
                      ','.join(transcript_ids.get(junction,["nan"])),
                      sep='\t',file=tsv)
    
        else:
            print("event\tmean1\tmean2\tmedian1\tmedian2\tdelta\tp-value\tcorrected",file=tsv)
            for n,event in enumerate(testedEvents):
                name,med1,med2,mean1,mean2,delta = event
                print(name,mean1,mean2,med1,med2,delta, pvals[n], corrected[n], sep='\t',file=tsv)
