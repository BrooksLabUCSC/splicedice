#!/usr/bin/env python3

import numpy as np
import os
    
    
def add_parser(parser):
    parser.add_argument("-i", "--inclusionCounts",
                        action="store",
                        help="")
    parser.add_argument("-c", "--clusters",
                        action="store",
                        help="allClusters.tsv file with mutually exclusive clusters for each junction.")

    parser.add_argument("-d", "--coverageDirectory",
                        action="store",
                        help="")
    parser.add_argument("-o", "--outputPrefix",
                        action="store",
                        help="")
    parser.add_argument("-r", "--makeRSDtable",
                        action="store_true",
                        help="Make a table of relative standard deviations in coverage across intron.")
    parser.add_argument("-s", "--singleJunctionCalculation",
                        action="store_true",
                        help="Calculate IR value using individual junction counts, and not count of all junctions in cluster.")
    parser.add_argument("-a", "--annotation",
                        action="store",
                        help="GTF file with gene annotation.")
    parser.add_argument("-j", "--allJunctions",
                        action="store_true",
                        help="Output IR values for all junctions above RSD threshold. Default: only annotated junctions")
    parser.add_argument("-t", "--RSDthreshold",
                        default=1.0, action="store",
                        help="RSD cutoff for inclusion. Default: 1.0.")
                        
    
    
def getAnnotated(annotation):
    genes = {}
    names = {}
    transcripts = {}
    with open(annotation) as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            row = line.rstrip().split('\t')
            if row[2] == "transcript":
                info = [x.split('"') for x in row[8].split(';')]
                tid = [x[1] for x in info if 'transcript_id' in x[0]][0]
                try:
                    gid = [x[1] for x in info if 'gene_name' in x[0]][0]
                except IndexError:
                    gid = [x[1] for x in info if 'gene_id' in x[0]][0]
                genes[tid] = gid
                transcripts[(tid, row[0], row[6])] = []
            if row[2] == "exon":
                info = [x.split('"') for x in row[8].split(';')]
                tid = [x[1] for x in info if 'transcript_id' in x[0]][0]
                transcripts[(tid, row[0], row[6])].append((int(row[3]), int(row[4])))
    #annotated = {}
    annotated = set()
    for transcript, exons in transcripts.items():
        tid, chromosome, strand = transcript
        for i in range(len(exons)-1):
            #annotated[(chromosome,exons[i][1],exons[i+1][0]-1,strand)] = genes[tid]
            annotated.add(f"{chromosome}:{exons[i][1]+1}-{exons[i+1][0]-1}:{strand}")
    return annotated


    
def getInclusionCounts(filename):
    with open(filename) as inclusionCounts:
        header = inclusionCounts.readline().strip().split("\t")
        counts = {sample:{} for sample in header[1:]}
        for line in inclusionCounts:
            row = line.rstrip().split("\t")
            for i, count in enumerate(row[1:]):
                counts[header[i+1]][row[0]] = float(count)
    return counts


def getClusters(filename):
    clusters = {}
    with open(filename) as allClusters:
        for line in allClusters:
            row = line.strip().split("\t")
            try:
                clusters[row[0]] = row[1].split(",")
            except IndexError:
                clusters[row[0]] = []
    return clusters


def getJunctions(samples, coverageDirectory, annotated, args):
    junctions = set()
    for sample in samples:
        filename = os.path.join(coverageDirectory, f"{sample}_intron_coverage.txt")
        with open(filename) as percentileCoverage:
            for line in percentileCoverage:
                row = line.strip().split("\t")
                cluster = f"{row[0]}:{int(row[1])+1}-{row[2]}:{row[5]}"
                if not args.allJunctions and cluster not in annotated:
                    continue
                junctions.add(cluster)
    return junctions


def getFilteredJunctions(samples, coverageDirectory, annotated, args):
    import time
    t = time.time()
    junctions = getJunctions(samples, coverageDirectory, annotated, args)
    print(f"getJunctions complete: {len(junctions)} junctions. {time.time()-t:.1f}s")
    filtered_junctions = set()
    for sample in samples:
        filename = os.path.join(coverageDirectory, f"{sample}_intron_coverage.txt")
        with open(filename) as percentileCoverage:
            for line in percentileCoverage:
                row = line.strip().split("\t")
                cluster = f"{row[0]}:{int(row[1])+1}-{row[2]}:{row[5]}"
                if cluster not in junctions:
                    continue
                covArray = np.array(row[-1].split(","), dtype=float)
                mean = np.mean(covArray)
                rsd = np.std(covArray) / mean if mean > 0 else np.nan
                if rsd < float(args.RSDthreshold):
                    filtered_junctions.add(cluster)
    print(f"RSD filtering complete: {len(filtered_junctions)} junctions retained. {time.time()-t:.1f}s")
    return filtered_junctions


def calculateIRforSample(sample, coverageDirectory, counts, clusters, junctions, args):
    IR = {}
    RSD = {}
    filename = os.path.join(coverageDirectory, f"{sample}_intron_coverage.txt")
    with open(filename) as percentileCoverage:
        for line in percentileCoverage:
            row = line.strip().split("\t")
            cluster = f"{row[0]}:{int(row[1])+1}-{row[2]}:{row[5]}"
            if cluster not in junctions:
                continue
            median = float(row[4])
            covArray = np.array(row[-1].split(","), dtype=float)
            mean = np.mean(covArray)
            RSD[cluster] = np.std(covArray) / mean if mean > 0 else np.nan
            try:
                intronCount = counts[sample][cluster]
                if not args.singleJunctionCalculation:
                    for mxCluster in clusters[cluster]:
                        try:
                            intronCount += counts[sample][mxCluster]
                        except KeyError:
                            print("mxCluster", sample, cluster, mxCluster)
                try:
                    IR[cluster] = median/(median+intronCount)
                except ZeroDivisionError:
                    IR[cluster] = np.nan
            except KeyError:
                pass
    return IR, RSD


def writeIRtable(samples, coverageDirectory, counts, clusters, outputPrefix, junctions, args):
    tab = "\t"
    total = len(samples)
    with open(f"{outputPrefix}_intron_retention.tsv", "w") as irTable:
        irTable.write(f"Junction\t{tab.join(samples)}\n")
        sampleData = {}
        for i, sample in enumerate(samples):
            IR, RSD = calculateIRforSample(sample, coverageDirectory, counts, clusters, junctions, args)
            sampleData[sample] = IR
            if (i+1) % 50 == 0 or (i+1) == total:
                print(f"IR calculated for {i+1}/{total} samples")
        for junction in sorted(junctions):
            irValues = [f"{sampleData[sample].get(junction, float('nan')):0.03f}" for sample in samples]
            irTable.write(f"{junction}\t{tab.join(irValues)}\n")


def writeRSDtable(samples, coverageDirectory, counts, clusters, outputPrefix, junctions, args):
    tab = "\t"
    total = len(samples)
    with open(f"{outputPrefix}_intron_retention_RSD.tsv", "w") as rsdTable:
        header = tab.join([f'{sample}_RSD' for sample in samples])
        rsdTable.write(f"Junction\t{header}\n")
        sampleData = {}
        for i, sample in enumerate(samples):
            IR, RSD = calculateIRforSample(sample, coverageDirectory, counts, clusters, junctions, args)
            sampleData[sample] = RSD
            if (i+1) % 50 == 0 or (i+1) == total:
                print(f"RSD calculated for {i+1}/{total} samples")
        for junction in sorted(junctions):
            rsd = [f"{sampleData[sample].get(junction, float('nan')):0.03f}" for sample in samples]
            rsdTable.write(f"{junction}\t{tab.join(rsd)}\n")


def run_with(args):
    """ """
    import time
    start = time.time()

    countFile = args.inclusionCounts
    clusterFilename = args.clusters
    coverageDirectory = args.coverageDirectory
    outputPrefix = args.outputPrefix
    annotation = args.annotation

    samples = [s.replace("_intron_coverage.txt", "") for s in os.listdir(coverageDirectory) if s.endswith("intron_coverage.txt")]
    print(f"Starting ir_table with {len(samples)} samples")
    print("Gathering inclusion counts and clusters...")
    counts = getInclusionCounts(countFile)
    clusters = None
    if not args.singleJunctionCalculation:
        clusters = getClusters(clusterFilename)
    print(f"Loaded {len(counts)} samples and {len(clusters) if clusters else 0} clusters. {time.time()-start:.1f}s")

    if not args.allJunctions:
        annotated = getAnnotated(annotation)
        print(f"Annotation loaded: {len(annotated)} annotated junctions. {time.time()-start:.1f}s")
    else:
        annotated = None

    print("Collecting junctions across all samples...")
    junctions = getFilteredJunctions(samples, coverageDirectory, annotated, args)
    print(f"Junction collection and RSD filtering complete: {len(junctions)} junctions retained. {time.time()-start:.1f}s")

    print("Writing IR table...")
    writeIRtable(samples, coverageDirectory, counts, clusters, outputPrefix, junctions, args)
    print(f"IR table written. {time.time()-start:.1f}s")

    if args.makeRSDtable:
        print("Writing RSD table...")
        writeRSDtable(samples, coverageDirectory, counts, clusters, outputPrefix, junctions, args)
        print(f"RSD table written. {time.time()-start:.1f}s")

    print(f"Done. Total runtime: {time.time()-start:.1f}s")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    add_parser(parser)
    args = parser.parse_args()
    run_with(args)

