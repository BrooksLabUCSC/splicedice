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
    parser.add_argument("-n", "--numThreads",
                        default=1, type=int, action="store",
                        help="Number of parallel workers for junction collection. Default: 1.")


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
        exons = sorted(exons, key=lambda x: x[0])
        for i in range(len(exons)-1):
            #annotated[(chromosome,exons[i][1],exons[i+1][0]-1,strand)] = genes[tid]
            annotated.add(f"{chromosome}:{exons[i][1]+1}-{exons[i+1][0]-1}:{strand}")
    return annotated


def getInclusionCounts(filename, annotated=None):
    # Load counts into a numpy float32 matrix rather than a nested dict.
    # For 495 samples x 800K junctions, float32 uses ~1.6 GB vs ~30+ GB
    # for a nested dict of Python floats.
    with open(filename) as f:
        samples = f.readline().rstrip().split("\t")[1:]
        junctions = []
        rows = []
        for line in f:
            row = line.rstrip().split("\t")
            junctions.append(row[0])
            rows.append(row[1:])
    matrix = np.array(rows, dtype=np.float32)  # shape: (junctions, samples)
    junction_index = {j: i for i, j in enumerate(junctions)}
    sample_index = {s: i for i, s in enumerate(samples)}
    return matrix, samples, junctions, junction_index, sample_index


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


def _getJunctionsForSample(args_tuple):
    sample, coverageDirectory, annotated, allJunctions = args_tuple
    junctions = set()
    filename = os.path.join(coverageDirectory, f"{sample}_intron_coverage.txt")
    with open(filename) as percentileCoverage:
        for line in percentileCoverage:
            row = line.strip().split("\t")
            cluster = f"{row[0]}:{int(row[1])+1}-{row[2]}:{row[5]}"
            if not allJunctions and cluster not in annotated:
                continue
            junctions.add(cluster)
    return junctions


def getJunctions(samples, coverageDirectory, annotated, args):
    from multiprocessing import Pool
    worker_args = [(sample, coverageDirectory, annotated, args.allJunctions) for sample in samples]
    with Pool(args.numThreads) as pool:
        junctions = set()
        for result in pool.imap_unordered(_getJunctionsForSample, worker_args):
            junctions |= result
    return junctions


def _filterJunctionsForSample(args_tuple):
    sample, coverageDirectory, junctions, RSDthreshold = args_tuple
    filtered = set()
    filename = os.path.join(coverageDirectory, f"{sample}_intron_coverage.txt")
    with open(filename) as percentileCoverage:
        for line in percentileCoverage:
            row = line.strip().split("\t")
            cluster = f"{row[0]}:{int(row[1])+1}-{row[2]}:{row[5]}"
            if cluster not in junctions:
                continue
            covArray = np.array(row[-1].split(","), dtype=float)
            mean = np.mean(covArray)
            std = np.sqrt(np.mean((covArray - mean)**2))
            rsd = std / mean if mean > 0 else np.nan
            if rsd < float(RSDthreshold):
                filtered.add(cluster)
    return filtered


def getFilteredJunctions(samples, coverageDirectory, annotated, args):
    import time
    from multiprocessing import Pool
    t = time.time()
    junctions = getJunctions(samples, coverageDirectory, annotated, args)
    print(f"getJunctions complete: {len(junctions)} junctions. {time.time()-t:.1f}s")
    worker_args = [(sample, coverageDirectory, junctions, args.RSDthreshold) for sample in samples]
    with Pool(args.numThreads) as pool:
        filtered_junctions = set()
        for result in pool.imap_unordered(_filterJunctionsForSample, worker_args):
            filtered_junctions |= result
    print(f"RSD filtering complete: {len(filtered_junctions)} junctions retained. {time.time()-t:.1f}s")
    return filtered_junctions


def calculateIRforSample(sample, coverageDirectory, counts, clusters, junctions, args):
    matrix, samples, all_junctions, junction_index, sample_index = counts
    IR = {}
    RSD = {}
    si = sample_index[sample]
    filename = os.path.join(coverageDirectory, f"{sample}_intron_coverage.txt")
    with open(filename) as percentileCoverage:
        for line in percentileCoverage:
            row = line.strip().split("\t")
            cluster = f"{row[0]}:{int(row[1])+1}-{row[2]}:{row[5]}"
            if cluster not in junctions:
                continue
            median = float(row[4])
            # covArray is computed and discarded immediately per line
            # rather than stored in a per-sample dict, avoiding O(samples x junctions)
            # memory overhead for large cohorts.
            covArray = np.array(row[-1].split(","), dtype=float)
            mean = np.mean(covArray)
            std = np.sqrt(np.mean((covArray - mean)**2))
            RSD[cluster] = std / mean if mean > 0 else np.nan
            if cluster in junction_index:
                intronCount = float(matrix[junction_index[cluster], si])
                if not args.singleJunctionCalculation:
                    for mxCluster in clusters[cluster]:
                        if mxCluster in junction_index:
                            intronCount += float(matrix[junction_index[mxCluster], si])
                try:
                    IR[cluster] = median / (median + intronCount)
                except ZeroDivisionError:
                	IR[cluster] = np.nan
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

    if not args.allJunctions:
        print("Loading annotation...")
        annotated = getAnnotated(annotation)
        print(f"Annotation loaded: {len(annotated)} annotated junctions. {time.time()-start:.1f}s")
    else:
        annotated = None

    print("Gathering inclusion counts and clusters...")
    counts = getInclusionCounts(countFile)
    matrix, samples, all_junctions, junction_index, sample_index = counts
    clusters = None
    if not args.singleJunctionCalculation:
        clusters = getClusters(clusterFilename)
    print(f"Loaded {len(samples)} samples and {len(clusters) if clusters else 0} clusters. {time.time()-start:.1f}s")

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
    