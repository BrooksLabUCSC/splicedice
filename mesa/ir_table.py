#!/usr/bin/env python3


def add_parser(parser):
    parser.add_argument("-i","--inclusionCounts",
                        action="store",
                        help="")
    parser.add_argument("-c","--clusters",
                        action="store",
                        help="allClusters.tsv file with mutually exclusive clusters for each junction")

    parser.add_argument("-d","--coverageDirectory",
                        action="store",
                        help="")
    parser.add_argument("-o","--outputFilename",
                        action="store",
                        help="")

def getInclusionCounts(filename):
    with open(filename) as inclusionCounts:
        header = inclusionCounts.readline().strip().split("\t")
        counts = {sample:{} for sample in header[:-1]}
            
        for line in inclusionCounts:
            row = line.rstrip().split("\t")
            for i,count in enumerate(row[:-1]):
                counts[header[i]][row[-1]] = float(count)
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

def calculateIR(samples,coverageDirectory,counts,clusters):
    IR = {}
    coverage = {}
    junctions = set()
    RSD = {}
    for sample in samples:
        filename = os.path.join(coverageDirectory,f"{sample}_intron_coverage.txt")
        
        IR[sample] = {}
        coverage[sample] = {}
        RSD[sample] = {}
        with open(filename) as percentileCoverage:
            for line in percentileCoverage:
                row = line.strip().split("\t")
                cluster = f"{row[0]}:{row[1]}-{row[2]}"
                junctions.add(cluster)
                median = float(row[4])
                coverage[sample][cluster] = row[-1].split(",")
                covArray = np.array(coverage[sample][cluster])
                RSD[sample][cluster] = np.stddev(covArray) / np.mean(covArray)
                try:
                    intronCount = counts[sample][cluster]
                    for mxCluster in clusters[cluster]:
                        try:
                            intronCount += counts[sample][mxCluster]
                        except KeyError:
                            print("mxCluster",sample,cluster,mxCluster)
                    try:
                        IR[sample][cluster] = median/(median+intronCount)
                    except ZeroDivisionError:
                        IR[sample][cluster] = np.nan
                except KeyError:
                    print("cluster",sample,cluster)
    return junctions, IR, RSD


def writeIRtable(samples, outfilename, junctions, IR, RSD):
    tab = "\t"
    with open(outfilename,"w") as irTable:
        header = tab.join([f'{sample}\t{sample}_RSD'])
        irTable.write(f"Junction\t{header}\n")
        for junction in sorted(clusters):
            irValues = [f"{IR[sample][junction]:0.03f}\t{RSD[sample][junction]}" for sample in samples]
            irTable.write(f"{junction}\t{tab.join(irValues)}\n")   


def run_with(args):
    """ """
    import numpy as np
    import os

    countFile = args.inclusionsCounts
    clusterFilename = args.clusters
    coverageDir = args.coverageDirectory
    outfilename = args.outputFilename

    samples = [s for s in os.listdir(coverageDir) if s.endswith("intron_coverage.txt")]

    counts = getInclusionCounts(countFilename)
    clusters = getClusters(clusterFilename)

    junctions, IR, RSD = calculateIR(samples,coverageDirectory,counts,clusters)

    writeIRtable(samples, outfilename, junctions, IR, RSD)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    add_parser(parser)
    args = parser.parse_args()
    run_with(args)