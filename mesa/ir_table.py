#!/usr/bin/env python3

import numpy as np
import os
    
    
def add_parser(parser):
    parser.add_argument("-i","--inclusionCounts",
                        action="store",
                        help="")
    parser.add_argument("-c","--clusters",
                        action="store",
                        help="allClusters.tsv file with mutually exclusive clusters for each junction.")

    parser.add_argument("-d","--coverageDirectory",
                        action="store",
                        help="")
    parser.add_argument("-o","--outputPrefix",
                        action="store",
                        help="")
    parser.add_argument("-r","--makeRSDtable",
                        action="store_true",
                        help="Make a table of relative standard deviations in coverage across intron.")
    parser.add_argument("-s","--singleJunctionCalculation",
                        action="store_true",
                        help="Calculate IR value using individual junction counts, and not count of all junctions in cluster.")
def getInclusionCounts(filename):
    with open(filename) as inclusionCounts:
        header = inclusionCounts.readline().strip().split("\t")
        counts = {sample:{} for sample in header[1:]}
        for line in inclusionCounts:
            row = line.rstrip().split("\t")
            for i,count in enumerate(row[1:]):
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

def calculateIR(samples,coverageDirectory,counts,clusters,args):
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
                cluster = f"{row[0]}:{row[1]}-{row[2]}:{row[5]}"
                junctions.add(cluster)
                median = float(row[4])
                coverage[sample][cluster] = row[-1].split(",")
                covArray = np.array(coverage[sample][cluster]).astype(np.float)
                if args.makeRSDtable:
                    RSD[sample][cluster] = np.std(covArray) / np.mean(covArray)
                try:
                    intronCount = counts[sample][cluster]
                    if not args.singleJunctionCalculation:
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
                    break
    
    return junctions, IR, RSD


def writeIRtable(samples, outputPrefix, junctions, IR):
    tab = "\t"
    with open(f"{outputPrefix}_intron_retention.tsv","w") as irTable:
        irTable.write(f"Junction\t{tab.join(samples)}\n")
        for junction in sorted(junctions):
            irValues = [f"{IR[sample][junction]:0.03f}" for sample in samples]
            irTable.write(f"{junction}\t{tab.join(irValues)}\n")   
            
def writeRSDtable(samples, outputPrefix, junctions, RSD):
    tab = "\t"
    with open(f"{outputPrefix}_intron_retention_RSD.tsv","w") as rsdTable:
        header = tab.join([f'{sample}_RSD' for sample in samples])
        rsdTable.write(f"Junction\t{header}\n")
        for junction in sorted(junctions):
            rsd = [f"{RSD[s][junction]:0.03f}" for s in samples]
            rsdTable.write(f"{junction}\t{tab.join(rsd)}\n") 


def run_with(args):
    """ """
    import time
    start = time.time()

    countFile = args.inclusionCounts
    clusterFilename = args.clusters
    coverageDirectory = args.coverageDirectory
    outputPrefix = args.outputPrefix

    samples = [s.rstrip("intron_coverage.txt") for s in os.listdir(coverageDirectory) if s.endswith("intron_coverage.txt")]
    print("Gathering inclusion counts and clusters...")
    counts = getInclusionCounts(countFile)
    if not args.singleJunctionCalculation:
        clusters = getClusters(clusterFilename)
    else:
        clusters = None
    print("Calculating IR values...")
    junctions, IR, RSD = calculateIR(samples,coverageDirectory,counts,clusters,args)
    print("Done",time.time()-start)
    print("Writing output...")
    writeIRtable(samples, outputPrefix, junctions, IR)
    if args.makeRSDtable:
        writeRSDtable(samples, outputPrefix, junctions, RSD)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    add_parser(parser)
    args = parser.parse_args()
    run_with(args)