#!/usr/bin/env python3

"""counts_to_ps.py
Calculate PS values from inclusion count file and cluster file.
"""
import numpy as np

def get_clusters(cluster_file):
    with open(cluster_file) as cf:
        clusters = {}
        for line in cf:
            junction,overlaps = line.rstrip('\n').split('\t')
            clusters[junction] = overlaps.split(',')
    return clusters
            
def determine_clusters(counts_file):
    with open(counts_file) as cfile:
        junctions = []
        cfile.readline()
        for line in cfile:
            chromosome,coords,strand = line.split('\t',1)[0].split(":")
            start,stop = [int(x) for x in coords.split("-")]
            junctions.append((chromosome,start,stop,strand))
        clusters = {}
        chromosome,junction = None,None
        for junction in sorted(junctions, key = lambda x: (x[0],x[3],x[1],x[2])):
            if junction[0] != chromosome or junction[3] != strand:
                chromosome = junction[0]
                strand = junction[3]
                potentialOverlaps = []
            j_string = f"{junction[0]}:{junction[1]}-{junction[2]}:{junction[3]}"
            clusters[j_string] = []
            newPotentialOverlaps = [junction]
            for priorJunction in potentialOverlaps:
                if priorJunction[2] >= junction[1]: 
                    pj_string = f"{priorJunction[0]}:{priorJunction[1]}-{priorJunction[2]}:{priorJunction[3]}"
                    clusters[pj_string].append(j_string)
                    clusters[j_string].append(pj_string)  
                    newPotentialOverlaps.append(priorJunction)
            potentialOverlaps = newPotentialOverlaps
    return clusters
    
def get_counts(count_file):
    with open(count_file) as cf:
        header = cf.readline()
        counts = {}
        for line in cf:
            row = line.rstrip().split('\t')
            junction = row[0]
            counts[junction] = np.array(row[1:],dtype=float)
    return header,counts

def junctionStringToTuple(string):
    chromosome,coords,strand = string.split(":")
    start,end = [int(x) for x in coords.split('-')]
    return (chromosome,start,end,strand)

def writePsValues(clusters,header,counts,output_prefix):
    with open(f"{output_prefix}_allPS.tsv","w") as psfile:
        psfile.write(header)
        j_list = sorted(clusters.keys(),key = junctionStringToTuple)
        for junction in j_list:
            exclusion = counts[junction].copy()
            for overlap in clusters[junction]:
                if overlap == '':
                    continue
                exclusion += counts[overlap]
            ps = counts[junction]/exclusion
            ps = '\t'.join([f"{x:0.3f}" for x in ps])
            psfile.write(f"{junction}\t{ps}\n")

            
def writeClusters(clusters,output_prefix):
    """ """
    with open(f"{output_prefix}_allClusters.tsv","w") as clusterFile:
        for junction in sorted(clusters):
            clusterFile.write(f"{junction}\t{','.join(clusters[junction])}\n")
            
              
            
def add_parser(parser):
    """ """
    parser.add_argument("--clusters","-c",
                        action="store",default=None,
                       help="allClusters.tsv file from SPLICEDICE")
    parser.add_argument("--recluster","-r",
                        action="store_true",
                        help="Determine clusters from splice junctions in counts file")
    parser.add_argument("--inclusion_counts","-i",
                        action="store",required=True,
                       help="inclusionCounts.tsv file from SPLICEDICE")
    parser.add_argument("--output_prefix","-o",
                        action="store",required=True,
                       help="output filename path and prefix") 
    
def run_with(args):
    """ Main program to calculate PS values"""
    
    
    if args.clusters:
        print("Gathering clusters...")
        clusters = get_clusters(args.clusters)
    elif args.recluster:
        print("Determining clusters from counts file...")
        clusters = determine_clusters(args.inclusion_counts)
        writeClusters(clusters,args.output_prefix)
              
        
    print("Gathering counts...")
    header,counts = get_counts(args.inclusion_counts)
    print("Calculating PS values...")
    writePsValues(clusters,header,counts,args.output_prefix)
    print("Done.")
    


if __name__ == "__main__":
    import argparse 
    parser = argparse.ArgumentParser(description='Calculate PS values with inclusion count file.')
    add_parser(parser)
    args = parser.parse_args()
    run_with(args)