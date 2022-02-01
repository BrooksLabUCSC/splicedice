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
            
def get_counts(count_file):
    with open(count_file) as cf:
        header = cf.readline()
        counts = {}
        for line in cf:
            row = line.rstrip().split('\t')
            junction = row[0]
            counts[junction] = np.array(row[1:],dtype=float)
    return header,counts

def write_PS_values(clusters,header,counts,output_file):
    with open(output_file,"w") as psfile:
        psfile.write(header)
        for junction,overlaps in clusters.items():
            exclusion = counts[junction].copy()
            for overlap in overlaps:
                if overlap == '':
                    continue
                exclusion += counts[overlap]
            ps = counts[junction]/exclusion
            ps = '\t'.join([f"{x:0.3f}" for x in ps])
            psfile.write(f"{junction}\t{ps}\n")

def add_parser(parser):
    """ """
    parser.add_argument("--clusters","-c",
                        action="store",required=True,
                       help="allClusters.tsv file from MESA")
    parser.add_argument("--inclusion_counts","-i",
                        action="store",required=True,
                       help="inclusionCounts.tsv file from MESA")
    parser.add_argument("--output","-o",
                        action="store",required=True,
                       help="output filename") 
    
def run_with(args):
    """ Main program to calculate PS values"""
    print("Gathering clusters...")
    clusters = get_clusters(args.clusters)
    print("Gathering counts...")
    header,counts = get_counts(args.inclusion_counts)
    print("Calculating PS values...")
    write_PS_values(clusters,header,counts,args.output)
    


if __name__ == "__main__":
    import argparse 
    parser = argparse.ArgumentParser(description='Calculate PS values with inclusion count file.')
    add_parser(parser)
    args = parser.parse_args()
    run_with(args)