#!/usr/bin/env python3

import pysam
from multiprocessing import Pool
import os
import argparse




filenames = ("/scratch/dennisrm/sequin/mixa/sequin_mixa.Aligned.sortedByCoord.out.bam",
            "/scratch/dennisrm/sequin/mixb/sequin_mixb.Aligned.sortedByCoord.out.bam")


def parseManifest(manifestFilename,bedDirectory,outputPrefix):
    bams = []
    newManifest = []
    with open(manifestFilename, "r") as manifestFile:
        for line in manifestFile:
            sampleName, filename, metadata, condition = line.rstrip().split("\t")

            if filename.endswith(".bam"):
                bedFilename = os.path.join(bedDirectory,f"{os.path.basename(filename)[:-4]}.junc.bed")
                bams.append([sampleName,filename,metadata,condition,bedFilename])
                newManifest.append("\t".join([sampleName,bedFilename,metadata,condition]))
            else:
                newManifest.append(line)

    return bams, newManifest


def getJunctionsFromBam(sample):
    """
    """
    sampleName, filename, metadata, condition, bedFilename = sample
    genome = pysam.AlignmentFile(filename)
    counts = {}
    k = 0
    for read in genome.fetch():
        k += 1
        if k == 100:
            break
        if read.is_reverse:
            strand = "-"
        else:
            strand = "+"
        blocks = read.get_blocks()
        for i in range(len(blocks)-1):
            junction = (read.reference_name,blocks[i][1],blocks[i+1][0],strand)
            overhang = min(blocks[i][1]-blocks[i][0],blocks[i+1][1]-blocks[i+1][0])
            try:
                counts[junction][0] += 1
                counts[junction][1] = max(overhang,counts[junction][1])
            except KeyError:
                counts[junction] = [1,overhang]
    
    with open(bedFilename,"w") as bedOut:
        for junction in sorted(junction):
            chromosome,left,right,strand = junction
            bedOut.write(f"{chromosome}\t{left}\t{right}\t.\t{counts[junction]}\t{strand}\n")
            
    return filename, bedFilename, len(counts)


def bamsToBeds(bams, nThreads):

    with Pool(nThreads) as pool:
        for info in pool.map(getJunctionsFromBam,bams):
            filename,bedFilename,nJunctions = info
            print("bam:", filename)
            print("number of junctions found:",nJunctions)
            print("saved to bed:", bedFilename)


def writeNewManifest(newManifest):
    newManifestPath = os.path.join(path,f"{outputPrefix}_manifest")
    with open(newManifestPath,"w") as manifestFile:
        for line in newManifest:
            manifestFile.write("\t".join(line)+"\n")
    print("new manifest written to:", newManifestPath)


def parseArgs():
    """ """
    parser = argparse.ArgumentParser(description='Extract splice junctions from bam file and save to bed')           
    parser.add_argument("--manifest","-m")
    parser.add_argument("--output_prefix","-o")   
    parser.add_argument("--max_length",type=int,default=50000)
    parser.add_argument("--min_length",type=int,default=50)
    parser.add_argument("--min_unique",type=int,default=5)
    parser.add_argument("--min_overhang",type=int,default=5)
    parser.add_argument("--drim",action="store_true")
    parser.add_argument("--no_multimap",action="store_true")
    parser.add_argument("--filter",choices=["gtag_only","all"],default="gtag_only" )
    parser.add_argument("--number_threads","-n",type=int,default=1
    )

    return parser.parse_args()
    
def main():

    args = parseArgs()
    manifestFilename = args.manifest
    outputPrefix = args.output_prefix
    nThreads = args.number_threads

    bedDirectory = os.path.join(os.getcwd(),f"{outputPrefix}_junction_beds")
    os.mkdir(bedDirectory)

    bams, newManifest = parseManifest(manifestFilename,bedDirectory,outputPrefix)

    bamsToBeds(bams,nThreads)

    writeNewManifest(newManifest)


if __name__ == "__main__":
    main()
