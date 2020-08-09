#!/usr/bin/env python3
"""
Process alignment files (BAMs) to determine coverage at positions along introns, as defined in junction file.

python3 intron_coverage.py -b bam_manifest.tsv  -m output_allPSI.tsv  -j mesa_junctions.bed
"""


def add_parser(parser):
    """
    """

    parser.add_argument('-b', '--bamManifest', 
                        action = 'store', required=True, 
                        help='tab-separated list of bam files for all samples')
    parser.add_argument('-m', '--mesaTable', 
                        action = 'store', required=True, 
                        help='mesa allPSI.tsv output from mesa')
    parser.add_argument('-j', '--junctionFile', 
                        action = 'store', required=True, 
                        help='mesa junction.bed output from mesa')
    parser.add_argument('-s', '--binSize', 
                        action='store', default=500000, 
                        help='chromosome bins for processing (must exceed intron length)')
    parser.add_argument('-n', '--numThreads', 
                        action='store', default=1, 
                        help='number of BAM-parsing threads to run concurrently')

    parser.add_argument('-o', '--outputDir',
                        action = 'store', required=True, 
                        help='directory for outputting intron coverage count files')

class IntronCoverage():
    def __init__(self, manifest,mesa,junctionFile,binSize,numThreads):
        """
        Instantiate object attributes
        """
        self.manifest = manifest
        self.mesa = mesa
        self.junctionFilename = junctionFile
        self.binSize = binsSize
        self.numThreads = numThreads

        self.start = time.time()


    def parseManifest(self):

        print('getting paths for bam files')

        self.bams = {}
        self.bam_order = []

        with open(self.manifest, 'r') as bams:
            for line in bams:
                line = line.strip().split('\t')
                sample = line[0]
                path = line[1]
                self.bams[sample] = path
                self.bam_order.append(sample)

        self.header_list = [key+'%IR' for key in self.bams.keys()]

    def getIntronPercentiles(self):
        """
        Create a transcriptome in fasta format
        """

        print('creating junction percentiles')
        
        self.percentiles = {}
        self.junctions = {}
        self.chromosomes = set()
        self.binMax = {}
        
        with open(self.junctionFilename, 'r') as junc_file:
            for line in junc_file:
                line = line.strip().split('\t')
                chromosome = line[0]
                start = int(line[1])
                end = int(line[2])
                length = end-start
                strand = line[-1]
                percs = [ start+1,
                                int(start+(length*0.25)),
                                int(start+(length*0.50)),
                                int(start+(length*0.75)),
                                int(start+(length*0.99))]
                
                startbin = start // self.binSize
                chrom = (strand,chromosome,startbin)
                self.chromosomes.add(chrom)
                
                try:
                    self.percentiles[chrom].append(percs)
                    self.junctions[chrom].append((chromosome,start,end,strand))
                    self.binMax[chrom] = max(self.binMax[chrom],percs[-1])

                except KeyError:
                    self.percentiles[chrom] = [percs]
                    self.junctions[chrom] = [(chromosome,start,end,strand)]
                    self.binMax[chrom] = percs[-1]


        self.percentiles = {chrom:np.array(a) for chrom,a in self.percentiles.items()}


    def getCoverage(self,sample):
        tab = "\t"
        filename = self.bams[sample]
        
        percentiles = self.percentiles.copy()
        
        counts = {chrom:np.zeros(a.shape,dtype=int) for chrom,a in percentiles.items()}
        
        print(sample,"starting",time.time()-self.start) 
        
        k = 0
                        
        bamFile = pysam.AlignmentFile(filename,"rb")
        
        lefts = {}
        rights = {}
       
        for read in bamFile.fetch(until_eof=True):
            #k += 1
            #if k == 10000:
            #    break
            
            chromosome = read.reference_name
            
            if read.is_reverse:
                strand = "-"
            else:
                strand = "+"
            for block in read.get_blocks():
                chrom = (strand,chromosome,block[0]//self.binSize)
                try:
                    lefts[chrom][block[0]] += 1
                except KeyError:
                    try:
                        lefts[chrom][block[0]] = 1
                    except KeyError:
                        lefts[chrom] = {}
                        lefts[chrom][block[0]] = 1

                try:
                    rights[chrom][block[1]] += 1
                except KeyError:
                    try:
                        rights[chrom][block[1]] = 1
                    except KeyError:
                        rights[chrom] = {}
                        rights[chrom][block[1]] = 1

        print(sample,"collected",time.time()-self.start) 

        counts = {chrom:np.zeros(shape=a.shape,dtype=int) for chrom,a in percentiles.items()}

        for chrom in lefts.keys():
            strand,chromosome,startbin = chrom

            for left,count in lefts[chrom].items():

                # Add to bin
                try:
                    counts[chrom] += count * (percentiles[chrom] > left)
                except KeyError:
                    pass

                # Previous bin
                try:
                    preChrom = (strand,chromosome,startbin-1)
                    if self.binMax[preChrom] >= left:
                        counts[preChrom] += count * (percentiles[preChrom] > left)
                except KeyError:
                    pass


            for right,count in rights[chrom].items():
                try:
                    counts[chrom] -= count * (percentiles[chrom] > right)
                except KeyError:
                    pass

                # Previous bin
                try:
                    preChrom = (strand,chromosome,startbin-1)
                    if self.binMax[preChrom] >= right:
                        counts[preChrom] -= count * (percentiles[preChrom] > right)
                except KeyError:
                    pass

                # Next bin
                try:
                    if right // self.binSize != startbin:
                        nextChrom = (strand,chromosome, right//self.binSize)
                        counts[nextChrom] += count * (percentiles[nextChrom] < right)  
                except KeyError:
                    pass

        ####
        medians = {chrom:np.median(a,axis=1) for chrom,a in counts.items()}
        print(sample,"counted",time.time()-self.start) 
        
        # Making list of junctions and array indices, sorting by coordinates
        junctions = []
        for chrom in self.junctions:
            junctions.extend(list(enumerate(self.junctions[chrom])))

        junctions.sort(key = lambda j: j[1])
        
        # Writing junctions and percentiles to bed-like text file
        with open(f"{sample}_intron_coverage.txt","w") as outfile:
            
            for i,junction in junctions:
                chromosome,left,right,strand = junction
                chrom = (strand, chromosome, left//self.binSize)
                juncPercentiles = ",".join(percentiles[chrom][i,:].astype(str))
                juncCounts = ','.join(counts[chrom][i,:].astype(str))
                median = medians[chrom][i]
                outfile.write(f"{chromosome}\t{left}\t{right}\t.\t{median}\t{strand}\t{juncPercentiles}\t{juncCounts}\n")
        print(sample,"done",time.time()-start) 
            
    def getCoveragePool(self):
        samples = list(self.bams.keys())
        with Pool(self.numThreads) as pool:
            for run in pool.imap_unordered(self.getCoverage, samples):
                pass
                     

def run_with(args):
    """  """
    import time
    import sys
    import os
    import pysam
    import numpy as np
    from multiprocessing import Pool

    manifest = args.bamManifest
    mesa = args.mesaTable
    junctionFilename = args.juncFile
    binSize = args.binSize
    numThreads = args.numThreads

    intronCoverage = IntronCoverage(manifest,mesa,junctionFile,binSize,numThreads)
    intronCoverage.parseManifest()
    intronCoverage.getIntronPercentiles()
    intronCoverage.getCoveragePool()
    print("Your runtime was %s seconds." % (time.time() - start))

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    add_parser(parser)
    args = parser.parse_args()
    run_with(args)

