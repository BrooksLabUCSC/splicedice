#!/usr/bin/env python3
"""
This program will process junction files to create a reference "transciptome" which will contain junctions that will process through kallisto

PATH=/private/groups/brookslab/bin/bedtools-2.25.0/bin/:$PATH

python3 sigma_psi.py -b sigma_psi_bam_manifest.tsv  -m DMSO_v_SSA100_allPSI.tsv  -j me_junctions.bed -c me_all_clusters2.tsv
"""

import time
start = time.time()

import pybedtools
import sys
import argparse
import os
import pysam
from multiprocessing import Pool
import numpy as np


class CommandLine() :
    """
    Takes command line options
    """
    def __init__(self, inOpts=None) :
            """
            CommandLine constructor.
            Implements a parser to interpret the command line argv string using argparse.
            """
            self.parser = argparse.ArgumentParser(description = 'runRTest.py arguments',
                                                 epilog = 'For more help contact: ',
                                                 add_help = True, #default is True
                                                 prefix_chars = '-',
                                                 usage = '%(prog)s [options] -option1[default]'
                                                 )
             # Add args
            self.parser.add_argument('-b', '--bamManifest', action = 'store', required=True, help='merged bam file for all samples')
            self.parser.add_argument('-m', '--mesaTable', action = 'store', required=True, help='mesa allPSI.tsv output from quantMESA.py')
            self.parser.add_argument('-j', '--juncFile', action = 'store', required=True, help='mesa allPSI.tsv output from constructMESA.py')
            #self.parser.add_argument('-o', '--output', action = 'store', required=True, help='output PSI table with % IR')
            #self.parser.add_argument('-c', '--clusters', action = 'store', required=True, help='me_clusters.tsv output from constructMESA.py')

            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)

class sigmaPSI():
    def __init__(self, myCommandLine):
        """
        Instantiate object attributes
        """
        self.manifest = myCommandLine.args.bamManifest
        self.mesa = myCommandLine.args.mesaTable
        self.junctionFilename = myCommandLine.args.juncFile
        self.binSize = 500000


    def getBAMFiles(self):

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

    def createPercentiles(self):
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


    def getIRValuesNew(self,sample):
        tab = "\t"
        filename = self.bams[sample]
        
        percentiles = self.percentiles.copy()
        
        counts = {chrom:np.zeros(a.shape,dtype=int) for chrom,a in percentiles.items()}
        
        print(sample,"starting",time.time()-start) 
        
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

        print(sample,"collected",time.time()-start) 

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
        print(sample,"counted",time.time()-start) 
        
        # Making list of junctions and array indices, sorting by coordinates
        junctions = []
        for chrom in self.junctions:
            junctions.extend(list(enumerate(self.junctions[chrom])))

        junctions.sort(key = lambda j: j[1])
        
        # Writing junctions and percentiles to bed-like text file
        with open(f"{sample}_sigmaPSI_percentiles.txt","w") as outfile:
            
            for i,junction in junctions:
                chromosome,left,right,strand = junction
                chrom = (strand, chromosome, left//self.binSize)
                juncPercentiles = ",".join(percentiles[chrom][i,:].astype(str))
                juncCounts = ','.join(counts[chrom][i,:].astype(str))
                median = medians[chrom][i]
                outfile.write(f"{chromosome}\t{left}\t{right}\t.\t{median}\t{strand}\t{juncPercentiles}\t{juncCounts}\n")
        print(sample,"done",time.time()-start) 
            
    def getIRValuesPool(self):
        numThreads = 8
        samples = list(self.bams.keys())
        with Pool(numThreads) as pool:
            for run in pool.imap_unordered(self.getIRValuesNew, samples):
                pass
                     

def main(myCommandLine=None):
    myCommandLine = CommandLine(myCommandLine)
    sigma_psi = sigmaPSI(myCommandLine)
    sigma_psi.getBAMFiles()
    sigma_psi.createPercentiles()
    sigma_psi.getIRValuesPool()
    print("Your runtime was %s seconds." % (time.time() - start))

main()
