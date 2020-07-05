#!/usr/bin/env python3
"""
This program will process junction files to create a reference "transciptome" which will contain junctions that will process through kallisto

PATH=/private/groups/brookslab/bin/bedtools-2.25.0/bin/:$PATH

python3 sigma_psi.py -b sigma_psi_bam_manifest.tsv  -m DMSO_v_SSA100_allPSI.tsv  -j me_junctions.bed -c me_all_clusters2.tsv
"""

import time
start_time = time.time()

import pybedtools
import sys
import argparse
import os
import pysam
from multiprocessing import Pool


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
        self.bam = myCommandLine.args.bamManifest
        self.mesa = myCommandLine.args.mesaTable
        #self.output = myCommandLine.args.output
        #self.cluster = myCommandLine.args.clusters
        self.junctionFilename = myCommandLine.args.juncFile


    def getBAMFiles(self):

        print('getting paths for bam files')

        self.bams = {}
        self.bam_order = []

        with open(self.bam, 'r') as bams:
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
        self.fwdJunctions = []
        self.revJunctions = []

        fwdPercentiles = []
        revPercentiles = []

        with open(self.junctionFilename, 'r') as junc_file:
            for line in junc_file:
                line = line.strip().split('\t')
                chromosome = line[0]
                start = int(line[1])
                end = int(line[2])
                key = chromosome + ':' + str(start) + '-' + str(end)
                length = end-start
                strand = line[-1]
                percentiles = [ start+1,
                                int(start+(length*0.25)),
                                int(start+(length*0.50)),
                                int(start+(length*0.75)),
                                int(start+(length*0.99))]
                if strand is "+":
                    fwdPercentiles.append(percentiles)
                    self.fwdJunctions.append((chromosome,start,end,strand))

                else:
                    revPercentiles.append(percentiles)
                    self.revJunctions.append((chromosome,start,end,strand)
                self.percentiles[key] = percentiles + [line[-1]]
                self.junctions.append((chromosome,start,end,strand))

        self.fwdPercentiles = np.array(fwdPercentiles)
        self.revPercentiles = np.array(revPercentiles)


    def getPercentileValues(self):
        #create tmp bed file
        'creating bed file for percentiles'

        os.system("touch tmp_junc.bed")

        with open('tmp_junc.bed', 'w') as bed_file:
            for junc in self.percentiles:
                chr = junc.split(':')[0]
                for percentile in self.percentiles[junc][:-1]:
                    print (chr + '\t' + str(percentile-1) + '\t' + str(percentile) + '\t' + self.percentiles[junc][-1], file=bed_file)


    def getIRValues(self):
        for bam in self.bams:
            print('now calculating IR for %s' % (bam))
            new_file_name = bam + '_sigmaPSI_percentiles.txt'
            bam_bed = pybedtools.BedTool(self.bams[bam])
            bam_bed = bam_bed.bam_to_bed(split=True)
            junc_bed = pybedtools.BedTool('tmp_junc.bed')
            final_bed = junc_bed.coverage(bam_bed, d=True, output=new_file_name)

    def getIRValuesNew(self,sample):
        tab = "\t"
        filename = self.bams[sample]
        bamFile = pysam.AlignmentFile(filename,"rb")
        fwdCounts = np.zeros(self.fwdPercentiles.shape)
        revCounts = np.zeros(self.revPercentiles.shape)
        
        for read in bamFile.fetch():
            for block in read.get_blocks():
                if read.is_reverse:
                    revCounts += (self.revPercentiles > block[0]) & (self.revPercentiles < block[1])  
                else:
                    fwdCounts += (self.fwdPercentiles > block[0]) & (self.fwdPercentiles < block[1])
        
        fwdMedians = np.median(fwdCounts,axis=1)
        revMedians = np.median(revCounts,axis=1)

        with open(f"{sample}_sigmaPSI_percetiles.txt") as outfile:
            for i,junc in enumerate(self.junctions):
                if junc[3] is "+":
                    outfile.write(f"{junc[0]}\t{junc[1]}\t{junc[2]}\t.\t{fedMedians[i]}\t{strand}\t{','.join(self.revPecentiles[i,:])}\t{','.join(self.fwdCounts[i,:])}\n")
                else:
                    outfile.write(f"{junc[0]}\t{junc[1]}\t{junc[2]}\t.\t{revMedians}\t{strand}\t{','.join(self.revPercentiles[i,:])}\t{','.join(self.revCounts[i,:])}\n")




            
    def getIRValuesPool(self):
        numThreads = 16
        samples = list(self.bams.keys())
        with Pool(numThreads) as pool:
            for run in pool.map(self.getIRValuesNew, samples)


        



def main(myCommandLine=None):
    myCommandLine = CommandLine(myCommandLine)
    sigma_psi = sigmaPSI(myCommandLine)
    sigma_psi.getBAMFiles()
    sigma_psi.createPercentiles()
    sigma_psi.getClusters()
    sigma_psi.getPercentileValues()
    sigma_psi.getIRValues()
    #os.system("rm tmp_junc.bed")
    print("Your runtime was %s seconds." % (time.time() - start_time))

main()
