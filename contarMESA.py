#!/usr/bin/env python3


########################################################################
# File: countMESA.py
#  executable: countMESA.py
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 06/20/2018 Created
#
########################################################################

########################################################################
# Hot Imports & Global Variable
########################################################################


import os, sys
import numpy as np
from multiprocessing import Pool, Manager
from scipy import stats
from tqdm import *

########################################################################
# CommandLine
########################################################################

class CommandLine(object) :
    '''
    Handle the command line, usage and help requests.
    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    and a standard usage and help,
    attributes:
    myCommandLine.args is a dictionary which includes each of the available command line arguments as
    myCommandLine.args['option'] 
    
    methods:
    
    '''
    
    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''
        import argparse
        self.parser = argparse.ArgumentParser(description = 'countMESA.py - Contar la mesa',
                                             epilog = 'Please feel free to forward any questions or concerns to /dev/null', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s -m junction_beds_manifest.tsv ')
        # Add args
        self.parser.add_argument('-m', '--bed_manifest', type=str, action = 'store', required=True, help='List of bed files containing stranded junctions.')
        self.parser.add_argument('-j', '--junctions_bed', type=str, action = 'store', required=True, help='Junction bed file.')
        self.parser.add_argument('-c', '--clusters_table', type=str, action = 'store', required=True, help='cluster table tsv.')
        self.parser.add_argument('-t', '--threads', type=int, action = 'store', required=False, default=1, help='threads.')
        
        
        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))

########################################################################
# Clusters & Junctions
########################################################################

class Junction(object):
    def __init__(self, chrom=None, c1=None, c2=None, strand=None):
        self.chrom = chrom 
        self.c1 = c1
        self.c2 = c2
        self.strand = strand
        self.mxes = list() 
        self.name = "%s:%s-%s" % (chrom, c1, c2)
        self.quant = list()

class Sample(object):
    def __init__(self, name=None):
        self.sid = name
        self.introns = dict()

########################################################################
# Helper Functions
# 
# 
########################################################################

def bedToEvents(f,sampleID):
    '''
    takes bed file and returns a dictionary
    '''
    sampleObj = Sample(sampleID)
    
    with open(f,'r') as lines: 
        for line in lines:
            cols = line.split()
            name = "%s:%s-%s" % (cols[0],cols[1],cols[2])
            if name in junctionSet:
                sampleObj.introns[name] = int(cols[4])
            
    return sampleObj

def runCMD1(x):
    '''
    mutiprocessing helper function to run
    multiple instances of bedToCoordSet
    '''
    root, bed= x
    return bedToEvents(bed,root)

def runCMD2(x):
    '''
    mutiprocessing helper function to run
    multiple instances of bedToCoordSet
    '''
    return clusterQuant(x)


def clusterQuant(intron):
    '''
    stuff
    '''
    intron, samples = intron
    for x in samples:
        inclusion = x.introns.get(intron.name,0)
        exclusion = sum([x.introns.get(y,0) for y in intron.mxes])
        intron.quant.append("%s:%s" % (inclusion,exclusion))
    return intron


########################################################################
# Main
# Here is the main program
# 
########################################################################

def main():
    '''
    TDB
    '''
    myCommandLine = CommandLine()
    
    bedList = myCommandLine.args['bed_manifest']
    junctionBed = myCommandLine.args['junctions_bed']
    clusters = myCommandLine.args['clusters_table']
    threads = myCommandLine.args['threads']
    


    # Get list of bed file names for multiprocessing step.
    beds = list()
    with open(bedList,'r') as fnames:
        for fdata in fnames:
            root,fname = fdata.split()
            beds.append((root,fname))

    


    # Make junctionSet global for later steps.
    global junctionSet, samples
    junctionSet = dict()

    # Make junction objects
    n = 0
    with open(junctionBed,'r') as lines:
        for line in tqdm(lines, total=864511, desc="Initializing junction counts"):
            j = line.rstrip().split()
            c1, c2 = map(int, [j[1], j[2]])
            jid = "%s:%s-%s" % (j[0],j[1],j[2])
            junctionSet[jid] = Junction(j[0], c1, c2,j[5])
            junctionSet[jid].quant = np.zeros(len(beds), dtype=int)
            n += 1


    #p = Pool(threads)
    #finalSet = set()

    # samples = list()
    # for result in tqdm(p.imap_unordered(runCMD1, beds[:]), total=len(beds[:]), desc="Reading sample bed files for intron counts"):
    #     samples.append(result)


    num = 0
    for bed in tqdm(beds[:], total=len(beds[:]), desc="Quantifying junction usage per sample"):
        root,fname = bed
        with open(fname,'r') as lines:
            for line in lines:
                cols = line.split()
                name = "%s:%s-%s" % (cols[0],cols[1],cols[2])
                if name in junctionSet:
                    obj = junctionSet[name]
                    obj.quant[num] = int(cols[4])
                else:
                    continue
        num += 1
                

    #sys.exit(1)
    #p.close()
    #p.join()

    lineNum = 0
    with open(clusters,'r') as lines:
        for i in lines:
            lineNum += 1

    with open("default_out_new",'w') as out1:
        print("\t".join([x[0] for x in beds]), "intron", "mxes", file=out1)
        with open(clusters,'r') as lines:
            for line in tqdm(lines, total=lineNum, desc="Compute inclusion/exclusion counts"):
                intron, mxes = line.rstrip().split()
                mxes = mxes.split(",")
                
                inclusions = junctionSet[intron].quant
                exclusions = np.asarray([junctionSet[name].quant for name in mxes])

                exclusions = exclusions.sum(axis=0)
                quants = ["%s:%s" % (x[0],x[1]) for x in zip(inclusions,exclusions)]
                print("\t".join(quants), intron, ",".join(mxes), file=out1)
    


if __name__ == "__main__":
    main()      