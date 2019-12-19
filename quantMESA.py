#!/usr/bin/env python3


########################################################################
# File: quantMESA.py
#  executable: quantMESA.py
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
        self.parser = argparse.ArgumentParser(description = 'quantMESA.py - Quantify mutually exclusive clusters.',
                                             epilog = 'Please feel free to forward any questions or concerns to /dev/null', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s -m junction_beds_manifest.tsv -j step1_junctions.bed -c step1_clusters.bed ')
        # Add args
        self.parser.add_argument('-m', '--bed_manifest',  type=str, action = 'store', required=True, help='List of bed files containing stranded junctions.')
        self.parser.add_argument('-j', '--junctions_bed', type=str, action = 'store', required=True, help='Junction bed file.')
        self.parser.add_argument('--inclusion_thresh',    type=int, action = 'store', required=False, default=10, help='Filter events with less than N reads (default 10)')
        self.parser.add_argument('--psi_thresh',          type=int, action = 'store', required=True, default=20, help='Do not compute PSI for samples with less than N reads (default 20)')
        

        self.parser.add_argument('-c', '--clusters_table', type=str, action = 'store', required=True, help='cluster table tsv.')
        self.parser.add_argument('-o', '--output_prefix', type=str, action = 'store', required=True, help='Output prefix')
        #self.parser.add_argument('-p', '--threads', type=int, action = 'store', required=False, default=1, help='threads.')
        self.parser.add_argument('--compressed', action = 'store_true', required=False, default=False, help='NPZ compression')
        
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
    
    bedList     = myCommandLine.args['bed_manifest']
    junctionBed = myCommandLine.args['junctions_bed']
    clusters    = myCommandLine.args['clusters_table']
    comp        = myCommandLine.args['compressed']
    outPrefix   = myCommandLine.args['output_prefix']
    incT = myCommandLine.args['inclusion_thresh']
    psiT = myCommandLine.args['psi_thresh']

    # Get list of bed file names for multiprocessing step.
    beds   = list()
    groups = dict()
    with open(bedList,'r') as fnames:
        #next(fnames)
        for num,fdata in enumerate(fnames,0):
            root,fname,group1,group2 = fdata.rstrip().split("\t")
            group1,group2 = group1.replace(" ",""), group2.replace(" ","")
            beds.append((root,fname))
            if (group1,group2) not in groups:
                groups[(group1,group2)] = list()
            
            groups[(group1,group2)].append(int(num))

    # Make junctionSet global for later steps.
    global junctionSet
    junctionSet = dict()

    # Make junction objects
    n = 0
    with open(junctionBed,'r') as lines:
        for line in tqdm(lines, desc="Initializing junction counts", position=1):
            j = line.rstrip().split()
            c1, c2 = map(int, [j[1], j[2]])
            jid = "%s:%s-%s" % (j[0],j[1],j[2])
            junctionSet[jid] = Junction(j[0], c1, c2,j[5])
            junctionSet[jid].quant = np.zeros(len(beds), dtype=int)
            n += 1


    num = 0
    for bed in tqdm(beds[:], total=len(beds[:]), desc="Quantifying junction usage per sample", position=1):
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
                

    dataQ = list()
    dataP = list()
    events = list()
    beds = np.array(beds)
    if not comp:
        incOut = open('%s_inclusionCounts.tsv' % outPrefix ,'w')
        psiOut = open('%s_allPSI.tsv' % outPrefix ,'w')
        print("\t".join(beds[:,0]),"cluster",sep="\t",file=incOut)
        print("\t".join(beds[:,0]),"cluster",sep="\t",file=psiOut)

    with open(clusters,'r') as lines:
        for line in tqdm(lines,  desc="Compute inclusion/exclusion counts", position=1):
            intron, mxes = line.rstrip().split()
            mxes = mxes.split(",")
            
            inclusions = junctionSet[intron].quant
            exclusions = np.asarray([junctionSet[name].quant for name in mxes])

            total = exclusions.sum(axis=0) + inclusions

            # lets filter...

            if np.max(inclusions)<incT:
                continue

            events.append(intron)
            quants = np.asarray(["%s:%s" % (x[0],x[1]) for x in zip(inclusions,total)])
            psiQuants = inclusions/total
            psiQuants[total < psiT] = np.nan 
            dataQ.append(quants)
            dataP.append(psiQuants)
            if not comp:
                print("\t".join("%.2f" % x for x in psiQuants), intron, sep="\t", file=psiOut)
                print("\t".join(quants), intron, sep="\t", file=incOut)
    if comp:
        cols,rows,data1 = np.array(beds), np.array(events), np.array(dataQ)
        data2 = np.array(dataP,dtype=np.float32)
        np.savez_compressed("%s_inclusionCounts.npz" % outPrefix,cols=cols,rows=rows,data=data1)
        np.savez_compressed("%s_allPSI.npz" % outPrefix ,cols=cols,rows=rows,data=data2)

    if not comp:
        incOut.close()
        psiOut.close()



    print()

if __name__ == "__main__":
    main()      