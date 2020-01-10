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
import psutil

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
        self.parser.add_argument('-m', '--bed_manifest', type=str, action = 'store', required=True, help='List of bed files containing stranded junctions.')
        self.parser.add_argument('-j', '--junctions_bed', type=str, action = 'store', required=True, help='Junction bed file.')
        self.parser.add_argument('-c', '--clusters_table', type=str, action = 'store', required=True, help='cluster table tsv.')
        self.parser.add_argument('-o', '--output_prefix', type=str, action = 'store', required=True, help='Output prefix')
        self.parser.add_argument('-p', '--threads', type=int, action = 'store', required=False, default=1, help='threads.')
        self.parser.add_argument('--uncompressed_only', action = 'store_false', required=False, default=True, help='NPZ compression')
        self.parser.add_argument('--drimTable', action = 'store_true', required=False, default=False, help='Generate a table to use with DRIM-Seq.')
        

        self.parser.add_argument('--inclusion_thresh',    type=int, action = 'store', required=False, default=10, help='Filter events with less than N reads (default 10)')
        self.parser.add_argument('--psi_thresh',          type=int, action = 'store', required=False, default=20, help='Do not compute PSI for samples with less than N reads (default 20)')
        

             



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

def memory_usage_psutil():
    # return the memory usage in MB
    process = psutil.Process(os.getpid())
    mem = process.memory_info()[0] / float(1024.0 ** 2)
    return mem


def runCMD2(x):
    '''
    mutiprocessing helper function to run
    multiple instances of bedToCoordSet
    '''
    return clusterQuant(x)


def clusterQuant(data):
    '''
    stuff
    '''
    root, fname, order = data

    order = np.load(order)['introns']
    d = {x:np.nan for x in order}
    with open(fname) as fin:
        for i in fin:

            cols = i.rstrip().split()
            intronID = "%s:%s-%s" % (cols[0],cols[1],cols[2])
            if intronID in d:
                d[intronID] = int(cols[4])
    return [root, np.array([d[x] for x in order], dtype=np.float32 )]


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
    comp        = myCommandLine.args['uncompressed_only']
    outPrefix   = myCommandLine.args['output_prefix']
    threads     = myCommandLine.args['threads']
    incT = myCommandLine.args['inclusion_thresh']
    psiT = myCommandLine.args['psi_thresh']
    drim = myCommandLine.args['drimTable']
    # Get list of bed file names for multiprocessing step.
    beds   = list()
    with open(bedList,'r') as fnames:
        for num,fdata in enumerate(fnames,0):
            root,fname,group1,group2 = fdata.rstrip().split("\t")
            group1,group2 = group1.replace(" ",""), group2.replace(" ","")
            beds.append((root,fname,'order.npz'))

    # Make junction objects
    n = 0
    order = list()
    orderDict = dict()
    with open(junctionBed,'r') as lines:
        for line in tqdm(lines, desc="Initializing junction counts", position=1):
            j = line.rstrip().split()
            c1, c2 = map(int, [j[1], j[2]])
            jid = "%s:%s-%s" % (j[0],j[1],j[2])
            order.append(jid)
            orderDict[jid] = n
            n += 1

    # this compressed order array will be used to re-order junction counts
    # from individual bed files
    order = np.array(order)
    np.savez_compressed('order.npz',introns=order)
    allData = list()
    samps  = list()
    with Pool(threads) as p:
        for i in tqdm(p.imap_unordered(runCMD2, beds), total=len(beds), desc="Computing counts", position = 1):
            allData.append(i[-1])
            samps.append(i[0])

    #print("current mem %s" % memory_usage_psutil())
    allData = np.nan_to_num(np.transpose(np.array(allData)))
    #print(allData.shape)
    
    dataQ  = list()
    dataP  = list()
    events = list()
    incOut  = open('%s_inclusionCounts.tsv' % outPrefix ,'w')
    psiOut  = open('%s_allPSI.tsv' % outPrefix ,'w')

    if drim: 
        drimOut = open('%s_drimTable.tsv' % outPrefix ,'w')
    
    print("\t".join(samps),"cluster",sep="\t",file=incOut)
    print("\t".join(samps),"cluster",sep="\t",file=psiOut)
    
    if drim:
        print("gene","feature_id","\t".join(samps),sep="\t",file=drimOut)

    cluster_num = 0
    with open(clusters,'r') as lines:
        for line in tqdm(lines,  desc="Compute inclusion/exclusion counts", position=1):
            intron, mxes = line.rstrip().split()
            mxes = mxes.split(",")
            
            intronPos = orderDict[intron]
            mxePos = [orderDict[x] for x in mxes]

            inclusions = allData[intronPos]
            exclusions = allData[mxePos]

            total = exclusions.sum(axis=0) + inclusions

            # lets filter...

            # if np.max(inclusions)<incT:
            #     continue

            events.append(intron)
            psiQuants = inclusions/total
            psiQuants[total < psiT] = np.nan 
            dataQ.append(inclusions)
            dataP.append(psiQuants)
            
            print("\t".join("%.2f" % x for x in psiQuants), intron, sep="\t", file=psiOut)
            print("\t".join(str(x) for x in inclusions), intron, sep="\t", file=incOut)

            if drim:
                print("cl_%s_%s\t%s_%s" % (cluster_num,intron,intron,cluster_num), "\t".join(str(x) for x in inclusions), sep="\t", file=drimOut)
                # print other intron values too
                for mxeNum, x in enumerate(exclusions):
                    print("cl_%s_%s\t%s_%s" % (cluster_num,intron,mxes[mxeNum],cluster_num), "\t".join(str(y) for y in x), sep="\t", file=drimOut)
                
            cluster_num += 1

    cols,rows,data1 = np.array(samps), np.array(events), np.array(dataQ)
    data2 = np.array(dataP,dtype=np.float32)
    
    if comp:
        np.savez_compressed("%s_inclusionCounts.npz" % outPrefix,cols=cols,rows=rows,data=data1)
        np.savez_compressed("%s_allPSI.npz" % outPrefix ,cols=cols,rows=rows,data=data2)

    incOut.close()
    psiOut.close()
    drimOut.close()


    #print()

if __name__ == "__main__":
    main()      