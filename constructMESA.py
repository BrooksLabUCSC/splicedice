#!/usr/bin/env python3


########################################################################
# File: constructMESA.py
#  executable: constructMESA.py
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 06/20/2018 Created
#
########################################################################

#IMPORT
import time
import os, sys
from tqdm import *
from multiprocessing import Pool
import tempfile
import uuid 
import pybedtools

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
        self.parser = argparse.ArgumentParser(description = 'constructMESA.py - Construct splicing table',
                                             #epilog = 'Please feel free to forward any questions or concerns to /dev/null', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s -m junction_beds_manifest.tsv -f genome.fa [other options]',
                                             formatter_class=argparse.RawTextHelpFormatter)

        # Add args
        self.parser.add_argument('-m', '--bed_manifest', action = 'store', required=True, 
                                help='List of junction bed files.')
        self.parser.add_argument('-f', '--fasta_in',     action = 'store', required=True, 
                                help='BedTools indexed genome gasta.')

        self.parser.add_argument('-o', '--out_prefix',   action = 'store', required=False, default='me',
                                help='Output file prefix.')
       

        self.parser.add_argument('--filter',             action = 'store', required=False, default='gtag_only', type=lambda x: self.is_valid_filter(self.parser, x),
                                help='Filter splice sites that\n'+
                                    'filter non GT-AG splice sites : "gtag_only" (default),\n'+
                                    'filter non GT/GC-AG and AT-AC introns: "gc_at",\n'+
                                    'keep all splice sites: "all"')
   
        self.parser.add_argument('--resolve_strands',    action = 'store_true', required=False, default=False, 
                                help='Resolve splice site strand based on dinucleotide.')

        self.parser.add_argument('--support_threshold',  action = 'store', required=False, default = 10, type=int,
                                help='Filter junctions < non-nomrlized read counts (10)')
        self.parser.add_argument('--max_length',         action = 'store', required=False, default = 50000, type=int, 
                                help='Filter junctions > N length (50,000)')
        self.parser.add_argument('--min_length',     action = 'store', required=False, default = 50, type=int, 
                                help='Filter junctions < N length (50)')
        self.parser.add_argument('--min_samp_threshold', action = 'store', required=False, default = 0.01, type=float, 
                                help='Filter junctions found in < N samples (0.1)')
        
        self.parser.add_argument('-p', '--threads',      action = 'store', required=False, type=int, default=2, 
                                help='Number of threads.')
        
        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))


    def is_valid_filter(self, parser, arg):
        if arg not in ['gtag_only','gc_at','all']:
            parser.error("Filter argument %s is invalid!" % arg)
            sys.exit(1)
        else:
            return arg

########################################################################
# Helper Functions
# 
# 
########################################################################

def bedToCoordSet(f):
    '''
    takes bed file and returns a set of junctions
    '''
    fname, tdir = f
    newOut = fname.split("/")[-1]
    outName = os.path.join(tdir,"%s.tmp" % newOut)
    jSet = set()
    with open(outName,'wb+') as fout:
        with open(fname,'r') as lines:
            for line in lines:
                cols = line.rstrip().split()
                chrom, c1, c2, i, score, strand = cols[:6]
                c1,c2,score = map(int, [c1,c2,score])
                if score<readThresh or c2 - c1 > lengthThresh:
                    continue
                cols[4] = '0'
                outString = str.encode("\t".join(cols[:6]) + "\n")
                fout.write(outString)


def runCMD(x):
    '''
    mutiprocessing helper function to run
    multiple instances of bedToCoordSet
    '''
    return bedToCoordSet(x)

def makeClusters(intersection):
    '''
    takes a bedtools intersection with -wa -wb options and returns
    a list of overlaps between -a and all -b

    e.g. of intersection

    GL000008.2  164884  170271  .   0   -   GL000008.2  164884  170271  .   0   -
    '''
    mxeDict = dict()
    for i in intersection:
        left = "%s:%s-%s" % (i[0],i[1],i[2])
        right = "%s:%s-%s" % (i[6],i[7],i[8])
        if left == right:
            continue
        else:
            if left not in mxeDict:
                mxeDict[left] = list()
            mxeDict[left].append(right)

    return mxeDict

def checkFileExist(f):
    if os.path.isfile(f):
        pass
    else:
        print("** ERR: %s file does not exist. **" % f, file=sys.stderr)
        sys.exit(1)

def juncFilesFromManifest(manifest,tdir):
    '''function takes in 2 column manifest
    and returns list of file names (second column)'''

    beds = list()
    with open(manifest,'r') as lines:
        for line in lines:
            data = line.rstrip().split("\t")
            beds.append((data[1],tdir))
    return beds


def filterJunctions(bed,filta,correct,genome):
    '''
    reads in junction bed file, and filters based on ss motif.
    also resolves if true
    '''
    finalBed = str()
    btJuncs = pybedtools.BedTool(bed)
    dinucSeq = btJuncs.sequence(fi=genome, s=True, tab=True)
    if correct:
        print("** Under construction **",file=sys.stderr)
        sys.exit(1)
        with open(dinucSeq.seqfn) as fileObj:
            for i in fileObj:
                head,seq = i.split()
                donor,acceptor = seq[:2].upper(),seq[-2:].upper()
    else:
        dinucSeq = btJuncs.sequence(fi=genome, s=True, tab=True)
        acceptableDinucs = set(['GTAG'])
        if filta == 'gc_at':
            acceptableDinucs.add(['GTAG'])
            acceptableDinucs.add(['ATAC'])

        with open(dinucSeq.seqfn) as fileObj:
            for i in fileObj:
                head,seq = i.split()
                l,strand = head.split("(")
                seqID = l.replace("-",":")
                seqID = seqID.split(":")
                if filta == "all":
                    finalBed += "%s\t%s\t%s\t%s\t%s\t%s\n" % (seqID[0],seqID[1],seqID[2],head,0,strand.rstrip(")"))
                else:
                    donor,acceptor = seq[:2].upper(),seq[-2:].upper()
                    if donor+acceptor in acceptableDinucs:
                        finalBed += "%s\t%s\t%s\t%s\t%s\t%s\n" % (seqID[0],seqID[1],seqID[2],head,0,strand.rstrip(")"))
    return finalBed

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
    
    global minLengthThresh    
    global lengthThresh
    global readThresh

    # filtering
    minLengthThresh = myCommandLine.args["min_length"]
    lengthThresh    = myCommandLine.args["max_length"]
    readThresh      = myCommandLine.args["support_threshold"]
    threads         = myCommandLine.args['threads']
    
    # required args
    fPrefix = myCommandLine.args['out_prefix']
    bedList = myCommandLine.args['bed_manifest']
    genome  = myCommandLine.args['fasta_in']
    ssFilter= myCommandLine.args['filter']
    correct = myCommandLine.args['resolve_strands']


    # initialize tempdir
    with tempfile.TemporaryDirectory() as tmpdirname:
        print(tmpdirname)
        # Collect junction bed files
        beds = juncFilesFromManifest(bedList, tmpdirname)
        # start worker pool
        with Pool(threads) as p:
            for fname in tqdm(p.imap_unordered(runCMD, beds), total=len(beds[:]), position=1, desc='1/2 Filtering sample junctions'):
                continue
        
        #get list of temp files
        flist = [os.path.join(tmpdirname, x) for x in os.listdir(tmpdirname)]
        junctions = set()
        with open('%s_junctions.bed' % fPrefix,'w') as fout:
            for fname in tqdm(flist, total=len(flist),position=1, desc='2/2 Concatenating junctions' ):
                with open(fname,'rb') as fin:
                    for i in fin:
                        if i in junctions:
                            continue
                        else:
                            junctions.add(i)
                            print(i.rstrip().decode(),file=fout)

    # Filter junctions
    filteredBed = filterJunctions('%s_junctions.bed' % fPrefix, ssFilter, correct, genome)
    bTool = pybedtools.BedTool(filteredBed, from_string=True)
    bTool.sort()
    intersection = bTool.intersect(bTool, s=True, wa=True, wb=True)
    bTool.saveas('%s_junctions.bed' % fPrefix)

    
    clusters = makeClusters(intersection)

    with open("%s_all_clusters2.tsv" % fPrefix,'w') as out:
        for i, c in clusters.items():
            print(i,",".join(c), file=out)
    print("done.", file=sys.stderr)


if __name__ == "__main__":
    main()      