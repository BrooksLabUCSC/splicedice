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

########################################################################
# Hot Imports & Global Variable
########################################################################

import time
import os, sys
from tqdm import *
from pybedtools import BedTool
from multiprocessing import Pool, Manager
from hurry.filesize import size

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
        self.parser = argparse.ArgumentParser(description = 'constructMESA.py - Construyo una mesa',
                                             #epilog = 'Please feel free to forward any questions or concerns to /dev/null', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s -m junction_beds_manifest.tsv ')
        # Add args
        self.parser.add_argument('-m', '--bed_manifest', type=str, action = 'store', required=True, help='List of bed files containing stranded junctions.')
        self.parser.add_argument('-p', '--threads', type=int, action = 'store', default=1, required=False, help='Number of threads.')
        self.parser.add_argument('-f', '--fasta_in', action = 'store', required=False, help='BedTools indexed Genome Fasta')
        self.parser.add_argument('--junction_support_threshold', action = 'store', type=int, default = 10, required=False, help='Junction = Real if read support > Minimum threshold (10)')
        self.parser.add_argument('--junction_length_max', action = 'store', type=int, default = 500000, required=False, help='Junction = Real if length < length threshold (500000)')
        self.parser.add_argument('--junction_length_min', action = 'store', type=int, default = 50, required=False, help='Junction = Real if length > length threshold (50)')
        self.parser.add_argument('--junction_min_samp_threshold', action = 'store', type=float, default = 0.1, required=False, help='Junction = Real if read support > Minimum threshold (0.1)')
        
        
        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))

########################################################################
# Helper Functions
# 
# 
########################################################################

def bedToCoordSet(f):
    '''
    takes bed file and returns a set of junctions
    '''
    jSet = set()

    with open(f,'r') as lines:
        for line in lines:
            cols = line.rstrip().split()
            chrom, c1, c2, i, score, strand = cols[:6]
            c1,c2,score = map(int, [c1,c2,score])
            if score<readThresh or c2 - c1 > lengthThresh:
                continue
            
            bedEntry = (chrom, c1, c2, ".", 0, strand)
            jSet.add(bedEntry)
            
    #print(f,sys.getsizeof(jSet))
    
    return jSet#s

def runCMD(x):
    '''
    mutiprocessing helper function to run
    multiple instances of bedToCoordSet
    '''
    return bedToCoordSet(x)

def makeClusters(intersection):
    '''
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

def juncFilesFromManifest(manifest):
    '''function takes in 2 column manifest
    and returns list of file names (second column)'''

    beds = list()
    try:
        with open(manifest,'r') as lines:
                for line in lines:
                    sample, file = line.rstrip().split()
                    beds.append(file)
    except:
        print("** ERR: %s expected 2 columns in manifest. Check format. **" % manifest, file=sys.stderr)
        sys.exit(1)
    return beds

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

    lengthThresh = myCommandLine.args["junction_length_max"]
    minLengthThresh = myCommandLine.args["junction_length_min"]
    readThresh = myCommandLine.args["junction_support_threshold"]
    bedList = myCommandLine.args['bed_manifest']
    threads = myCommandLine.args['threads']
    genome = myCommandLine.args['fasta_in']
    


    # Collect junction bed files
    beds = juncFilesFromManifest(bedList)
    p = Pool(threads)

    finalJunctionSet = set()
    l = list()
    for result in tqdm(p.imap_unordered(runCMD, beds[:]), total=len(beds[:])):
        j = result
        l.extend(j)
        if sys.getsizeof(l, set()) > 2500000000:
            finalJunctionSet = finalJunctionSet | set(l)
            l.clear()



    finalJunctionSet = finalJunctionSet | set(l)
    print("\n%s bed files read successfully... %s unique junctions collected." % (len(beds),len(finalJunctionSet)), file=sys.stderr)


    finalJunctions = list(finalJunctionSet)
    bedObj = BedTool(finalJunctions).sort()
    
    intersection = bedObj.intersect(bedObj, s=True, wa=True, wb=True)

    bedObj.saveas('all_junctions2.bed')
    intersection.saveas('all_junctions_intersection2.bed')

    clusters = makeClusters(intersection)

    with open("all_clusters2.tsv",'w') as out:
        for i, c in clusters.items():
            print(i,",".join(c), file=out)
    print("done.", file=sys.stderr)


if __name__ == "__main__":
    main()      