#!/usr/bin/env python3


########################################################################
# File: splitMESA.py
#  executable: splitMESA.py
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 11/05/2018 Created
#
########################################################################

########################################################################
# Hot Imports & Global Variable
########################################################################


import os, sys
from multiprocessing import Pool
import numpy as np

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
        self.parser = argparse.ArgumentParser(description = 'splitMESA.py - A tool to split MESA into sample groups.',
                                             epilog = 'Please feel free to forward any questions or concerns to /dev/null', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s -m sample_manifest.tsv -i mesa_table.tsv -d output_dir')
        # Add args
        self.parser.add_argument('-i', '--input_mesa', action = 'store', required=True, help='MESA Table from quantMESA.')
        self.parser.add_argument('-m', '--manifest', type=str, action = 'store', required=True, help='List of bed files containing stranded junctions.')
        self.parser.add_argument('-d', '--working_dir', type=str, action = 'store', required=True, help='Working Dir.')
        
        
        
        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))

########################################################################
# Helper Functions
# 
# 
########################################################################

def readMESA(f):

    with f as lines:

        header = next(lines)
        headerDict = {val:pos for pos,val in header.split()}
        yield headerDict

        for l in lines:
            cols = l.split()
            yield [x.split(":") for x in cols[:-2]] + [cols[-2],cols[-1]] 


def splitData(events,mani,k):

    tempF = dict()
    manifestRoot = "_manifest_temp.tsv"
    sampleByGroup = dict()
    with open(mani) as lines:
        for line in lines:
            c = line.rstrip().split("\t")
            sample, juncs, group1, group2 = c
            if group1 not in tempF:
                sampleByGroup[group1] = set()
                tempF[group1] = open(k + "/" + group1 + manifestRoot,'w')

            print(line.rstrip(), file=tempF[group1])
            sampleByGroup[group1].add(c[0])
    [x.close() for y,x in tempF.items()]

    tempF = dict()
    eventsRoot = "_events_temp.tsv"
    with open(events) as lines:
        header = next(lines).rstrip().split()
        groups = {y:[pos for pos,x in enumerate(header) if x in sampleByGroup[y]] for y in sampleByGroup}
        tempF  = {y: open(k + "/" + y + eventsRoot,'w') for y in groups}
        [print("\t".join(header[i] for i in groups[x]), "intron","mxe", sep="\t", file=y) for x,y in tempF.items()]
        
        for line in lines:
            c = line.rstrip().split()
            [print("\t".join(c[i] for i in groups[x]), c[-2], c[-1], sep="\t", file=y) for x,y in tempF.items()]

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
    
  
    inputMesa = myCommandLine.args["input_mesa"]
    manifest = myCommandLine.args["manifest"]
    bedList = myCommandLine.args['threads']
    wdir = myCommandLine.args['working_dir']
    
    try:
        os.mkdir(wdir)
    except:
        pass

    files = splitData(inputMesa,manifest,wdir)


if __name__ == "__main__":
    main()      