#!/usr/bin/env python3


########################################################################
# File: pairwise_fisher.py
#  executable: 
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 01/08/2020 Created
#
########################################################################

########################################################################
# Hot Imports & Global Variable
########################################################################

import os, sys
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
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
        self.parser = argparse.ArgumentParser(description = 'TBD',
                                             epilog = 'Please feel free to forward any usage questions or concerns', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s --inclusionMESA psis.npz')
        # Add args
        self.parser.add_argument('--inclusionMESA', type=str, action = 'store', required=True, help='Compressed NPZ formatted Inclusion count matrix from quantMESA.') 
        self.parser.add_argument('-c','--clusters', type=str, action = 'store', required=True, help='Clusters table.') 
        self.parser.add_argument('--chi2', action = 'store_true', default=False,  help='Use X^2 instead of fishers. Quicker, not as sensitive.') 
        
        
        if inOpts is None :
            self.args = vars(self.parser.parse_args())
        else :
            self.args = vars(self.parser.parse_args(inOpts))

########################################################################
# Helper Functions
# 
# 
########################################################################

def loadNPZ(x):
    '''
    takes in npz formatted matrix.
    '''

    try:
        data = np.load(x)
    except:
        print("ERR ** Cannot load matrix %s. Check path or format." % x)
        sys.exit(1)
    return data

def getColIndexFromArray(x,y):
    '''
    takes in list of strings = x
    and finds list index in array = y
    '''

    return np.nonzero(np.isin(y,x))
    


def returnSamplesFromManifest(x):
    '''
    reads in mesa formatted manifest
    returns list of samples
    '''
    s = list()
    with open(x) as fin:
        for i in fin:
            s.append(i.split()[0])

    return s

def getClust(fname):
    data = dict()
    with open(fname) as fin:
        for i in fin:
            left,right = i.rstrip().split()
            mxes = right.split(",")
            data[left] = mxes + [left]
    return data
########################################################################
# MAINE
# 
# 
########################################################################

def main():
    '''
    A workflow to compute the significance difference
    between two distributions of PSI values.
    Values are assumed to not be normall distributed, thus
    we invoke the wilcoxon ranksum test as the statistical analysis.
    '''

    myCommandLine = CommandLine()

    # args
    pmesa  = myCommandLine.args["inclusionMESA"]
    cmesa  = myCommandLine.args["clusters"]
    x  = myCommandLine.args["chi2"]
    
    #load psi
    data = loadNPZ(pmesa)

    clusters = getClust(cmesa)

    #table has 3 arrays, cols, rows and data
    cols, rows, matrix = data['cols'], data['rows'], data['data']


    comparisons = set()
    for i,v in enumerate(cols):
        for j,v in enumerate(cols):
            if i == j:
                continue
            comparisons.add(tuple(sorted([i,j])))

    
    # do the math
    comps = list(comparisons)

    print("clusterID","\t".join("%s_%s" % (cols[j[0]],cols[j[1]]) for j in comps), sep="\t")

    pvals = list()
    testedEvents = list()
    if not x:
        for n,vals in enumerate(matrix):
            eventID = rows[n]
            mxes = matrix[np.isin(rows,clusters[eventID])]
            

            inc = vals
            exc = np.sum(mxes,axis=0)
            
            tempPvals = list()
            
            for i in comps:
                left,right = i
                table = [[inc[left],inc[right]],
                         [exc[left],exc[right]]]

                data = fisher_exact(table)[-1]
                tempPvals.append(data)

            
            print(eventID,"\t".join(str(x) for x in tempPvals),sep="\t")
    else:
        for n,vals in enumerate(matrix):
            eventID = rows[n]
            mxes = matrix[np.isin(rows,clusters[eventID])]
            

            inc = vals
            exc = np.sum(mxes,axis=0)
            
            tempPvals = list()
            
            for i in comps:
                left,right = i
                
                table = [[inc[left],inc[right]],
                         [exc[left],exc[right]]]
                try:
                    data = chi2_contingency(table)[1]
                except:
                    # not enough data for text
                    data = np.nan
                tempPvals.append(data)

            
            print(eventID,"\t".join(str(x) for x in tempPvals),sep="\t")


if __name__ == "__main__":
    main()

