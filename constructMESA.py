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


import os, sys
import numpy as np
import pysam
from pybedtools import BedTool
from multiprocessing import Pool

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
											 epilog = 'Please feel free to forward any questions or concerns to /dev/null', 
											 add_help = True, #default is True 
											 prefix_chars = '-', 
											 usage = '%(prog)s -m junction_beds_manifest.tsv ')
		# Add args
		self.parser.add_argument('-m', '--bed_manifest', type=str, action = 'store', required=True, help='List of bed files containing stranded junctions.')
		self.parser.add_argument('-p', '--threads', type=int, action = 'store', default=1, required=False, help='Number of threads.')
		
		
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
			chrom, c1, c2, i, score, strand = line.rstrip().split()[:6]
			if int(score)<5:
				continue
			bedEntry = (chrom, c1, c2,".",0, strand)
			jSet.add(bedEntry)
	#print(jSet)
	return jSet

def runCMD(x):
	'''
	mutiprocessing helper function to run
	multiple instances of bedToCoordSet
	'''
	
	return bedToCoordSet(x)

def makeClusters(intersection):
	'''
	GL000008.2	164884	170271	.	0	-	GL000008.2	164884	170271	.	0	-
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
	threads = myCommandLine.args['threads']
	beds = list()
	try:
		with open(bedList,'r') as lines:
			for line in lines:
				beds.append(line.rstrip().split()[-1])
	except:
		print("No such file or incorrectly formatted file %s" % bedList, file=sys.stderr)
		sys.exit(1)

	p = Pool(threads)
	results = list()

	for i, result in enumerate(p.imap_unordered(runCMD, beds), 1):
		sys.stderr.write('\rReading Bed files ... %s/%s' % (i,len(beds)))#{0:%}'.format(i/len(beds)))
		results.append(result)

	setList = results
	finalSet = set.union(*setList)
	print("\n%s bed files read successfully... %s unique junctions collected." % (len(results),len(finalSet)), file=sys.stderr)

	finalBed = list(finalSet)
	bedObj = BedTool(finalBed).sort()
	
	intersection = bedObj.intersect(bedObj, s=True, wa=True, wb=True)

	bedObj.saveas('all_junctions.bed')
	intersection.saveas('all_junctions_intersection.bed')

	clusters = makeClusters(intersection)

	with open("all_clusters.tsv",'w') as out:
		for i, c in clusters.items():
			print(i,",".join(c), file=out)
	print("done.", file=sys.stderr)


if __name__ == "__main__":
	main()      