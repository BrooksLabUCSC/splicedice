#!/usr/bin/env python3


########################################################################
# File: clusterMESA.py
#  executable: 
# Purpose: 
#
#          
# Author: Cameron M. Soulette
# History:      cms 12/05/2019 Created
#
########################################################################

########################################################################
# Hot Imports & Global Variable
########################################################################


import os, sys
import numpy as np
from sklearn.decomposition import PCA
from sklearn import manifold
from time import time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import umap
import numpy.ma as ma
import scipy.cluster.hierarchy as shc
import pandas as pd
import seaborn as sns

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
		self.parser.add_argument('--quant_mesa', action = 'store', required=True, 
								help='all PSI quant matrix')

		self.parser.add_argument('-m', '--manifest' action = 'store', required=True, 
								help='Sample manifest')
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
# Functions
########################################################################

def quantsToMatrix(f):

	samples = list()
	returnList = list()
	for i in f:
		print(i)
		tempList = list()
		with open(i.rstrip()) as mat:
			samples = np.asarray(next(mat).rstrip().split()[:-1])
			for j in mat:
				# for JB
				if "NA" in j:
					j = j.replace("NA",'nan')
				cols = j.split()
				tempList.append(cols[:-1])
				
		returnList.append(np.asarray(tempList,dtype=float))

	shapes = set([x.shape[0] for x in returnList])
	if len(shapes)>1:
		print("something went wrong with number of events",file=sys.stderr)
		sys.exit(1)

	mergedMatrix = np.concatenate(returnList, axis=1)
	groupIndices = list()
	counter = 0
	for i in returnList:
		temp = list()
		for num,j in enumerate(i[0,:]):
			temp.append(num+counter)
		counter = counter + len(temp)
		groupIndices.append(temp)
	return mergedMatrix, groupIndices, samples
	

def makePCA(dm,groups,groupIndices,colorD,markerD):
	
	X_pca = PCA(n_components=50).fit_transform(dm)
	#first and second component
	for num,g in enumerate(groups):
		indices = groupIndices[num]
		plt.scatter(X_pca[indices,0],X_pca[indices,1],label="%s_%s" % (g[0],g[1]), color=colorD[g[0]], marker=markerD[g[1]], alpha=0.2)

	plt.legend()
	plt.tight_layout()
	plt.savefig("mesaPCA.pdf")


def makeUMAP(dm, groups, groupIndices, colorD, markerD):

	# sub ploting
	(fig, subplots) = plt.subplots(1, 3, figsize=(16,4))
	# number of components for PCA
	components = [25,50,100]

	for num,s in enumerate(components):
		X_pca = PCA(n_components=s).fit_transform(dm)
		reducer = umap.UMAP()
		embedding = reducer.fit_transform(X_pca)
		for num2,g in enumerate(groups):
			indices = groupIndices[num2]
			subplots[num].set_title("UMAP Comp=%d" % (s))
			subplots[num].scatter(embedding[indices,0],embedding[indices,1],label="%s_%s" % (g[0],g[1]), color=colorD[g[0]], marker=markerD[g[1]], alpha=0.2)
	plt.legend()
	plt.tight_layout()
	plt.savefig("mesaUMAP.pdf")
	
def makeTSNE(dm, groups, groupIndices, colorD, markerD):

	# sub ploting
	(fig, subplots) = plt.subplots(4, 3, figsize=(12,8))
	# number of components for PCA
	components = [25,50,100]
	# num of perplexities
	perplexities = [5,30, 50, 100]

	for num1,perp in enumerate(perplexities):
		tsne = manifold.TSNE(n_components=2, init='random',
                         random_state=0, perplexity=perp)
		
		for num2,c in enumerate(components):
			X_pca = PCA(n_components=c).fit_transform(dm)
			Y = tsne.fit_transform(X_pca)
			subplots[num1][num2].set_title("Perp=%d Comp=%d" % (perp,c))
		
			for num3,g in enumerate(groups):
				indices = groupIndices[num3]
				subplots[num1][num2].scatter(Y[indices,0],Y[indices,1],label="%s_%s" % (g[0],g[1]), color=colorD[g[0]], marker=markerD[g[1]], alpha=0.2)
	plt.legend()
	plt.tight_layout()
	plt.savefig("mesaTSNE.pdf")
	

def makeSHC(dm, groups,groupIndices, colorD, makrerD, samps):


	df = pd.DataFrame(data=dm, index=samps)
	sns.clustermap(df)
	# plt.figure(figsize=(10, 7))
	# plt.title("Customer Dendograms")
	# linked = shc.linkage(dm, method='ward')
	# labelList = samps
	# dend = shc.dendrogram(linked, labels=labelList)
	# # dend = shc.dendrogram(shc.linkage(dm, method='ward'))
	# # labels = [int(x) for x in dend['ivl']]
	# # # labels = [item.get_text() for item in plt.xticklabels()]
	# # # labels[1] = 'Testing'
	# # plt.xticks(labels, samps[[int(x) for x in labels]])
	plt.legend()
	plt.tight_layout()
	plt.savefig("mesaSHC.pdf")
	

def getColorMarkers(g1,g2):
	cmap = plt.get_cmap('viridis')
	colors = cmap(np.linspace(0, 1, len(g1)))
	colorDict = dict()
	for num,i in enumerate(g1):
		colorDict[i] = colors[num]

	markers=['.', '+', '^', 'v', '1', '2','3','4','*']
	markerDict = dict()
	for num,i in enumerate(g2):
		markerDict[i] = markers[num]

	return colorDict, markerDict



########################################################################
# Main
########################################################################

def main():

	'''
	TDB
	'''

	myCommandLine = CommandLine()


	files = myCommandLine.args["quant_files"]

	# formatting
	groups = list()
	g1 = set()
	g2 = set()
	flist = list()
	with open(files) as fin:
		for i in fin:
			info = i.rstrip().split()
			group1,group2 = info[-2],info[-1]
			groups.append((group1,group2))
			g1.add(group1)
			g2.add(group2)
			flist.append(info[0])

	g1,g2 = sorted(list(g1)), sorted(list(g2))

	mergedMatrix,groupIndices,samps = quantsToMatrix(flist)
	#print(mergedMatrix.shape,groupIndices)
	
	# filter events and fill in missing data, otherwise it's all the same?
	N = 0.5 
	mergedMatrix = mergedMatrix[np.sum(~np.isnan(mergedMatrix),axis=1) / len(mergedMatrix[0,:]) > N]
	
	#filter events with no variability
	N = 0.01
	finalM =  mergedMatrix[ np.where(np.nanstd(mergedMatrix,axis=1) > N) ]
	
	# Fill in missing data with the average, not sure how this works
	finalM = np.where(np.isnan(finalM), ma.array(finalM, mask=np.isnan(finalM)).mean(axis=1)[:, np.newaxis], finalM)
	
	dm = np.transpose(finalM)
	print(dm.shape)
	# plotting

	colorDict,markerDict = getColorMarkers(g1,g2)
	#heirch
	#makeSHC(dm,groups,groupIndices,colorDict,markerDict,samps)
	#PCA
	makePCA(dm,groups,groupIndices,colorDict,markerDict)
	# #UMAP
	makeUMAP(dm,groups,groupIndices,colorDict,markerDict)
	# #TSNE
	makeTSNE(dm,groups,groupIndices,colorDict,markerDict)


if __name__ == "__main__":
	main()