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


import sys
import numpy as np
from sklearn.decomposition import PCA
from sklearn import manifold
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import umap
import numpy.ma as ma
import pandas as pd
from scipy.cluster.hierarchy import linkage


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
                    j = j.replace("NA", "nan")
                cols = j.split()
                tempList.append(cols[:-1])

        returnList.append(np.asarray(tempList, dtype=float))

    shapes = set([x.shape[0] for x in returnList])
    if len(shapes) > 1:
        print("something went wrong with number of events", file=sys.stderr)
        sys.exit(1)

    mergedMatrix = np.concatenate(returnList, axis=1)
    groupIndices = list()
    counter = 0
    for i in returnList:
        temp = list()
        for num, j in enumerate(i[0, :]):
            temp.append(num + counter)
        counter = counter + len(temp)
        groupIndices.append(temp)
    return mergedMatrix, groupIndices, samples


def makePCA(dm, groups, colorD, markerD):

    X = PCA(n_components=2)
    X_pca = X.fit_transform(dm)
    X_variance = X.explained_variance_ratio_
    # first and second component
    for k, v in groups.items():
        v = np.array(v)
        indices, label = np.array(v[:, 0], dtype=int), v[:, 1]
        plt.scatter(
            X_pca[indices, 0],
            X_pca[indices, 1],
            label=label[0],
            color=colorD[k[0]],
            marker=markerD[k[1]],
            alpha=0.2,
        )

    plt.xlabel("PCA 1 (%.2f variance)" % X_variance[0])
    plt.ylabel("PCA 2 (%.2f variance)" % X_variance[1])
    plt.legend()
    plt.tight_layout()
    plt.savefig("mesaPCA.pdf")


def makeUMAP(dm, groups, colorD, markerD):
    # sub ploting
    (fig, subplots) = plt.subplots(1, 3, figsize=(16, 4))
    # number of components for PCA
    components = [25, 50, 100]
    # components = [1,2,4] # for samps smaller than 100

    for num, s in enumerate(components):
        X_pca = PCA(n_components=s).fit_transform(dm)
        reducer = umap.UMAP()
        embedding = reducer.fit_transform(X_pca)
        for k, v in groups.items():
            v = np.array(v)
            indices, label = np.array(v[:, 0], dtype=int), v[:, 1]
            subplots[num].set_title("UMAP Comp=%d" % (s))
            subplots[num].scatter(
                embedding[indices, 0],
                embedding[indices, 1],
                label=label[0],
                color=colorD[k[0]],
                marker=markerD[k[1]],
                alpha=0.2,
            )
    plt.legend()
    plt.tight_layout()
    plt.savefig("mesaUMAP.pdf")


def makeTSNE(dm, groups, colorD, markerD):

    # sub ploting
    (fig, subplots) = plt.subplots(4, 3, figsize=(12, 8))
    # number of components for PCA
    components = [25, 50, 100]
    # components = [1,2,4] # for samps smaller than 100
    # num of perplexities
    perplexities = [5, 30, 50, 100]

    for num1, perp in enumerate(perplexities):
        tsne = manifold.TSNE(
            n_components=2, init="random", random_state=0, perplexity=perp
        )

        for num2, c in enumerate(components):
            X_pca = PCA(n_components=c).fit_transform(dm)
            Y = tsne.fit_transform(X_pca)
            subplots[num1][num2].set_title("Perp=%d Comp=%d" % (perp, c))

            for k, v in groups.items():
                v = np.array(v)
                indices, label = np.array(v[:, 0], dtype=int), v[:, 1]
                subplots[num1][num2].scatter(
                    Y[indices, 0],
                    Y[indices, 1],
                    label=label[0],
                    color=colorD[k[0]],
                    marker=markerD[k[1]],
                    alpha=0.2,
                )
    plt.legend()
    plt.tight_layout()
    plt.savefig("mesaTSNE.pdf")


# def makeSHC(dm, groups,groupIndices, colorD, makrerD, samps):
def makeSHC(dm, groups, colorD, markerD, samps):

    df = pd.DataFrame(data=dm, index=samps)
    # sns.clustermap(df)
    Z = linkage(df, "ward")

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


def getColorMarkers(data, cols):

    # samps plots will be a dict
    # where key is marker / color
    # and value is list of samples and their column index in
    # the data matrix
    sampsPlots = dict()

    # group1 first
    g1 = np.unique(data[:, 0])
    if len(g1) <= 10:
        colors = plt.get_cmap("tab10")(np.arange(10, dtype=int))
        # colors = cmap(np.linspace(0, 1, 10))
    elif len(g1) <= 20:
        cmap = plt.get_cmap("tab20")(np.arange(20, dtype=int))
        # colors = cmap(np.linspace(0, 1, 20))
    else:
        colors = plt.get_cmap("viridis")(np.arange(len(g1), dtype=int))
        # colors = cmap(np.linspace(0, 1, len(g1)))

    # colors = ['red','green']

    # now group2 ( i expect less unique values in this group)
    g2 = np.unique(data[:, 1])
    markers = [".", "+", "^", "v", "1", "2", "3", "4", "*"]

    for i in data:
        val1, val2, sampID = i
        val1Ind = tuple(np.argwhere(g1 == val1)[0])[0]
        val2Ind = tuple(np.argwhere(g2 == val2)[0])[0]
        matrixColumnInd = np.argwhere(cols == sampID)[0][0]
        key = (val1Ind, val2Ind)
        if key not in sampsPlots:
            sampsPlots[key] = list()
        sampsPlots[key].append((matrixColumnInd, "%s_%s" % (g1[val1Ind], g2[val2Ind])))

    return colors, markers, sampsPlots


def loadNPZ(x):
    """
    takes in npz formatted matrix.
    """

    try:
        data = np.load(x)
    except:
        print("ERR ** Cannot load matrix %s. Check path or format." % x)
        sys.exit(1)
    return data


########################################################################
# Main
########################################################################

def add_parser(parser):
    parser.add_argument(
        "--psiMESA",
        type=str,
        action="store",
        required=False,
        help="Compressed NPZ formatted PSI matrix from quantMESA.",
    )
    parser.add_argument("-m", "--manifest", required=True,)


def run_with(args):
    pmesa = args.psiMESA
    manifest = args.manifest

    data = loadNPZ(pmesa)
    cols, rows, matrix = data["cols"], data["rows"], data["data"]
    mergedMatrix, samps = matrix, cols
    print("Pre filtering shape", mergedMatrix.shape)

    # filter events and fill in missing data, otherwise it's all the same?
    N = 0.5
    mergedMatrix = mergedMatrix[
        np.sum(~np.isnan(mergedMatrix), axis=1) / len(mergedMatrix[0, :]) > N
    ]

    # filter events with no variability
    N = 0.01
    finalM = mergedMatrix[np.where(np.nanstd(mergedMatrix, axis=1) > N)]

    # Fill in missing data with the average, not sure how this works
    finalM = np.where(
        np.isnan(finalM),
        ma.array(finalM, mask=np.isnan(finalM)).mean(axis=1)[:, np.newaxis],
        finalM,
    )
    print("Post filtering shape", finalM.shape)
    dm = np.transpose(finalM)

    # get covariates
    # formatting
    groupData = list()
    with open(manifest) as fin:
        for i in fin:
            info = i.rstrip().split()
            group1, group2 = info[-2], info[-1]
            groupData.append((group1, group2, info[0]))
    groupData = np.array(groupData)
    colors, markers, sampsPlots = getColorMarkers(groupData, samps)

    # plotting
    # heirch - in dev
    # makeSHC(dm,sampsPlots,colors,markers, groupData[:,-1])
    # PCA
    makePCA(dm, sampsPlots, colors, markers)

    if len(groupData) < 100:
        sys.exit()

    # #UMAP
    makeUMAP(dm, sampsPlots, colors, markers)

    # #TSNE
    makeTSNE(dm, sampsPlots, colors, markers)
