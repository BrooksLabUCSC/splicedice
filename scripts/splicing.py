# For analyzing alternative splicing events
# Written by Dennis Mulligan
import numpy as np

from tools import Annotation


class AlternativeSplicing:

    def __init__(self,intervals,annotation):
        self.alt5,self.alt3,self.skips = self.classify(intervals,annotation)

    def classify(self,intervals,annotation):
        alt5,alt3 = [],[]
        skips = []
        for contig,spans in intervals.contigs.items():
            if contig[1] == "+":
                five,three = 0,1
            elif contig[1] == "-":
                five,three = 1,0
            exon_spans = annotation.exons[contig]
            i = 0
            n = len(exon_spans)
            for span in spans:
                if (contig,span[five]) in alt5:
                    alt5[(contig,span[five])].append(span[three])
                else:
                    alt5[(contig,span[five])] = [span[three]]
                if (contig,span[three]) in alt3:
                    alt3[(contig,span[three])].append(span[five])
                else:
                    alt3[(contig,span[three])] = [span[five]]
                while i < n and exon_spans[i][0] < span[0]:
                    i += 1
                if exon_spans[i][1] <= span[1]:
                    skips.append(contig,span)
                else:
                    while exon_spans[i-1][0] > span[0] and i>0:
                        i -= 1
        return alt5,alt3,skips
            
    #### Intervals class

class Intervals(list):

    @staticmethod
    def parse_interval(string):
        name,span,strand = string.split(":")
        left,right = (int(s) for s in span.split("-"))
        return (name,left,right,strand)

    @staticmethod
    def get_string(name,left,right,strand):
        return f"{name}:{left}-{right}:{strand}"
    
    def __init__(self,intervals=[],contigs={}):
        if contigs:
            self.contigs = contigs
            if not intervals:
                intervals = self.get_intervals(contigs)
        list.__init__(self,intervals)
        if not contigs:
            self.contigs = self.get_contigs
        # Store properties
        self.index = {interval:i for i,interval in enumerate(self)}
        self.n = len(self)
        self.exclusions = {}

    def get_contigs(self,intervals):
        contigs = {}
        for i,interval in enumerate(intervals):
            contig,left,right,strand = self.parse_interval(interval)
            try:
                contigs[(contig,strand)].append((left,right))
            except KeyError:
                contigs[(contig,strand)] = [(left,right)]
        for span_list in contigs.values():
            span_list.sort()
        return contigs

    def get_intervals(self,contigs):
        intervals = []
        for contig,spans in contigs.items():
            name,strand = contig
            for left,right in spans:
                intervals.append(f"{name}:{left}-{right}:{strand}")
        return intervals


    def read_exclusion(self,exclusion_file):
        exclusions = {}
        with open(exclusion_file) as tsv:
            for line in tsv:
                interval,exclusive = line.rstrip().split("\t")
                exclusions[interval] = exclusive.split(",") # <<<--------- CHECK THIS CHARACTER
        return exclusions

    def find_exclusion(self,contigs):
        exclusions = {}
        for contig,span_list in contigs.items():
            span_list.sort()
            exclusions[contig] = {}
            active = []
            for new_span in span_list:
                left = new_span[0]
                exclusions[contig][new_span] = []
                new_active = []
                for span in active:
                    if left <= span[1]:
                        exclusions[contig][span].append(new_span)
                        exclusions[contig][new_span].append(span)
                        new_active.append(span)
                new_active.append(new_span)
                active = new_active
                exclusions[contig][new_span].append(new_span)
        return exclusions
    
    def combine(self,other,keep="union"):
        if keep == "union":
            contigs = self.contigs.items() ^ other.contigs.items()
            for contig in self.contigs.keys():
                if contig in other.contigs:
                    contigs[contig] = sorted(set(self.contigs[contig]).intersection(set(other.contigs[contig])))
        elif keep == "intersection":
            contigs = {}
            for contig in self.contigs.keys():
                if contig in other.contigs:
                    contigs[contig] = sorted(set(self.contigs[contig]).intersection(set(other.contigs[contig])))
        return Intervals(contigs=contigs)
    
    def intersection(self,other):
        contigs = {}
        for contig in self.contigs.keys():
            if contig in other.contigs:
                contigs[contig] = sorted(set(self.contigs[contig]).intersection(set(other.contigs[contig])))
        return Intervals(contigs=contigs)
    
    def union(self,other):
        contigs = self.contigs.items() ^ other.contigs.items()
        for contig in self.contigs.keys():
            if contig in other.contigs:
                contigs[contig] = sorted(set(self.contigs[contig]).intersection(set(other.contigs[contig])))
        return Intervals(contigs=contigs)



def cluster_intervals(intervals):
    clusters = {}
    active_contig = None
    for contig,left,right,strand in intervals:
        if contig != active_contig:
            active_contig = contig
            clusters[contig] = {s:[] for s in "+-"}
            active_cluster = {s:[] for s in "+-"}
        if not active_cluster[strand] or left <= active_cluster[strand][-1][1]:
            active_cluster[strand].append((left,right))
        else:
            clusters[contig][strand].append(active_cluster[strand])
            active_cluster[strand] = [(left,right)]
    for strand in "+-":
        if active_cluster[strand]:
            clusters[contig][strand].append(active_cluster[strand])
    return clusters

def cluster_spans(spans):
    clusters = []
    active_cluster = []
    for left,right in spans:
        active_cluster = {s:[] for s in "+-"}
        if not active_cluster or left <= active_cluster[-1][1]:
            active_cluster.append((left,right))
        else:
            clusters.append(active_cluster)
            active_cluster = [(left,right)]
    if active_cluster:
        clusters.append(active_cluster)
    return clusters

def find_exclusion(cluster):
    n = len(cluster)
    exclusions = [[] for i in range(n)]
    for i in range(n):
        right = cluster[i][1]
        exclusions[i].append(i)
        j = i+1
        while j < n and cluster[j][0] <= right:
            exclusions[i].append(j)
            exclusions[j].append(i)
            j += 1
    return exclusions

def calculate_ps(exclusions,counts):
    for i,ex_index in enumerate(exclusions):
        counts[i,:] / np.sum(counts[ex_index,:])

def read_counts(filename):
    with open(filename) as tsv:
        header = tsv.readline()
        for line in tsv:
            row = line.rstrip().split("\t")
            yield row[0],np.array(row[1:],dtype=int)

def something(counts):
    for interval,row in counts:
        pass

def find_skipping(self,contigs,exons):
        skips = []
        for contig,spans in contigs.items():
            exon_spans = exons[contig]
            i = 0
            n = len(exon_spans)
            for span in spans:
                while i < n and exon_spans[i][0] < span[0]:
                    i += 1
                if exon_spans[i][1] <= span[1]:
                    skips.append(contig,span)
                else:
                    while exon_spans[i-1][0] > span[0] and i>0:
                        i -= 1

def binary_search(query,intervals):
    low = 0
    high = len(intervals)
    i = (low+high)//2
    old_i = None
    while i != old_i:
        if query[0] < intervals[i][0]:
            low = i
        elif query[1] > intervals[i][1]:
            high = i
        else:
            return True
        old_i = i
        i = (low+high)//2
    return False



#### ####               #### ####
    
def main():

    pass

if __name__=="__main__":
    main()

