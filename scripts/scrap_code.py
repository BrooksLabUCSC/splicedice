
# For handling percent-spliced value tables
# Written by Dennis Mulligan



import numpy as np


class Samples:

    def __init__(self,samples=None,groups=None,manifest=None):

        if samples and groups:
            self.samples = samples
            self.groups = groups
        elif samples:
            self.samples = samples
            self.groups = {}
        elif manifest:
            self.samples,self.groups = self.parse_manifest()
        else:
            pass

        self.index = {sample:i for i,sample in enumerate(self.samples)}
        self.m = len(self.samples)

    def parse_manifest(self,manifest):
        """Read sample manifest file. Single column (name) or tab separated (name,group,*)."""
        with open(manifest) as tsv:
            row = tsv.readline().rstrip().split("\t")
            if len(row) == 1:
                samples = [row[0]]
                groups = {}
                for line in tsv:
                    row = line.rstrip().split("\t")
                    samples.append(row[0])
            else:
                samples = [row[0]]
                groups = {row[1]:[row[0]]}
                for line in tsv:
                    row = line.rstrip().split("\t")
                    samples.append(row[0])
                    try:
                        groups[row[1]].append(row[0])
                    except KeyError:
                        groups[row[1]] = [row[0]]
        return samples,groups
    
    def combine(self,other):
        for sample in other.samples:
            if sample in self.index:
                raise ValueError('Overlapping sample names, cannot combine.')
        new_groups = self.groups.items() ^ other.groups.items()
        for group in self.groups.keys() & other.groups.keys():
            new_groups[group] = self.groups[group] + other.groups[group]
        new_samples = self.samples + other.samples
        return Samples(samples=new_samples,groups=new_groups)
    

            


class Intervals:

    @staticmethod
    def parse_interval(string):
        name,span,strand = string.split(":")
        left,right = (int(s) for s in span.split("-"))
        return (name,left,right,strand)

    @staticmethod
    def interval_string(name,left,right,strand):
        return f"{name}:{left}-{right}:{strand}"
    
    def __init__(self,intervals=None,contigs=None):

        if intervals and contigs:
            self.intervals = intervals
            self.contigs = contigs
        elif intervals:
            self.intervals = intervals
            self.contigs = self.get_contigs(intervals)
        elif contigs:
            self.contigs = contigs
            self.intervals = self.get_intervals(contigs)

        self.index = {interval:i for i,interval in enumerate(self.intervals)}
        self.n = len(self.intervals)

        self.exclusions = {}

    def get_contigs(self,intervals):
        contigs = {}
        for i,interval in enumerate(intervals):
            name,left,right,strand = self.parse_interval(interval)
            try:
                contigs[(name,strand)].append((left,right,i))
            except KeyError:
                contigs[(name,strand)] = [(left,right,i)]
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


class Table:
    def __init__(self,filename=None,intervals=None,samples=None,data=None):
        if filename:
            self.samples,self.intervals,self.data = self.read_table(filename)
        elif intervals and samples and data:
            self.samples = samples
            self.intervals = intervals
            self.data = data
        else:
            self.samples = []
            self.intervals = []
            self.data = {}

    def read_table(self,table_file):
        with open(table_file) as table:
            samples = table.readline().rstrip().split("\t")[1:]
            intervals = []
            data = []
            for line in table:
                row = line.rstrip().split("\t")
                intervals.append(row[0])
                data.append([float(x) for x in row[1:]])
            intervals = Intervals(intervals)
        return samples,intervals,data

    def write_table(self,out_filename):
        with open(out_filename,'w') as tsv:
            tab = "\t"
            tsv.write(f"splice_interval\t{tab.join(self.samples.samples)}\n")
            for interval,row in zip(self.intervals.intervals,self.data):
                tsv.write(f"{interval}\t{tab.join(str(val) for val in row)}\n")
        return None
    
    def read_table_filter(self,table_file,intervals):
        interval_set = set(intervals)
        intervals = []
        with open(table_file) as table:
            samples = table.readline().rstrip().split("\t")[1:]
            data = []
            for line in table:
                row = line.rstrip().split("\t")
                if row[0] in interval_set:
                    intervals.append(row[0])
                    data.append([float(x) for x in row[1:]])
        return samples,intervals,data
    
    def get_sample(self,name):
        return [self.data[i][self.sample_index[name]] for i in range(self.intervals.n)]
    
    def combine(self,other,keep="both"):
        if keep == "union":
            intervals = self.intervals.union(other.intervals)
        elif keep == "intersection":
            intervals = self.intervals.intersection(other.intervals)
        elif keep == "left":
            intervals = self.intervals
        elif keep == "right":
            intervals = other.intervals
        samples = self.samples.combine(other.samples)
        data = [[] for i in range(intervals.n)]
        for i,interval in enumerate(intervals):
            data[i] = self.get_row(interval) + other.get_row(interval)
        return Table(intervals=intervals,samples=samples,data=data)
    
    def reset_na(self,splice_interval,counts):
        
        name,left,right,strand = Intervals.parse_interval(splice_interval)
        contig = (name,strand)
        span_i = (left,right,self.intervals.index[splice_interval])

        ex_count = np.zeros(range(self.samples.m))
        for exclusion in self.exclusions[contig][span_i]:
            ex_count += counts[exclusion]
            splice_interval =  self.intervals.intervals[span_i[2]]

 
            



    def get_row(self,interval,get_exclusion=False):
        try:
            return self.data[self.intervals.index[interval]]
        except KeyError:
            if get_exclusion:

                return ["nan" for i in range(self.samples.m)]
            else:
                return ["nan" for i in range(self.samples.m)]
        
        

        



