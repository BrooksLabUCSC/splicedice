from scipy.stats import beta as stats_beta
import numpy as np

class Manifest:
    def __init__(self,filename):
        self.samples = []
        self.get_group = {}
        self.groups = {}
        with open(filename) as tsv:
            for line in tsv:
                row = line.rstrip().split("\t")
                sample_name,group_name = row[0:2]
                self.samples.append(sample_name)
                self.get_group[sample_name] = group_name
                if group_name not in self.groups:
                    self.groups[group_name] = []
                self.groups[group_name].append(sample_name)

    def get_group_indices(self,samples,anti=False):
        group_indices = {}
        for i,s in enumerate(samples):
            if s in self.get_group:
                try:
                    group_indices[self.get_group[s]].append(i)
                except KeyError:
                    group_indices[self.get_group[s]] = [i]
        if anti:
            anti_indices = {}
            for name,group in group_indices.items():
                anti_indices[name]  = [i for i,s in enumerate(samples) if s not in group]
            return group_indices, anti_indices
        else:
            return group_indices
        
        
#### Beta Class ####        
class Beta:
    def __init__(self,exclude_zero_ones=False):
        self.exclude = exclude_zero_ones
        if self.exclude:
            self.loc = 0
            self.scale = 1
        else:
            self.loc = -0.001
            self.scale = 1.002

    def cdf(self,x,m,a,b):
        if x == m:
            return 1
        else:
            if x > m:
                return 1 - stats_beta.cdf(x,a,b,loc=self.loc,scale=self.scale)
            else:
                return stats_beta.cdf(x,a,b,loc=self.loc,scale=self.scale)

    def fit_beta(self,values):
        values = [x for x in values if not np.isnan(x)]
        if not values:
            return (None,None,(None,None,None))
        if self.exclude:
            values = [x for x in values if x != 0 and x != 1]
        try:
            a,b,l,s = stats_beta.fit(values,floc=self.loc,fscale=self.scale)
        except:
            a,b = None,None 
        return (np.median(values),a,b)
    
#### Multi Class ####        
class Multi:

    @staticmethod
    def mp_reader(read_function,info,q,n): 
        for item in read_function(info):
            q.put(item)

        for i in range(n):
            q.put("DONE")
        return None

    @staticmethod
    def mp_do_rows(q,f,info,o):
        while True:
            item = q.get()
            if item == "DONE":
                o.put("DONE")
                break
            o.put(f(item,info))
        return None

#### Table Class ####        
class Table:
    def __init__(self,filename=None,samples=None,intervals=None,data=None,store=None):
        self.store = store
        if intervals and samples and data:
            self.samples = samples
            self.intervals = intervals
            self.data = data
        else:
            self.filename = filename
            self.samples = None
            self.intervals = None
            self.data = None

    def get_samples(self):
        if self.samples:
            return self.samples
        else:
            with open(self.filename) as tsv:
                return tsv.readline().rstrip().split('\t')[1:]
            
    def get_rows(self,interval_set=None):
        if self.store == None:
            with open(self.filename) as data_file:
                header = data_file.readline().rstrip().split('\t')[1:]
                
                if interval_set == None:
                    for line in data_file:
                        row = line.rstrip().split("\t")
                        yield (row[0],[float(x) for x in row[1:]])
                else:
                    for line in data_file:
                        interval,row = line.rstrip().split("\t",1)
                        if interval in interval_set:
                            yield (interval,[float(x) for x in row.split('\t')])
                        
#### Annotation class ####
class Annotation:
    def __init__(self,gtf_filename):
        self.intervals,self.contigs,self.exons,self.genes = self.getAnnotated(gtf_filename)

    def getAnnotated(self,gtf_filename):
        transcripts = {}
        exons = {}
        with open(gtf_filename) as gtf:
            for line in gtf:
                if line.startswith("#"):
                    continue
                row = line.rstrip().split('\t')
                info = {}
                for x in row[8].split(';'):
                    item = x.split('"')
                    info[item[0].rstrip()] = item[1]
                contig = row[0]
                strand = row[6]
                start = int(row[3])
                stop = int(row[4])-1
                if row[2] == "transcript":
                    transcripts[info['trascript_id']] = [(contig,strand,info['gene_name']),[]]
                elif row[2] == "exon":
                    try:
                        exons[(contig,strand)].add((start,stop))
                    except KeyError:
                        exons[(contig,strand)] = set()
                        exons[(contig,strand)].add((start,stop))
                    transcripts[info['trascript_id']][1].append((start,stop))
        
        for key,values in exons.items():
            exons[key] = sorted(values)
            
        intervals = {}
        contigs = {}
        genes = {}
        for transcript_id,info in transcripts.items():
            contig,strand,gene_name = info[0]
            for i in range(len(info[1])-1):
                left = info[1][i][1]
                right = info[1][i+1][0]
                interval = f"{contig}:{left}-{right}:{strand}"
                contigs[(contig,strand)].add(interval)
                if interval in intervals:
                    intervals[interval].append(transcript_id)
                    if gene_name not in genes[interval]:
                        genes[interval].append[gene_name]
                else:
                    intervals[interval] = [transcript_id]
                    genes[interval] = [gene_name]
        return intervals,contigs,exons,genes
    

