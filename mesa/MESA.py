#!/usr/bin/env python3

"""
Mutually Exclusive Splicing Analysis (MESA) main quantification step

"""

import numpy as np
from time import time

class Sample:
    
    sampleList = []
    groups = {}
    
    def __init__(self,manifestLine):
        
        self.name = manifestLine[0]
        self.filename = manifestLine[1]
        Sample.sampleList.append(self)
        
        # Check filetype
        if self.filename.upper().endswith(".BED"):
            self.type = "bed"
        elif self.filename.upper().endswith("SJ.OUT.TAB"):
            self.type = "SJ"
        elif self.filename.upper().endswith(".BAM"):
            self.type = "bam"
        else:
            self.type = "unknown"
            
        

        self.metadata = manifestLine[2]
        self.condition = manifestLine[3]

        if self.condition in Sample.groups:
            Sample.groups[self.condition].append(self)
        else:
            Sample.groups[self.condition] = [self]
            
class Timer:
    def __init__(self):
        self.start = time()
        self.checkpoint = self.start
        
    def total(self):
        passed = time() - self.start
        hours = int(passed // 3600)
        minutes = int((passed % 3600) // 60)
        seconds = passed % 60
        return f"[{hours}:{minutes:02d}:{seconds:02.2f}]"
            
    def check(self):
        passed = time() - self.checkpoint
        self.checkpoint = time()
        hours = int(passed // 3600)
        minutes = int((passed % 3600) // 60)
        seconds = passed % 60
        return f"[{hours}:{minutes:02d}:{seconds:02.2f}]"
        
        
                        
        
class MESA:
    """Main algorithm for Mutually Exclusive Splicing Analysis"""
    def __init__(self,manifestFilename,outputPrefix,args):
        """Call methods to get input, process, and write output"""
        # Parsed Arguments: 
        self.args = args
        
        #
        self.manifestFilename = manifestFilename
        self.outputPrefix = outputPrefix
                
        # Construct MESA
        timer = Timer()
        print("Parsing manifest...")
        self.manifest = self.parseManifest()
        print("\tDone",timer.check())
        
        print(f"Getting all junctions from {len(self.manifest)} files...")
        #print(f"Getting all junctions from {len(self.manifest)} files...",end=' ')
        self.junctions = self.getAllJunctions()
        print("\tDone",timer.check())

        
        print(f"Finding clusters from {len(self.junctions)} junctions...")
        self.clusters = self.getClusters()
        self.junctionIndex = {junction:i for i,junction in enumerate(sorted(self.clusters))}
        print("\tDone",timer.check())
        
        # Write tsv file
        print("Writing cluster file...")
        self.writeClusters()
        print("\tDone",timer.check())
        
        # Write bed file
        print("Writing junction bed file...")
        self.writeJunctionBed()
        print("\tDone",timer.check())
        
        # Quantify MESA
        print("Gathering junction counts...")
        self.counts, self.low = self.getJunctionCounts()
        print("\tDone",timer.check())
        
        print("Writing inclusion counts...")
        self.writeInclusions()
        print("\tDone",timer.check())
        
        print("Calculating PS values...")
        self.psi = self.calculatePsi()
        print("\tDone",timer.check())
        
        print("Writing PS values...")
        self.writeAllpsi()
        print("\tDone",timer.check())
        
        if self.args.drim:
            print("Writing drim table...")
            self.writeDrimTable()
            print("\tDone",timer.check())
        
        print("All done",timer.total())
        
        
    def parseManifest(self):
        """Get sample info and paths from manifest"""
        manifest = []
        with open(self.manifestFilename,"r") as manifestFile:
            for line in manifestFile:
                row = line.rstrip().split("\t")
                if len(row) != 4:
                    pass # improperly formatted manifest
                sample = Sample(row)
                manifest.append(sample)
        return manifest
                

    def getAllJunctions(self):
        """Read junctions from SJ_out_tab"""
        
        strandSymbol = {"0":"0", "1":"+", "2":"-", "+":"+", "-":"-"}
        plusminus = {"+","-"}
        
        filters = {"gtag_only":{1,2}, "gc_at":{1,2,3,4,5}, "all":{0,1,2,3,4,5,6}}
        validMotifs = filters[self.args.filter]
                
        junctions = set()
        
        # Read all sample files from manifest
        for sample in self.manifest:
            with open(sample.filename,"r") as junctionFile:
                
                if sample.type is "SJ":
                    for line in junctionFile:
                        row = line.rstrip().split("\t")
                        chromosome= row[0]
                        left = int(row[1]) - 1
                        right = int(row[2])
                        strand = strandSymbol[row[3]]
                        intronMotif = int(row[4])
                        #annotation = int(row[5])
                        #overhang = int(row[8])
                        if self.args.noMultimap:
                            score = int(row[6])
                        else:
                            score = int(row[6]) + int(row[7]) 
                        if (right-left < self.args.maxLength and 
                            right-left > self.args.minLength and
                            strand != "0" and
                            score >= self.args.minUnique and
                            intronMotif in validMotifs):
                            junctions.add((chromosome,left,right,strand))
                            
                elif sample.type is "bed":
                    for line in junctionFile:
                        row = line.rstrip().split("\t")

                        score = int(row[4])
                        if score < self.args.minUnique:
                            continue
                             
                        left = int(row[1])
                        right = int(row[2])
                        length = right-left
                        if length > self.args.maxLength or length < self.args.minLength:    
                            continue
                            
                        strand = row[5]
                        if strand in plusminus:
                            chromosome = row[0]
                            junctions.add((chromosome,left,right,strand))
                            
        return junctions
        
    def getClusters(self):
        """Read all junctions from *.SJ.out.tab file  """
        
        chromosome = None
        strand = None
        clusters = {}
        
        for junction in sorted(self.junctions, key = lambda x: (x[0],x[3],x[1],x[2])):
            
            # Reset 
            if junction[0] != chromosome or junction[3] != strand:
                chromosome = junction[0]
                strand = junction[3]
                potentialOverlaps = []
                
            clusters[junction] = []
                
            newPotentialOverlaps = [junction]
            
            for priorJunction in potentialOverlaps:
                if priorJunction[2] >= junction[1]: 
                    clusters[priorJunction].append(junction)
                    clusters[junction].append(priorJunction)  
                    newPotentialOverlaps.append(priorJunction)
            potentialOverlaps = newPotentialOverlaps
        return clusters
            
    def getJunctionCounts(self):
        """ """
        counts = np.zeros((len(self.clusters),len(self.manifest)),dtype='float32')
        low = []
        
        for sampleIndex,sample in enumerate(self.manifest):
            
            with open(sample.filename,"r") as sampleFile:
                
                if sample.type is "bed":
                    for line in sampleFile:
                        row = line.rstrip().split("\t")

                        junction = (row[0], int(row[1]), int(row[2]), row[5])
                            
                        if junction in self.junctionIndex:
                            score = int(row[4])
                            counts[self.junctionIndex[junction],sampleIndex] = score
                            if self.args.lowCoverageNan and score < self.args.minUnique:
                                low.append((self.junctionIndex[junction],sampleIndex))
                    
                elif sample.type is "SJ":
                                        
                    strandSymbol = {'0':'0', '1':'+', '2':'-'}
                    
                    for line in sampleFile:
                        row = line.rstrip().split("\t")

                        junction = (row[0], int(row[1])-1, int(row[2]), strandSymbol[row[3]])

                                                    
                        if junction in self.junctionIndex:
                            if self.args.noMultimap:
                                counts[self.junctionIndex[junction],sampleIndex] = int(row[6])
                            else:
                                counts[self.junctionIndex[junction],sampleIndex] = int(row[6]) + int(row[7]) 
                          

        return counts, low
                        
    def calculatePsi(self):
        """ """
        psi = np.zeros((len(self.clusters),len(self.manifest)),dtype='float32')
        for junction in sorted(self.clusters):
            chromosome,left,right,strand = junction
            inclusions = self.counts[self.junctionIndex[junction],:]
            exclusions = np.zeros(len(self.manifest))
            for excluded in self.clusters[junction]:
                exclusions += self.counts[self.junctionIndex[excluded],:]
            psi[self.junctionIndex[junction],:] = inclusions / (inclusions + exclusions)
        if self.args.lowCoverageNan:
            for junctionIndex,sampleIndex in self.low:
                psi[junctionIndex,sampleIndex] = np.nan          
        return psi

    def junctionString(self,junction):
        """ """
        return f"{junction[0]}:{junction[1]}-{junction[2]}"
        
    def writeJunctionBed(self):
        with open(f"{self.outputPrefix}_junctions.bed", "w") as outbed:
            for junction in sorted(self.junctions):
                chromosome,left,right,strand = junction
                name = f"{chromosome}:{left}-{right}({strand})"
                outbed.write(f"{chromosome}\t{left}\t{right}\t{name}\t0\t{strand}\n")
            
        
    def writeClusters(self):
        """"""
        with open(f"{self.outputPrefix}_allClusters.tsv","w") as clusterFile:
            for junction in sorted(self.clusters):
                line = f"{junction[0]}:{junction[1]}-{junction[2]}\t"
                line += ",".join([f"{j[0]}:{j[1]}-{j[2]}" for j in self.clusters[junction]])
                print(line, file=clusterFile)
                
    def writeInclusions(self):
        """ """
        tab = '\t'
        with open(f"{self.outputPrefix}_inclusionCounts.tsv","w") as inclusionTsv:
            
            inclusionTsv.write("cluster\t")
            inclusionTsv.write("\t".join([s.name for s in self.manifest])+"\n")
            
            for i,junction in enumerate(sorted(self.clusters)):
                inclusionTsv.write(f"{self.junctionString(junction)}\t{tab.join([f'{x:.0f}' for x in self.counts[i,:]])}\n")
                
                
    def writeAllpsi(self):
        """ """
        tab = '\t'
        with open(f"{self.outputPrefix}_allPS.tsv","w") as allpsTsv:
            
            samples = "\t".join([s.name for s in self.manifest])
            allpsTsv.write(f"cluster\t{samples}\n")

            
            for i,junction in enumerate(sorted(self.clusters)):
                allpsTsv.write(f"{self.junctionString(junction)}\t{tab.join([f'{x:.3f}' for x in self.psi[i,:]])}\n")
                
    def writeDrimLine(self,i,junction,other,file):
        """Format and output line for drim table"""
        print(f"cl_{i}_{self.junctionString(junction)}",
              f"{self.junctionString(other)}_{i}",
              "\t".join(self.counts[self.junctionIndex[other],:].astype("str")),
              sep="\t", file=file)
        
    def writeDrimTable(self):
        """Write tsv file to use for DRIMSeq"""
        with open(f"{self.outputPrefix}_drimTable.tsv", "w") as drimTable:
            sampleNames = "\t".join([s.name for s in self.manifest])
            drimTable.write(f"gene\tfeature_id\t{sampleNames}\n")
            for i,junction in enumerate(sorted(self.clusters)):
                self.writeDrimLine(i,junction,junction,drimTable)
                for excludedJunction in self.clusters[junction]:
                    self.writeDrimLine(i,junction,excludedJunction,drimTable)
        
                
                
                

        
def add_parser(parser):
    """ """
    parser.add_argument("--manifest","-m",
                        action="store",required=True,
                       help="tab-separated list of samples with file paths")
    parser.add_argument("--output_prefix","-o",
                        action="store",required=True,
                       help="prefix for output filenames") 
    parser.add_argument("--maxLength",type=int,default=50000,
                       help="maximum splice junction size")
    parser.add_argument("--minLength",type=int,default=50,
                       help="minimum splice junction size")
    parser.add_argument("--minOverhang",type=int,default=5,
                       help="minimum overlap on reads to support splice junction")
    parser.add_argument("--drim",action="store_true",
                       help="create table for use by DRIMSeq")
    parser.add_argument("--noMultimap",action="store_true",
                       help="use only reads that uniquely map to one location")
    parser.add_argument("--filter",default="gtag_only",choices=["gtag_only"],
                       help="donor and acceptor intron sequences to include.")
    parser.add_argument("--minUnique",type=int,default=5,
                        help="minimum number of unique reads to support splice junction")
    parser.add_argument("--lowCoverageNan",action="store_true",
                        help="Report NaN for splicing events with coverage below minUnique")
    
def run_with(args):
    """ Main program which calls MESA algorithm class"""
    manifestFilename = args.manifest
    outputPrefix = args.output_prefix

    MESA(manifestFilename,outputPrefix,args)

if __name__ == "__main__":
    import argparse 
    parser = argparse.ArgumentParser(description='Mutually Exclusive Splicing Analysis.')           
    add_parser(parser)
    args = parser.parse_args()
    run_with(args)



