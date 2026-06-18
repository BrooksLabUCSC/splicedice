#!/usr/bin/env python3

"""
Main quantification step

Precondition: All junctions reported to this script have passed necessary filters.

Expected input bed formats:
1. chromosome
2. left coordinate (0-based half open)
3. right coordinate (0-based half open)
4. name (not used here)
5. score (used as read count for quantification)
6. strand (+ or -)
"""

import numpy as np
from time import time

class Sample:
    
    sampleList = []
    
    def __init__(self,manifestLine):
        self.name = manifestLine[0]
        self.filename = manifestLine[1]
        Sample.sampleList.append(self)
            
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

        
class SPLICEDICE:
    """Main algorithm for Mutually Exclusive Splicing Analysis"""
    def __init__(self,manifestFilename,outputPrefix,args):
        """Call methods to get input, process, and write output"""
        # Parsed Arguments: 
        self.args = args
        
        self.manifestFilename = manifestFilename
        self.outputPrefix = outputPrefix
                
        # Construct SPLICEDICE
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
        
        # Quantify SPLICEDICE
        print("Gathering junction counts...")
        self.counts = self.getJunctionCounts()

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
        """
        Build a union of junctions from all samples. 
        Iterate the input manifest and read each input .bed file.

        BED files are expected to have the following column format:
        1. chromosome
        2. left coordinate (0-based half open)
        3. right coordinate (0-based half open)
        4. Any (not used here)
        5. Any (not used here)
        6. strand (+ or -)

        Only junctions with a valid strand (+ or -) are included.

        Returns:
            A set of tuples (chromosome, left, right, strand) representing all junctions observed in the input files.
        """
    
        plusminus = {"+","-"}
        junctions = set()
        # Read all sample files from manifest
        for sample in self.manifest:
            with open(sample.filename,"r") as junctionFile:
                for line in junctionFile:
                    row = line.rstrip().split("\t")                    
                    chrom = row[0]
                    start = int(row[1])
                    end = int(row[2])
                    strand = row[5]
                    if strand in plusminus:
                        junctions.add((chrom, start, end, strand))
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
        """
        Build a matrix of junction counts. This matrix allows any sample to report a read count for 
        all junctions in the union of each sample's junctions. Junctions not observed in a sample will have a count of 0.

        Returns:
            A 2D numpy array of shape (number of junctions across all samples, number of samples) where 
            each entry [i, j] contains the read count for junction i in sample j.
        """
        counts = np.zeros((len(self.clusters),len(self.manifest)),dtype='float32')
        for sampleIndex,sample in enumerate(self.manifest):
            with open(sample.filename,"r") as sampleFile:
                for line in sampleFile:
                    row = line.rstrip().split("\t")
                    junction = (row[0], int(row[1]), int(row[2]), row[5])
                    if junction in self.junctionIndex:
                        score = int(row[4])
                        counts[self.junctionIndex[junction],sampleIndex] = score
        return counts
                        
    def calculatePsi(self):
        """ """
        psi = np.zeros((len(self.clusters),len(self.manifest)),dtype='float32')
        for junction in sorted(self.clusters):
            chromosome,left,right,strand = junction
            inclusions = self.counts[self.junctionIndex[junction],:]
            exclusions = np.zeros(len(self.manifest))
            for excluded in self.clusters[junction]:
                exclusions += self.counts[self.junctionIndex[excluded],:]
            denom = inclusions + exclusions
            psi[self.junctionIndex[junction],:] = np.where(
                denom > 0,
                np.divide(inclusions, denom, where=denom > 0),
                np.nan
            )
        return psi

    def junctionString(self, junction, one_based=False):
        """
        Convert a junction tuple to a string representation.

        Args:
            junction (tuple): A tuple containing (chromosome, left, right, strand).
            one_based (bool): If True, convert coordinates to 1-based. Default is False (0-based half open).
        """        
        chromosome = junction[0]
        start = junction[1] + 1 if one_based else junction[1]
        end = junction[2]
        strand = junction[3]
        return f"{chromosome}:{start}-{end}:{strand}"
        
    def writeJunctionBed(self):
        with open(f"{self.outputPrefix}_junctions.bed", "w") as outbed:
            for junction in sorted(self.junctions):
                chromosome,left,right,strand = junction
                name = self.junctionString(junction, True)
                outbed.write(f"{chromosome}\t{left}\t{right}\t{name}\t0\t{strand}\n")
            
    def writeClusters(self):
        """"""
        with open(f"{self.outputPrefix}_allClusters.tsv","w") as clusterFile:
            for junction in sorted(self.clusters):
                line = f"{self.junctionString(junction, True)}\t"
                line += ",".join([f"{self.junctionString(j, True)}" for j in self.clusters[junction]])
                print(line, file=clusterFile)
                
    def writeInclusions(self):
        """ """
        tab = '\t'
        with open(f"{self.outputPrefix}_inclusionCounts.tsv","w") as inclusionTsv:
            inclusionTsv.write(f"cluster\t{tab.join([s.name for s in self.manifest])}\n")

            
            for i,junction in enumerate(sorted(self.clusters)):
                inclusionTsv.write(f"{self.junctionString(junction, True)}\t{tab.join([f'{x:.0f}' for x in self.counts[i,:]])}\n")
                
                
    def writeAllpsi(self):
        """ """
        tab = '\t'
        with open(f"{self.outputPrefix}_allPS.tsv","w") as allpsTsv:
            
            samples = "\t".join([s.name for s in self.manifest])
            allpsTsv.write(f"cluster\t{samples}\n")

            
            for i,junction in enumerate(sorted(self.clusters)):
                allpsTsv.write(f"{self.junctionString(junction, True)}\t{tab.join([f'{x:.3f}' for x in self.psi[i,:]])}\n")
                
    def writeDrimLine(self,i,junction,other,file):
        """Format and output line for drim table"""
        print(f"cl_{i}_{self.junctionString(junction, True)}",
              f"{self.junctionString(other, True)}_{i}",
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
    parser.add_argument("--drim",action="store_true",
                       help="create table for use by DRIMSeq")


def run_with(args):
    """ Main program which calls SPLICEDICE algorithm class"""
    manifestFilename = args.manifest
    outputPrefix = args.output_prefix

    SPLICEDICE(manifestFilename,outputPrefix,args)


if __name__ == "__main__":
    import argparse 
    parser = argparse.ArgumentParser(description='Mutually Exclusive Splicing Analysis.')           
    add_parser(parser)
    args = parser.parse_args()
    run_with(args)



