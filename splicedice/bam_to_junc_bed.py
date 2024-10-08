#!/usr/bin/env python3

"""
Read BAM files to extract all splice junctions

Takes manifest of BAM false
Outputs bed files with 

"""

import os    
import pysam
from multiprocessing import Pool
import numpy as np
    
class BamToJuncBed:
    
    def __init__(self,args):
        self.args = args
        
        
        self.manifest_filename = args.manifest
        self.output_prefix = args.output_prefix
        self.n_threads = args.number_threads

        self.bed_directory = os.path.join(os.getcwd(),f"{self.output_prefix}_junction_beds")
        try:
            os.mkdir(self.bed_directory)
        except FileExistsError:
            pass

        bams, new_manifest = self.parseManifest()
        
        self.annotated = self.getAnnotated()

        print(f"Finding junctions from {len(bams)} BAM files...")

        self.bamsToBeds(bams)

        self.writeNewManifest(new_manifest)
        
    def parseManifest(self):
        """
        """
        bams = []
        new_manifest = []
        with open(self.manifest_filename, "r") as manifestFile:
            for line in manifestFile:
                samplename, filename, metadata, condition = line.rstrip().split("\t")

                if filename.endswith(".bam") or filename.endswith(".sam"):
                    bedfilename = os.path.join(self.bed_directory,f"{os.path.basename(filename)[:-4]}.junc.bed")
                    bams.append([samplename,filename,metadata,condition,bedfilename])
                    new_manifest.append(f"{samplename}\t{bedfilename}\t{metadata}\t{condition}")
                else:
                    new_manifest.append(line)

        return bams, new_manifest
    
    def getAnnotated(self):
        genes = {}
        names = {}
        transcripts = {}
        with open(self.args.annotation) as gtf:
            for line in gtf:
                if line.startswith("#"):
                    continue
                row = line.rstrip().split('\t')
                if row[2] == "transcript":
                    info = [x.split('"') for x in row[8].split(';')]
                    tid =  [x[1] for x in info if 'transcript_id' in x[0]][0]
                    try:
                        gid = [x[1] for x in info if 'gene_name' in x[0]][0]
                    except IndexError:
                        gid = [x[1] for x in info if 'gene_id' in x[0]][0]
                    genes[tid] = gid
                    transcripts[(tid,row[0],row[6])] = []
                if row[2] == "exon":
                    info = [x.split('"') for x in row[8].split(';')]
                    tid =  [x[1] for x in info if 'transcript_id' in x[0]][0]
                    transcripts[(tid,row[0],row[6])].append((int(row[3]),int(row[4])))
        annotated = {}
        for transcript,exons in transcripts.items():
            tid,chromosome,strand = transcript
            for i in range(len(exons)-1):
                annotated[(chromosome,exons[i][1],exons[i+1][0]-1,strand)] = genes[tid]
        return annotated

    def getJunctionsFromBam(self,sample):
        """
        """
        min_length = self.args.min_length
        max_length = self.args.max_length
        min_reads = self.args.min_reads
        fasta = self.args.genome
        
        
        samplename, filename, metadata, condition, bedfilename = sample
        
        #genome = pysam.AlignmentFile(filename)
        #############
        old_verbosity = pysam.set_verbosity(0)
        try:
            genome = pysam.AlignmentFile(filename)
        except ValueError:
            print("Using: pysam.AlignmentFile(filename,check_sq=False) with",filename)
            genome = pysam.AlignmentFile(filename,check_sq=False)
        pysam.set_verbosity(old_verbosity)
        ##############
        counts = {}
        leftDiversity = {}
        rightDiversity = {}
        overhangs = {}
        for read in genome.fetch(until_eof=True):
            if True: #read.is_read2:
                if read.is_reverse:
                    strand = "-"
                else:
                    strand = "+"
            else:
                if read.is_reverse:
                    strand = "+"
                else:
                    strand = "-"

            blocks = read.get_blocks()
            try:
                read_start = blocks[0][0]
            except IndexError:
                continue
            read_end = blocks[-1][1]
            for i in range(len(blocks)-1):
                junction = (read.reference_name,blocks[i][1],blocks[i+1][0],strand)
                length = junction[2] - junction[1]
                if length >= min_length and length <= max_length:
                    leftOH = blocks[i][1]-blocks[i][0]
                    rightOH = blocks[i+1][1]-blocks[i+1][0]
                    overhang = min(leftOH,rightOH)
                    try:
                        counts[junction] += 1
                        overhangs[junction] = max(overhang,overhangs[junction])
                        try:
                            leftDiversity[junction][read_start] += 1
                            rightDiversity[junction][read_end] += 1
                        except KeyError:
                            leftDiversity[junction][read_start] = 1
                            rightDiversity[junction][read_end] = 1
                    except KeyError:
                        counts[junction] = 1
                        overhangs[junction] = overhang
                        leftDiversity[junction] = {read_start:1}
                        rightDiversity[junction] = {read_end:1}

        filteredJunctions = []
        leftEntropy = {}
        rightEntropy = {}

        if genome:
            leftMotif = {}
            rightMotif = {}
            genome = pysam.FastaFile(fasta)
        for junction in sorted(counts):
            chromosome,left,right,strand = junction
            
            if genome:
                if (chromosome,left) not in leftMotif:
                    try:
                        leftMotif[(chromosome,left)] = genome.fetch(chromosome,left,left+2)
                    except KeyError:
                        leftMotif[(chromosome,left)] = "NN"
                if (chromosome,right) not in rightMotif:
                    try:
                        rightMotif[(chromosome,right)] = genome.fetch(chromosome,right-2,right)
                    except KeyError:
                        rightMotif[(chromosome,right)] = "NN"
            leftEntropy[junction] = 0
            total = sum(leftDiversity[junction].values())
            for species,count in leftDiversity[junction].items():
                prop = count/total
                leftEntropy[junction] -= (prop) * np.log(prop)
            rightEntropy[junction] = 0
            total = sum(rightDiversity[junction].values())
            for species,count in rightDiversity[junction].items():
                prop = count/total
                rightEntropy[junction] -= (prop) * np.log(prop)

            filteredJunctions.append(junction)

            
        #
        if self.args.strands in ("inferCombine", "inferOnly"):
            
            opposite = {"+":"-", "-":"+"}
            plus_motifs = {"GT_AG","GC_AG","AT_AC"}
            minus_motifs = {"CT_AC","CT_GC","GT_AT"}
            
            firstFiltered = filteredJunctions
            filteredJunctions = []
                        
            for junction in firstFiltered:
                
                chromosome,left,right,strand = junction
                motif = f"{leftMotif[(chromosome,left)]}_{rightMotif[(chromosome,right)]}"


                complement = (chromosome,left,right,opposite[strand])
                
                if complement not in counts:
                    filteredJunctions.append(junction)
                elif (junction in self.annotated or
                     (strand == "+" and motif in plus_motifs) or
                     (strand == "-" and motif in minus_motifs)):

                    filteredJunctions.append(junction)
                    
                    if self.args.strands == "inferCombine":
                        counts[junction] += counts[complement]
                        
                elif (complement in self.annotated or
                     (strand == "-" and motif in plus_motifs) or
                     (strand == "+" and motif in minus_motifs)):
                    pass
                else:
                    filteredJunctions.append(junction)
                        


        
        with open(bedfilename,"w") as bedOut:
            for junction in filteredJunctions:
                chromosome,left,right,strand = junction
                name = f"e:{leftEntropy[junction]:0.02f}:{rightEntropy[junction]:0.02f};o:{overhangs[junction]};m:{leftMotif[(chromosome,left)]}_{rightMotif[(chromosome,right)]};a:{self.annotated.get(junction,'?')}"
                bedOut.write(f"{chromosome}\t{left}\t{right}\t{name}\t{counts[junction]}\t{strand}\n")

        return filename, bedfilename, len(filteredJunctions)


    def bamsToBeds(self,bams):
        """
        """
        with Pool(self.n_threads) as pool:
            for info in pool.imap(self.getJunctionsFromBam,bams):
                filename,bedfilename,num_junctions = info
                print("bam:", filename)
                print("number of junctions found:",num_junctions)
                print("saved to bed:", bedfilename)


    def writeNewManifest(self,new_manifest):
        """
        """
        #newManifestPath = os.path.join(path,f"{outputPrefix}_manifest.txt")
        new_manifest_path = f"{self.output_prefix}_manifest.txt"
        with open(new_manifest_path,"w") as manifest_file:
            for line in new_manifest:
                manifest_file.write(line+"\n")
        print("new manifest written to:", new_manifest_path)


def add_parser(parser):
    """ """
    parser.add_argument("--manifest","-m",required=True,
                       help="Tab-separated file of samples and bam file paths.")
    parser.add_argument("--output_prefix","-o",required=True,
                       help="Prefix for junction file directory and manifest.")
    parser.add_argument("--genome","-g",
                        help="Optional. Genome fasta file to report splice site motifs")
    parser.add_argument("--annotation","-a",
                        help="Optional. Gene annotation gtf file to report annotation status")
    parser.add_argument("--max_length",type=int,default=100000,
                       help="Maximum distance between ends of junction [Default 100000]")
    parser.add_argument("--min_length",type=int,default=50,
                       help="Minimum distance between ends of junction [Default 50]")
    parser.add_argument("--min_reads",type=int,default=5,
                       help="Minimum number of reads required to report junction [Default 5]")
    
#    parser.add_argument("--no_multimap",action="store_true",
#                       help="STILL IN PROGRESS")

    parser.add_argument("--number_threads","-n",type=int,default=1,
                       help="Number of bam files to search concurrently [Default 1]")
    
    parser.add_argument("--strands","-s",default="inferCombine",
                        choices=["keepBoth","inferOnly","inferCombine"],
                        help="How to handle junctions with same coordinates on opposite strands. [Default 'inferCombine'. Options 'keepBoth','inferOnly']")

    
def run_with(args):
    """
    """
    
    BamToJuncBed(args)
        

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    add_parser(parser)
    args = parser.parse_args()
    run_with(args)
