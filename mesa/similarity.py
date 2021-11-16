"""
MESA similarity
"""

def readVsFile(vs_filename):
    with open(vs_filename) as vs_file:
        vs_file.readline()

        midpoints = {}
        deltas = {}
        for line in vs_file:
            row = line.rstrip().split("\t")
            event = row[0]
            median1 = float(row[3])
            delta = float(row[5])
            midpoints[event] = median1-(delta/2)
            deltas[event] = delta
    return midpoints,deltas
        
    
def scoreSamples(allps_filename,midpoints,deltas):    
    with open(allps_filename) as allps_file:
        samples = allps_file.readline().rstrip().split("\t")[1:]
        scores = [0 for sample in samples]
        count = 0
        for line in allps_file:
            row = line.rstrip().split("\t")
            event = row[0]
            if event not in deltas or deltas[event] == 0:
                continue
            count += 1
            if deltas[event] < 0:
                for i,ps in enumerate(row[1:]):
                    if float(ps) < midpoints[event]:
                        scores[i] += 1
            elif deltas[event] > 0:
                for i,ps in enumerate(row[1:]):
                    if float(ps) > midpoints[event]:
                        scores[i] += 1
    return samples,scores,count
        
def getGroups(manifest_filename):
    groups = {}
    with open(manifest_filename) as manifest:
        for line in manifest:
            name,path,group1,group2 = line.rstrip().split("\t")
            groups[name] = (group1,group2)
    return groups

def writeScores(output_filename,scores,samples,count,groups=None):
    with open(output_filename,"w") as score_file:
        for score,sample in sorted(zip(scores,samples)):
            if groups:
                group1,group2 = groups[sample]
                score_file.write(f"{sample}\t{group1}\t{group2}\t{score/count:0.03f}\n")
            else:
                score_file.write(f"{sample}\t{score/count:0.03f}\n")
        


def add_parser(parser):
    """ """
    parser.add_argument("--manifest","-m",
                        action="store",default=None,
                       help="tab-separated list of samples for group names")
    parser.add_argument("--comparison","-c",
                        action="store",required=True,
                       help="Output table from compare_sample_sets")
    parser.add_argument("--allps","-a",
                        action="store",required=True,
                       help="Allps table from mesa quant")
    parser.add_argument("--output","-o",
                        action="store",required=True,
                       help="Output filename")
    
def run_with(args):
    """ """
    vs_filename = args.comparison
    allps_filename = args.allps
    manifest_filename = args.manifest
    output_filename = args.output
        
    midpoints,deltas = readVsFile(vs_filename)
    samples,scores,count = scoreSamples(allps_filename,midpoints,deltas)
    
    if manifest_filename:
        groups = getGroups(manifest_filename)
        writeScores(output_filename,scores,samples,count,groups)
    else:
        writeScores(output_filename,scores,samples,count)


if __name__ == "__main__":
    import argparse 
    parser = argparse.ArgumentParser(description='Scores similarity of samples to a condition from mesa compare_sample_sets.')           
    add_parser(parser)
    args = parser.parse_args()
    run_with(args)
