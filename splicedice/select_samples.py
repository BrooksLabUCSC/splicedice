"""
SPLICEDICE select
"""

def parseManifest(manifest_filename):
    with open(manifest_filename) as manifest:
        filenames = []
        for line in manifest:
            filenames.append(line.rstrip())
    return filenames
            
def getJunctions(filenames):
    """ """
    all_junctions = {}
    for filename in filenames:
        with open(filename) as file:
            junctions = set()
            file.readline()
            for line in file:
                junctions.add(line.split("\t")[0])
        all_junctions[filename] = junctions
    return all_junctions
    


def decideJunctions(all_junctions,which):
    if which == "union":
        return sorted(set.union(*all_junctions.values()))
    if which == "intersection":
        return sorted(set.intersection(*all_junctions.values()))
    

def writeOutput(filenames,junctions,outfilename):
    
    header = {}
    table_rows = {j:{} for j in junctions}
    for filename in filenames:
        with open(filename) as file:
            header[filename] = file.readline().rstrip().split("\t")[1:]
            for line in file:
                row = line.rstrip().split("\t")
                if row[0] in junctions:
                    table_rows[row[0]][filename] = row[1:]
                        
    with open(outfilename,"w") as output:
        output.write("cluster\t")
        
        output.write("\t".join(["\t".join(header[f]) for f in filenames]) + "\n")
        
        for junction in junctions:
            table_row = table_rows[junction]
            out_row = [junction]
            for filename in filenames:
                if filename in table_row:
                    out_row.extend(table_row[filename])
                else:
                    out_row.extend(["nan"]*len(header[filename]))
            output.write("\t".join(out_row)+"\n")
        



def add_parser(parser):
    """ """
    parser.add_argument("--allps1","-a1",
                        action="store",
                       help="First allPS file output from splicedice quant")
    parser.add_argument("--allps2","-a2",
                        action="store",
                       help="Second allPS file output from splicedice quant")
    parser.add_argument("--manifest","-m",
                        action="store",
                        help="File with list of allPS file paths (supercedes allPS1 and allPS2)")
    parser.add_argument("--output","-o",
                        action="store",required=True,
                       help="Output filename")
    parser.add_argument("--join","-j",
                        choices=["intersection","union"],
                        default="intersection",
                        help="Which junctions to select from each allPS table (Default AND CURRENTLY ONLY OPTION: intersection)")
    
def run_with(args):
    """ """
    if args.manifest:
        filenames = parseManifest(args.manifest)
    elif args.allps1 and args.allps2:
        filenames = [args.allps1,args.allps2]
    else:
        import sys
        sys.exit("Requires --manifest file or --allPS1 and --allPS2 files.")
                 
    all_junctions = getJunctions(filenames)
    
    junctions = decideJunctions(all_junctions,args.join)
    
    writeOutput(filenames,junctions,args.output)

    
    
if __name__ == "__main__":
    import argparse 
    parser = argparse.ArgumentParser(description='Combines two allPS tables into one file.')           
    add_parser(parser)
    args = parser.parse_args()
    run_with(args)
