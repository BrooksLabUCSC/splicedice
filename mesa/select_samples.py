"""
MESA select
"""

def getJunctions(filename1,filename2):
    with open(filename1) as file1:
        with open(filename2) as file2:

            header1 = file1.readline().rstrip().split("\t")
            header2 = file2.readline().rstrip().split("\t")


            junctions1 = set()
            for line in file1:
                row = line.split("\t")
                junctions1.add(row[0])

            junctions2 = set()
            for line in file2:
                row = line.split("\t")
                junctions2.add(row[0])
    return junctions1,junctions2

def writeOutput(filename1,filename2,junctions,outfilename):
    with open(outfilename,"w") as output:
        with open(filename1) as file1:
            with open(filename2) as file2:

                header1 = file1.readline().rstrip().split("\t")
                header2 = file2.readline().rstrip().split("\t")

                print(f"Writing table with {len(header1)+len(header2)-2} samples and {len(junctions)} junctions")
                header = "\t".join(header1+header2[1:])
                output.write(f"{header}\n")

                to_write = {}
                for line in file1:
                    row = line.rstrip().split("\t")
                    if row[0] in junctions:
                        to_write[row[0]] = row[1:]

                for line in file2:
                    row = line.rstrip().split("\t")
                    if row[0] in junctions:
                        new_line = "\t".join(row[0:1] + to_write[row[0]] + row[1:])
                        output.write(f"{new_line}\n")



def add_parser(parser):
    """ """
    parser.add_argument("--allps1","-a1",
                        action="store",required=True,
                       help="First allPS file output from mesa quant")
    parser.add_argument("--allps2","-a2",
                        action="store",required=True,
                       help="Second allPS file output from mesa quant")
    parser.add_argument("--output","-o",
                        action="store",required=True,
                       help="Output filename")
    parser.add_argument("--junctions","-j",
                        choices=["intersection","union","first","second"],
                        default="intersection",
                        help="Which junctions to select from each allPS table (Default AND CURRENTLY ONLY OPTION: intersection)")
    
def run_with(args):
    """ """
    junctions1,junctions2 = getJunctions(args.allps1,args.allps2)
    
    junctions = junctions1.intersection(junctions2)
    
    writeOutput(args.allps1,args.allps2,junctions,args.output)

    
    
if __name__ == "__main__":
    import argparse 
    parser = argparse.ArgumentParser(description='Combines two allPS tables into one file.')           
    add_parser(parser)
    args = parser.parse_args()
    run_with(args)
