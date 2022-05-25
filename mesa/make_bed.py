#!/usr/bin/env python3
"""
Make a bed file from MESA's compare_sample_sets output file
"""

def deltaColor(delta):
    if delta > 0:
        return "0,60,120"
    elif delta < 0:
        return "100,50,0"
    else:
        return "0,0,0"

def makeBed(bedfilename,vsfilename,use_corrected,alpha):
    """Reads comparison file, writes relevant information to bed"""
    with open(bedfilename,'w') as bed:
        track_name = bedfilename.replace('bed','')
        bed.write(f'track name="{track_name}" visibility=2 itemRgb=On\n')
        with open(vsfilename) as vsfile:
            vsfile.readline()
            for line in vsfile:
                row = line.rstrip().split('\t')
                chromosome,coords,strand = row[0].split(':')
                start,end = coords.split('-')
                delta = float(row[4])
                if use_corrected:
                    pval = float(row[7])
                else:
                    pval = float(row[6])
                score = str(int(abs(delta)*999)+1)
                name = f"d:{delta}/p:{pval:.1e}"
                thickstart = start
                thickend = end
                rgb = deltaColor(float(delta))
                if pval <= alpha:
                    bedline = '\t'.join([chromosome,start,end,name,score,strand,thickstart,thickend,rgb])
                    bed.write(f'{bedline}\n')      
                    
def get_args():
    import argparse
    parser = argparse.ArgumentParser(description='Make bed file from output of mesa compare_sample_sets')
    parser.add_argument('--alpha','-a',default=0.05,
                        help='Significance threshold for filtering based on p-value. Default 0.05')
    parser.add_argument('--use_corrected','-c',action='store_true',
                        help='Flag to use corrected p-value in filtering and output label. Default uses raw p-value')
    parser.add_argument('--input','-i',required=True,
                        help='Comparison file from mesa compare_sample_sets')
    parser.add_argument('--output','-o',required=True,
                        help='Output filename')
    return parser.parse_args()

def main():
    """Gets arguments, runs makeBed"""
    args = get_args()
    alpha = args.alpha
    use_corrected = args.use_corrected
    vsfilename = args.input
    bedfilename = args.output
           
    makeBed(bedfilename,vsfilename,use_corrected,alpha)    

if __name__ == "__main__":
    main()