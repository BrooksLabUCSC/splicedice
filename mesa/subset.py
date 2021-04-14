#!/usr/bin/env python3
"""
For subsetting a percent-spliced table output from mesa (filename ending with allPS.tsv).
"""


def add_parser(parser):
    """
    """

    parser.add_argument('-j', '--junctionSubset', 
                        action = 'store', required=False, default=None,
                        help='Text file of junctions to include in subset, one per line [Default includes all junctions]')
    parser.add_argument('-a', '--spliceSites', 
                        action = 'store', required=False, default=None,
                        help='Text file of alternative splice sites to include in subset, one per line []')
    parser.add_argument('-s', '--sampleSubset', 
                        action = 'store', required=False, default=None,
                        help='Text file of samples to include in subset, one per line [Default includes all samples]')
    parser.add_argument('-i', '--inputFilename',
                        required=True,
                        help='Table of PS values from MESA quant (*_allPS.tsv)')
    parser.add_argument('-o', '--outputFilename',
                        required=True,
                        help='Table of PS values for selected samples and junctions')
    parser.add_argument('-f','--psRangeFilter',
                        type=float,default=0.1,
                        help='Do not include junctions that have a difference in min and max PS values below cutoff [Default 0.1]')


def run_with(args):
    """ """
    input_filename = args.inputFilename
    output_filename = args.outputFilename
    
    if args.sampleSubset:
        sample_subset = set()
        with open(args.sampleSubset) as ssfile:
            for line in ssfile:
                sample_subset.add(line.strip())
        all_samples = False
    else:
        all_samples = True
        
    if args.junctionSubset:
        junction_subset = set()
        with open(args.junctionSubset) as jsfile:
            for line in jsfile:
                #junction = line.strip().split('-')
                #junction[1] = str(int(junction[1])-1)
                #junction_subset.add('-'.join(junction))
                junction_subset.add(line.strip())
        all_junctions = False
        find_splice_sites = False

    elif args.spliceSites:
        splice_site_subset = set()
        with open(args.spliceSites) as siteFile:
            for line in siteFile:
                chromosome,a = line.strip().split(':')
                splice_site_subset.add((chromosome,a))
                splice_site_subset.add((chromosome,str(int(a)-1)))
                splice_site_subset.add((chromosome,str(int(a)+1)))

        all_junctions = False
        find_splice_sites = True
    else:
        all_junctions = True
        find_splice_sites = False
    
    
    if all_junctions and sample_junctions:
        print("All junctions and all samples selected. Exiting...")
        return None
                  
    with open(input_filename) as full_table:
        with open(output_filename,"w") as new_table:
            if all_samples:
                new_table.write(full_table.readline())
                for line in full_table:
                    row = line.strip().split('\t')
                    if find_splice_sites:
                        chromosome,a,b = row[0].replace('-',':').split(':')
                        if (chromosome,a) in splice_site_subset or (chromosome,b) in splice_site_subset:
                            ps_values = [float(x) for x in row[1:]]
                            
                            #if max(ps_values) - min(ps_values) >= args.psRangeFilter:         
                            new_table.write(line)
                            
                    elif row[0] in junction_subset:
                        ps_values = [float(x) for x in row[1:]]
                        
                        #if max(ps_values) - min(ps_values) >= args.psRangeFilter:    
                        new_table.write(line)
                        
            else:
                tab = '\t'
                sample_indices = [0]
                for i,sample in enumerate(full_table.readline().strip().split('\t')):
                    if sample in sample_subset:
                        sample_indices.append(i)
                        sample_order.append(sample)
                new_table.write(f"cluster\t{tab.join(sample_order)}\n")
                for line in full_table:
                    row = line.strip().split('\t')
                    if find_splice_sites:
                        chromosome,a,b = line.split("\t",maxsplit=1)[0].replace('-',':').split(':')
                        if (chromosome,a) in splice_site_subset or (chromosome,b) in splice_site_subset:
                            ps_values = [float(x) for x in row[1:]]
                            
                            #if max(ps_values) - min(ps_vales) >= args.psRangeFilter:    
                            new_table.write(f"{tab.join([row[i] for i in sample_indices])}\n")
                            
                    elif all_junctions or row[0] in junction_subset:
                        ps_values = [float(x) for x in row[1:]]
                        
                        #if max(ps_values) - min(ps_vales) >= args.psRangeFilter:    
                        new_table.write(f"{tab.join([row[i] for i in sample_indices])}\n")
         
            
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    add_parser(parser)
    args = parser.parse_args()
    run_with(args)