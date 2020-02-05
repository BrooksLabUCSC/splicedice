import sys, os

star_tab = sys.argv[1]
min_overhang = 5
min_unique = 5
add_multimap = True

with open("%s.bed" % sys.argv[1], 'w') as fout:
    with open(sys.argv[1]) as fin:
        for i in fin:
            cols = i.rstrip().split()
            chrom,c1,c2,strand = cols[:4]
            c1,c2,strand = int(c1),int(c2),int(strand)
            unique = int(cols[6])
            multi = int(cols[7])
            overhang = int(cols[-1])

            strand = "+" if strand == 1 else "-" 
            if overhang < min_overhang or unique < min_unique or strand == 0:
                continue
            if add_multimap:
                print(chrom,c1-1,c2,".",unique+multi,strand,sep="\t",file=fout)
            else:
                print(chrom,c1-1,c2,".",unique,strand,sep="\t",file=fout)
