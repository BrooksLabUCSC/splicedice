import numpy as np
import os, sys

good = set()
bad = set()
with open(sys.argv[1]) as lines:
    for line in lines:
        c = line.rstrip().split()
        g1, g2 = c[1], c[2]
        
        if "Normal" in line:
            good.add(c[0])
        elif "Tumor" in line:
            bad.add(c[0])

o1 = open(sys.argv[1].split("_")[0] + "_" + "psi.out",'w')
o2 = open(sys.argv[1].split("_")[0] + "_" + "iqr.out",'w')
o3 = open(sys.argv[1].split("_")[0] + "_" + "data.out",'w')
o4 = open(sys.argv[1].split("_")[0] + "_" + "zscore.out",'w')

with open(sys.argv[2]) as lines:
    header = next(lines).rstrip().split()
    norm = [pos for pos,x in enumerate(header) if x in good]
    tumo = [pos for pos,x in enumerate(header) if x in bad]
    
    for line in lines:
        cols = line.rstrip().split()
        vals = np.asarray([x.split(":") for x in cols[:-2]],dtype=np.float32)
     
        
        sums = vals.sum(axis=1)
        if np.nanmax(sums)<35:
            continue

        sums[sums < 35] = np.nan



        psi = vals[:,0] / sums

        psiT = psi[tumo]
        psiN = psi[norm]

        tCounts = np.count_nonzero(np.isnan(psiT))
        nCounts = np.count_nonzero(np.isnan(psiN))

        if tCounts/len(tumo) > 0.25 or nCounts/len(norm) > 0.25 or len(psiN)<10:
            continue

        print("\t".join(("%.2f" % x for x in psiT)),"\t".join(("%.2f" % x for x in psiN)), cols[-2], cols[-1], sep="\t", file=o1)

        iqr = np.nanpercentile(psiN, 75) - np.nanpercentile(psiN, 25)
        
        if iqr<0.01:
            continue
        
        nmedian = np.nanmedian(psiN)
        tiqr = (psiT - nmedian) / iqr
        niqr = (psiN - nmedian) / iqr

        outliers = np.absolute(tiqr) > 3
        numOut = outliers.sum()

        avgm   = np.nanmean(psiT[outliers])

        print(cols[-2],"%.2f" % iqr, "%.2f" %  nmedian, numOut, "%.2f" % avgm,sep="\t", file=o3)
        print("\t".join(("%.2f" % x for x in tiqr)), "\t".join(("%.2f" % x for x in niqr)), cols[-2], cols[-1], sep="\t", file=o2)

        std = np.nanstd(psiN)
        if std<0.01: continue
        tz = (psiT - nmedian) / std
        nz = (psiN - nmedian) / std

        print("\t".join(("%.2f" % x for x in tz)), "\t".join(("%.2f" % x for x in nz)), cols[-2], cols[-1], sep="\t", file=o4)



o1.close()
o2.close()
o3.close()