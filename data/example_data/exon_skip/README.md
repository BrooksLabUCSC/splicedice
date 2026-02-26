# Exon Skipping BED Files

These files model a single exon skip event on `chr1` (`+` strand) using IP BED6 format:

chrom | start | end | name (sj#) | score (read count) | strand

## Coordinate Setup (same for all files)

Upstream exon ends at: 1000  
Variable exon: 1100–1200  
Downstream exon starts at: 1300  

Splice junctions:
- sj1: 1000–1100  (upstream → variable)
- sj2: 1200–1300  (variable → downstream)
- sj3: 1000–1300  (skipping)

## Files

- `ps_0.bed`  
  Only sj3 has reads → PS = 0%

- `ps_50.bed`  
  sj1/sj2 equal to sj3 → PS = 50%

- `ps_100.bed`  
  Only sj1/sj2 have reads → PS = 100%