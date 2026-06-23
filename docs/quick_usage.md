---
layout: default
---

## Quick usage guide

### Step 1.

#### Option 1: Convert BAM to a BED with junctions using leafcutter

``` bash
~/leafcutter/scripts/bam2junc.sh file.bam file.junc
```

#### Option 2: Convert STAR junction file to junction BED

``` bash

```
### Step 2.

#### Create junction file manifest
Create a tsv with the four following columns:
`sample_id, junction_file_location, group1, group2`



### Step 3.

#### Build MESA events

```bash
python3 ~/mesa/constructMESA.py -m manifest.txt -f genome.fa -o output_prefix
```
output - `output_prefix_juncs.bed, output_prefix_all_clusters2.tsv`

### Step 4.

#### Quantify MESA events

```bash
python3 ~/mesa/quantMesa.py -m manifest.txt -j step1_junctions.bed -c step1_clusters.bed -o outPrefix -o output_prefix
```



### [Go Home](./) or [Go to Installation](./installation.html)
