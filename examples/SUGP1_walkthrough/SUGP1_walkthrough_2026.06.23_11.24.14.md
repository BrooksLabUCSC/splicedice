# SUGP1_walkthrough_2026.06.23_11.24.14

Assumes you've downloaded and aligned RNA-Seq data from SRP286876, representing experiments comparing HEK293T cells transfected with control siRNA to those transfected with siSUGP1.  There are three  biological replicates of each condition. 

SRR ids: SRR12801019 SRR12801020 SRR12801023 SRR12801024 SRR12801027 SRR12801028



# Setup

## Server requirements

Docker 

## define location

```bash
this_commit=bc373ae
this_description=SUGP1_walkthrough
this_datestamp=2026.06.23_11.14.13
this_dockerfile=Dockerfile_for_${this_description}.txt
```

```bash
this_base_dir=/mnt/sd/${this_description}_${this_commit}_${this_datestamp}/
code_base=${this_base_dir}/git_code/splicedice/examples/SUGP1_walkthrough/
mkdir -p ${this_base_dir}/git_code/ ${this_base_dir}/analysis/ ${this_base_dir}/intron_beds

```


## get code

```bash
cd ${this_base_dir}/git_code/
git clone --depth 5 \
https://github.com/BrooksLabUCSC/splicedice.git 
cd splicedice
SHA1_splicedice="bc373ae02c20ce61412b2f9ba3229844b684e2dd"
git reset --hard $SHA1_splicedice

```

## build docker

```bash

docker build --build-arg CACHE_BUST=$(date +%s) -t splicedice_analysis:latest -f ${code_base}/${this_dockerfile} .
bash ~/alert_msg.sh "docker build complete"

```

aside: ~/alert_msg.sh is a personal convenience script that notifies me when the command is run. it's included here so I don't forget to include it when I use this as a template, but it can be ignored or replaced with your own notification script



## create manifests for this run

```bash
primary_manifest=${code_base}/primary_manifest.tsv
bed_manifest=${primary_manifest/primary/bed}
bam_manifest=${primary_manifest/primary/bam}

cat ${code_base}/generic_primary_manifest.tsv | \
sed "s|replace_with_bam_base|$this_base_dir|" | \
sed "s|replace_with_bed_base|$this_base_dir|" \
> $primary_manifest

cat $primary_manifest | cut -f1,2,4 > $bam_manifest
cat $primary_manifest | cut -f1,3,4 > $bed_manifest
```

(note; if your bams aren't already in the location in the manfest, link to them, e.g. `for i in `ls /mnt/output/star_2.7.11b_2026.04.16/`; do echo $i; ln -s /mnt/output/star_2.7.11b_2026.04.16/${i}/${i}.bam /mnt/sd/SUGP1_walkthrough_bc373ae_2026.06.23_11.14.13//bams/; done`)



## view manifest contents

bed_manifest

```bash
cat $bed_manifest
```

```bash
SRR12801019     /mnt/sd/SUGP1_walkthrough_bc373ae_2026.06.23_11.14.13//intron_beds/SRR12801019.bed      control
SRR12801020     /mnt/sd/SUGP1_walkthrough_bc373ae_2026.06.23_11.14.13//intron_beds/SRR12801020.bed      SUGP1_kd
SRR12801023     /mnt/sd/SUGP1_walkthrough_bc373ae_2026.06.23_11.14.13//intron_beds/SRR12801023.bed      control
SRR12801024     /mnt/sd/SUGP1_walkthrough_bc373ae_2026.06.23_11.14.13//intron_beds/SRR12801024.bed      SUGP1_kd
SRR12801027     /mnt/sd/SUGP1_walkthrough_bc373ae_2026.06.23_11.14.13//intron_beds/SRR12801027.bed      control
SRR12801028     /mnt/sd/SUGP1_walkthrough_bc373ae_2026.06.23_11.14.13//intron_beds/SRR12801028.bed      SUGP1_kd
```



bam_manifest

```bash
cat $bam_manifest
```

```bash

SRR12801019     /mnt/sd/SUGP1_walkthrough_bc373ae_2026.06.23_11.14.13//bams/SRR12801019.bam     control
SRR12801020     /mnt/sd/SUGP1_walkthrough_bc373ae_2026.06.23_11.14.13//bams/SRR12801020.bam     SUGP1_kd
SRR12801023     /mnt/sd/SUGP1_walkthrough_bc373ae_2026.06.23_11.14.13//bams/SRR12801023.bam     control
SRR12801024     /mnt/sd/SUGP1_walkthrough_bc373ae_2026.06.23_11.14.13//bams/SRR12801024.bam     SUGP1_kd
SRR12801027     /mnt/sd/SUGP1_walkthrough_bc373ae_2026.06.23_11.14.13//bams/SRR12801027.bam     control
SRR12801028     /mnt/sd/SUGP1_walkthrough_bc373ae_2026.06.23_11.14.13//bams/SRR12801028.bam     SUGP1_kd
```





# Identify introns with intron-prospector

```bash

genome=/mnt/ref/GRCh38.primary_assembly.genome.fa

date
cat $primary_manifest | while read id bam bed phenotype; 
do
echo $id

docker run --rm \
-v /mnt/:/mnt \
splicedice_analysis:latest \
intron-prospector \
--genome-fasta=$genome \
--intron-bed6=$bed \
$bam

done
date
bash ~/alert_msg.sh intron_prospector_complete 

```

std out

```bash
Tue Jun 23 20:22:03 UTC 2026
SRR12801019
SRR12801020
SRR12801023
SRR12801024
SRR12801027
SRR12801028
Tue Jun 23 20:28:59 UTC 2026
```





# Identify and quantify splice junctions

```bash
date
time docker run --rm \
-v /mnt/:/mnt \
splicedice_analysis:latest \
splicedice quant -m ${bed_manifest} \
-o ${this_base_dir}/analysis/
date
 ~/alert_msg.sh "quant run complete"

```

std out

```bash
Tue Jun 23 20:42:10 UTC 2026
Parsing manifest...
        Done [0:00:0.00]
Getting all junctions from 6 files...
        Done [0:00:3.92]
Finding clusters from 333069 junctions...
        Done [0:00:3.38]
Writing cluster file...
        Done [0:00:3.39]
Writing junction bed file...
        Done [0:00:2.28]
Gathering junction counts...
        Done [0:00:5.27]
Writing inclusion counts...
        Done [0:00:4.56]
Calculating PS values...
        Done [0:00:14.32]
Writing PS values...
        Done [0:00:4.72]
All done [0:00:41.85]

real    0m45.759s
user    0m0.044s
sys     0m0.056s
Tue Jun 23 20:42:55 UTC 2026

```


# Calculate intron coverage

about 25 minutes

```bash
n_threads=8
mkdir ${this_base_dir}/coverage_output/
date
time docker run --rm \
    -v /mnt/:/mnt \
    splicedice_analysis:latest \
    splicedice intron_coverage \
    -b $bam_manifest  \
    -j ${this_base_dir}/analysis/_junctions.bed \
    -n "${n_threads}" \
    -o ${this_base_dir}/coverage_output/
date
~/alert_msg.sh "intron_coverage run complete"


```

```bash
Tue Jun 23 20:49:56 UTC 2026
getting paths for bam files
creating junction percentiles
SRR12801019 starting 5.3634960651397705
SRR12801019 collected 628.1782286167145
SRR12801019 counted 1267.9633173942566
SRR12801019 done 1282.1344695091248
SRR12801024 starting 7.4829628467559814
SRR12801024 collected 795.3424100875854
SRR12801024 counted 1460.5932395458221
SRR12801024 done 1474.2689065933228
SRR12801023 starting 6.865115642547607
SRR12801023 collected 690.1928341388702
SRR12801023 counted 1351.9379482269287
SRR12801023 done 1365.303383588791
SRR12801028 starting 9.051880359649658
SRR12801028 collected 846.2235708236694
SRR12801028 counted 1570.6379480361938
SRR12801028 done 1584.6428875923157
SRR12801027 starting 8.246123313903809
SRR12801027 collected 808.1194612979889
SRR12801027 counted 1570.9109573364258
SRR12801027 done 1584.7905428409576
Your runtime was 1608.330510854721 seconds.

real    26m52.394s
user    0m0.133s
sys     0m0.055s
Tue Jun 23 21:16:48 UTC 2026

```

# Create intron table

about 2 minutes

```bash
genes=/mnt/ref/gencode.v47.primary_assembly.annotation.gtf
date
time docker run --rm \
-v /mnt/:/mnt \
splicedice_analysis:latest \
splicedice ir_table \
--annotation $genes \
-i ${this_base_dir}/analysis/_inclusionCounts.tsv \
-c ${this_base_dir}/analysis/_allClusters.tsv \
-d ${this_base_dir}/coverage_output \
-n 8 \
-o ${this_base_dir}/analysis/
date
~/alert_msg.sh intron_table_creation_complete 

```

```bash
Tue Jun 23 21:37:33 UTC 2026
Starting ir_table with 6 samples
Loading annotation...
Annotation loaded: 528735 annotated junctions. 60.2s
Gathering inclusion counts and clusters...
Loaded 6 samples and 333069 clusters. 69.2s
Collecting junctions across all samples...
getJunctions complete: 210470 junctions. 5.3s
RSD filtering complete: 96399 junctions retained. 23.0s
Junction collection and RSD filtering complete: 96399 junctions retained. 92.3s
Writing IR table...
IR calculated for 6/6 samples
IR table written. 149.7s
Done. Total runtime: 149.7s

real    2m35.469s
user    0m0.057s
sys     0m0.058s
Tue Jun 23 21:40:08 UTC 2026
```

# Compare to previous output

# 

```bash
prev_dir=/mnt/splicedice_ir_example_archives/2026.05.28_21.25.54/analysis/
current_dir=${this_base_dir}/analysis
for i in `ls $current_dir`; do echo $i; diff $current_dir/$i $prev_dir/${i}; done
```

```bash
ubuntu@hbeale-mesa:/mnt/sd/SUGP1_walkthrough_bc373ae_2026.06.23_11.14.13/git_code/splicedice$ prev_dir=/mnt/splicedice_ir_example_archives/2026.05.28_21.25.54/analysis/
current_dir=${this_base_dir}/analysis
for i in `ls $current_dir`; do echo $i; diff $current_dir/$i $prev_dir/${i}; done
_allClusters.tsv
_allPS.tsv
_inclusionCounts.tsv
_intron_retention.tsv
_junctions.bed
ubuntu@hbeale-mesa:/mnt/sd/SUGP1_walkthrough_bc373ae_2026.06.23_11.14.13/git_code/splicedice$
```

make sure my code works:

```bash
echo 1 > $prev_dir/not_matching
echo 2 > $current_dir/not_matching
for i in `ls $current_dir`; do echo $i; diff $current_dir/$i $prev_dir/${i}; done
```

```bash
ubuntu@hbeale-mesa:/mnt/sd/SUGP1_walkthrough_bc373ae_2026.06.23_11.14.13/git_code/splicedice$ for i in `ls $current_dir`; do echo $i; diff $current_dir/$i $prev_dir/${i}; done
_allClusters.tsv
_allPS.tsv
_inclusionCounts.tsv
_intron_retention.tsv
_junctions.bed
not_matching
1c1
< 2
---
> 1
ubuntu@hbeale-mesa:/mnt/sd/SUGP1_walkthrough_bc373ae_2026.06.23_11.14.13/git_code/splicedice$ 

```

