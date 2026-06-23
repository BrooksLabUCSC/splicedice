# SUGP1_walkthrough_2026.06.23_11.24.14

Assumes you've downloaded and aligned RNA-Seq data from SRP286876, representing experiments comparing HEK293T cells transfected with control siRNA to those transfected with siSUGP1.  There are three  biological replicates of each condition. 

SRR ids: SRR12801019 SRR12801020 SRR12801023 SRR12801024 SRR12801027 SRR12801028



# Setup

## Server requirements

Docker 

## define location

```bash
this_commit=
this_description=SUGP1_walkthrough
this_datestamp=2026.06.23_11.14.13
this_branch=SUGP1_walkthrough
this_dockerfile=Dockerfile_for_${this_description}.txt
```

```bash
this_base_dir=/mnt/sd/${this_description}_${this_commit}_${this_datestamp}/
code_base=${this_base_dir}/git_code/splicedice/examples/SUGP1_walkthrough/
mkdir -p ${this_base_dir}/git_code/ ${this_base_dir}/analysis/ ${this_base_dir}/intron_beds

```

git clone --depth 2  --branch Remove-untested-code  \
      https://github.com/BrooksLabUCSC/splicedice.git splicedice_Remove-untested-code

## get code

```bash
cd ${this_base_dir}/git_code/
git clone --depth 2  --branch $this_branch \
https://github.com/BrooksLabUCSC/splicedice.git 
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

cat $primary_manifest | \
sed "s/replace_with_bam_base/$this_base_dir/" | \
cut -f1,2,4 > $bam_manifest

cat $primary_manifest | \
sed "s/replace_with_bed_base/$this_base_dir/"  | \
cut -f1,3,4 > $bed_manifest
```



## view manifest contents

bed_manifest

```bash
cat $bed_manifest
```

```bash
SRR12801019     /mnt/data/intron_prospector_runs/2026-05-05_20-28-22/SRR12801019.bed    control
SRR12801020     /mnt/data/intron_prospector_runs/2026-05-05_20-28-22/SRR12801020.bed    SUGP1_kd
SRR12801023     /mnt/data/intron_prospector_runs/2026-05-05_20-28-22/SRR12801023.bed    control
SRR12801024     /mnt/data/intron_prospector_runs/2026-05-05_20-28-22/SRR12801024.bed    SUGP1_kd
SRR12801027     /mnt/data/intron_prospector_runs/2026-05-05_20-28-22/SRR12801027.bed    control
SRR12801028     /mnt/data/intron_prospector_runs/2026-05-05_20-28-22/SRR12801028.bed    SUGP1_kd

```



bam_manifest

```bash
cat $bam_manifest
```

```bash
SRR12801019     /mnt/output/star_2.7.11b_2026.04.16/SRR12801019/SRR12801019.bam control control
SRR12801020     /mnt/output/star_2.7.11b_2026.04.16/SRR12801020/SRR12801020.bam SUGP1_kd        SUGP1_kd
SRR12801023     /mnt/output/star_2.7.11b_2026.04.16/SRR12801023/SRR12801023.bam control control
SRR12801024     /mnt/output/star_2.7.11b_2026.04.16/SRR12801024/SRR12801024.bam SUGP1_kd        SUGP1_kd
SRR12801027     /mnt/output/star_2.7.11b_2026.04.16/SRR12801027/SRR12801027.bam control control
SRR12801028     /mnt/output/star_2.7.11b_2026.04.16/SRR12801028/SRR12801028.bam SUGP1_kd        SUGP1_kd

```





# Identify introns with intron-prospector

```bash

genome=/mnt/ref/GRCh38.primary_assembly.genome.fa

cat $primary_manifest | while read id bam bed phenotype; 
do

docker run --rm \
-v /mnt/:/mnt \
splicedice_analysis:latest \
intron-prospector \
--genome-fasta=$genome \
--intron-bed6=$bed \
$bam

done
bash ~/alert_msg.sh intron_prospector_complete 

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

