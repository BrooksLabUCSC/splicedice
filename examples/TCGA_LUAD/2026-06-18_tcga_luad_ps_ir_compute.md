# 2026-06-18_tcga_luad_ps_ir_compute

## Background

**Author**

Holly Beale

**Repos**

https://github.com/BrooksLabUCSC/splicedice

https://github.com/hbeale/splicedice_analysis

**Output**

Zenodo



**Context**

This notebook documents a splicedice run that generated percent spliced and intron retention values from TCGA lung adenocarcinoma (LUAD) samples.

This approach is designed for limited disk space. Bams are downloaded one at a time (for identifying introns with intron-prospector) or in small batches (for calculating intron coverage with splicedice intron_coverage). 

It uses tools from Holly Beale's splicedice_analysis repo. Those scripts have been copied into the same directory as this file.



# Setup per instance

## copy gdc file

```bash
ls ~/gdc-user-token.2026-05-28T20_33_35.481Z.txt
cp /mnt/git_code/gdc-user-token.2026-05-28T20_33_35.481Z.txt ~
```

make sure gdc-client is in the path

```bash
gdc-client
```





# Setup per run

## define location

```bash
this_commit=970d652
this_full_SHA_hash=970d6525e50f01bfd06d695ffd5ad6c41fffaabf
this_description=tcga_luad
timestamp=2026.06.18_10.26.39
this_base_dir=/mnt/sd/${this_description}_${this_commit}_${timestamp}/
code_base=${this_base_dir}/git_code/splicedice_analysis/2026_06_ps_ir_pipeline
working_files=${this_base_dir}/git_code/splicedice_analysis/2026-06_TCGA_IP_splicedice_PS_compute
mkdir -p ${this_base_dir}/git_code/ ${this_base_dir}/analysis/ ${this_base_dir}/intron_beds ${this_base_dir}/bams
```

## get code

```bash
cd ${this_base_dir}/git_code/
git clone https://github.com/hbeale/splicedice_analysis.git
```

## build docker

a few minutes

```bash
this_dockerfile=${working_files}/Dockerfile_$this_commit
cat $code_base/Dockerfile_splicedice_by_hash | sed "s/replace_with_hash/${this_full_SHA_hash}/" > $this_dockerfile
docker build --build-arg CACHE_BUST=$(date +%s) -t splicedice_analysis:latest -f $this_dockerfile .
bash ~/alert_msg.sh "docker build complete"
```



## make manifests

```bash
this_manifest=${working_files}/manifests/primary_manifest.txt
mkdir -p `dirname $this_manifest`
cat ${code_base}/manifests/primary_manifest.txt | sed "s|/mnt/data/tcga|${this_base_dir}/bams|" | \
sed "s|/mnt/data/intron_prospector_runs/common|${this_base_dir}/intron_beds|" > $this_manifest

cat $this_manifest | grep -v dataset_id | cut -f1,3,4 > ${this_base_dir}/analysis/quant_manifest.txt

```



# Run pipeline

## run intron-prospector

```bash
date
time bash ${code_base}/scripts/run_intron-prospector.sh \
    --manifest $this_manifest \
    --genome /mnt/ref/GRCh38.primary_assembly.genome.fa \
    --disk-constraint yes
date
~/alert_msg.sh "ip_done"

```

std out

```bash
Thu Jun 18 17:48:58 UTC 2026

processing TCGA-86-8074-01A
neither bed or bam file exists
downloading bam file...
using token: /home/ubuntu/gdc-user-token.2026-05-28T20_33_35.481Z.txt
...
processing TCGA-78-8660-01A
neither bed or bam file exists
downloading bam file...
using token: /home/ubuntu/gdc-user-token.2026-05-28T20_33_35.481Z.txt
100% [#######################################################################################################################] Time:  0:01:38  53.7 MiB/s
100% [#######################################################################################################################] Time:  0:00:04   1.0 MiB/s 
ERROR: ('Connection aborted.', ConnectionResetError(104, 'Connection reset by peer'))
WARNING: Unable to download annotations for 4ee7ff21-a0ae-4885-92d2-a088d6f87cf0: 'NoneType' object has no attribute 'raise_for_status'
Successfully downloaded: 1
bam file exists but bed file does not
running intron-prospector...
...
Warning: genomic sequence not found for chrUn_JTFH01001084v1_decoy splice junctions not available
Warning: genomic sequence not found for chrUn_JTFH01001171v1_decoy splice junctions not available
Warning: genomic sequence not found for chrUn_JTFH01001241v1_decoy splice junctions not available
both bed and bam files now exist
deleting bam file

real    2456m17.853s
user    830m22.349s
sys     306m55.285s
Sat Jun 20 10:45:16 UTC 2026
{"status":"OK","nsent":2,"apilimit":"0\/1000"}

```

41 hours

i don't know if the error was temporary or not. probably a good reason to check that all the bed files exist and have the expected chromosomes present. 



### check outputs

```bash
bash ${code_base}/scripts/check_intron_prospector_outputs.sh \
    --manifest $this_manifest
```

std out

```bash
...
TCGA-44-3398-01A: OK
TCGA-64-5774-01A: OK

SUMMARY: 495 datasets | 495 OK | 0 missing BED | 0 missing chromosomes
```



## run quant

```bash
date
time docker run --rm \
-v /mnt/:/mnt \
splicedice_analysis:latest \
splicedice quant -m ${this_base_dir}/analysis/quant_manifest.txt \
-o ${this_base_dir}/analysis/
date
 ~/alert_msg.sh "quant run complete"

```



std out

```bash
ubuntu@hbeale-clin-validation:/mnt$ date
time docker run --rm \
-v /mnt/:/mnt \
splicedice_analysis:latest \
splicedice quant -m ${this_base_dir}/analysis/quant_manifest.txt \
-o ${this_base_dir}/analysis/
date
 ~/alert_msg.sh "quant run complete"
Sat Jun 20 21:20:42 UTC 2026
Parsing manifest...
        Done [0:00:0.00]
Getting all junctions from 495 files...
        Done [0:03:28.53]
Finding clusters from 800108 junctions...
        Done [0:00:35.54]
Writing cluster file...
        Done [0:01:28.09]
Writing junction bed file...
        Done [0:00:5.55]
Gathering junction counts...
        Done [0:05:29.80]
Writing inclusion counts...
        Done [0:06:29.63]
Calculating PS values...
        Done [0:06:21.90]
Writing PS values...
        Done [0:07:28.58]
All done [0:31:27.63]

real    31m33.487s
user    0m0.183s
sys     0m0.063s
Sat Jun 20 21:52:15 UTC 2026
{"status":"OK","nsent":2,"apilimit":"2\/1000"}
ubuntu@hbeale-clin-validation:/mnt$ 


```





## run intron_coverage

```bash
batch_size=16

date
time bash ${code_base}/scripts/run_intron_coverage_pipeline.sh \
    --manifest $this_manifest \
    --analysis-base ${this_base_dir}/analysis/ \
    --disk-constraint yes \
    --batch-size $batch_size
date
~/alert_msg.sh "intron_coverage run complete"
 
 
```



std out

```bash
Sat Jun 20 21:57:05 UTC 2026
=== 495 samples still need intron_coverage ===

========================================================
BATCH 1: 16 samples
========================================================

--- downloading BAMs for batch 1 ---
using token: /home/ubuntu/gdc-user-token.2026-05-28T20_33_35.481Z.txt
TCGA-86-8074-01A: queuing 567c5d5f-2b27-4070-86c3-3905d06ed02b for download
TCGA-62-8402-01A: queuing cae0680e-f7bf-4742-aeca-8fac6d4f4934 for download
TCGA-86-8358-01A: queuing e5976aee-2a56-457c-80e5-00824254f6f8 for download
TCGA-86-8056-01A: queuing 6dee9448-b65a-498e-9490-7c282fb3b07d for download
TCGA-78-7158-01A: queuing 4d6609e2-6ad0-43a1-9bda-fbc4710e1da0 for download
TCGA-49-4507-01A: queuing a6f82885-da83-49d5-ad43-966b9dff4ea5 for download
TCGA-49-AARO-01A: queuing 1f8160a9-4e65-4d8f-a83d-39b6e03b38f3 for download
TCGA-91-8499-01A: queuing 073e7f71-5583-48fc-a037-4f799ec2d811 for download
TCGA-55-6983-01A: queuing af2f19ca-dc08-43b6-ae40-811cb952887a for download
TCGA-62-A46Y-01A: queuing 0fd61331-8363-402b-8d87-88bb39f467d0 for download
TCGA-53-7624-01A: queuing 888e7914-72dc-46ab-b34f-84e30cf4ede8 for download
TCGA-L9-A7SV-01A: queuing 0bd949c3-479b-4fec-a3e3-83a0e0f437a0 for download
TCGA-62-8395-01A: queuing 7bd81e9a-6853-4c8e-ba87-849957a8f015 for download
TCGA-73-4662-01A: queuing 1c0c1b78-5622-4da2-a6ce-d00c37c51915 for download
TCGA-97-7547-01A: queuing 93c69bd1-3426-4475-a59a-b5f6c25222b1 for download
TCGA-35-5375-01A: queuing a5b878f4-5fe7-4f53-a61f-d653d435b24c for download
Downloading 16 BAM(s)...
Sat Jun 20 21:57:09 UTC 2026
  5% [######              ] ETA:   0:05:45  28.9 MiB/s 
...
...
...
100% [#######################################################################################################################] Time:  0:03:27  28.1 MiB/s 
100% [#######################################################################################################################] Time:  0:04:19  26.7 MiB/s 
100% [#######################################################################################################################] Time:  0:03:29  28.8 MiB/s 
Successfully downloaded: 16
Downloads complete.
Mon Jun 22 15:56:19 UTC 2026

--- running intron_coverage for batch 12 ---
Running intron_coverage on 16 samples with 8 threads...
Mon Jun 22 15:56:19 UTC 2026
getting paths for bam files
creating junction percentiles
[W::hts_idx_load3] The index file is older than the data file: /mnt/sd/tcga_luad_970d652_2026.06.18_10.26.39//bams/8f82532e-91ca-469b-a1a5-9bca66844679/04c4affd-ca5d-4ea4-95d1-a6502e80d1af.rna_seq.genomic.gdc_realn.bam.bai
[W::hts_idx_load3] The index file is older than the data file: /mnt/sd/tcga_luad_970d652_2026.06.18_10.26.39//bams/abac4b21-6420-4324-95be-869ffaadc4d9/8c55ee80-c8c3-4196-9b67-81732fffe9f4.rna_seq.genomic.gdc_realn.bam.bai
...
...
...
TCGA-75-6206-01A: output confirmed, deleting BAM
TCGA-55-7907-01A: output confirmed, deleting BAM
TCGA-55-6972-01A: output confirmed, deleting BAM
TCGA-55-7283-01A: output confirmed, deleting BAM
TCGA-97-A4M1-01A: output confirmed, deleting BAM
TCGA-78-7160-01A: output confirmed, deleting BAM
TCGA-49-4494-01A: output confirmed, deleting BAM
TCGA-44-A47G-01A: output confirmed, deleting BAM
TCGA-55-7570-01A: output confirmed, deleting BAM
TCGA-78-7155-01A: output confirmed, deleting BAM
TCGA-44-3398-01A: output confirmed, deleting BAM
TCGA-64-5774-01A: output confirmed, deleting BAM

batch 31 complete

========================================================
All batches complete.
========================================================

real    6368m17.960s
user    789m47.782s
sys     262m20.919s
Thu Jun 25 08:05:23 UTC 2026

```

took 4.4 days



## Create intron table

After several issues with out of memory errors, this code worked.

build docker from branch (note; this branch was incoporated into the main splicedice branch on 7/2/2026, commit id 1f6ffe4)

```bash
cd /mnt/git_code/
git clone https://github.com/BrooksLabUCSC/splicedice.git splicedice
cd splicedice
git checkout ir-table-low-memory # commit 62773ae
git pull

docker build --build-arg CACHE_BUST=$(date +%s) \
    -t splicedice_analysis:ir_table_low_memory \
    -f /mnt/git_code/splicedice/validation/processes/Validation_ir-table-low-memory_branch/Dockerfile_ir_table_low_memory.txt .
~/alert_msg.sh "docker build complete"
```



#### start ir_table

```bash
genes=/mnt/ref/gencode.v47.primary_assembly.annotation.gtf
date
time docker run --rm \
-e PYTHONUNBUFFERED=1 \
-v /mnt/:/mnt \
splicedice_analysis:ir_table_low_memory \
splicedice ir_table \
--annotation $genes \
-i ${this_base_dir}/analysis/_inclusionCounts.tsv \
-c ${this_base_dir}/analysis/_allClusters.tsv \
-d ${this_base_dir}/analysis/coverage_output \
-n 8 \
-o ${this_base_dir}/analysis/
date
~/alert_msg.sh intron_table_creation_complete 


```

std out

```bash
Fri Jun 26 00:57:52 UTC 2026
Loading annotation...
Annotation loaded: 528735 annotated junctions. 55.2s
Gathering inclusion counts and clusters...
Loaded 495 samples and 800108 clusters. 203.5s
Collecting junctions across all samples...
getJunctions complete: 327831 junctions. 331.3s
RSD filtering complete: 251731 junctions retained. 1623.2s
Junction collection and RSD filtering complete: 251731 junctions retained. 1826.7s
Writing IR table...
IR calculated for 50/495 samples
IR calculated for 100/495 samples
IR calculated for 150/495 samples
IR calculated for 200/495 samples
IR calculated for 250/495 samples
IR calculated for 300/495 samples
IR calculated for 350/495 samples
IR calculated for 400/495 samples
IR calculated for 450/495 samples
IR calculated for 495/495 samples
IR table written. 11080.5s
Done. Total runtime: 11080.5s

real    184m49.284s
user    0m0.858s
sys     0m0.160s
Fri Jun 26 04:02:41 UTC 2026


```



## Copy files

```bash
cd ${this_base_dir}/analysis/
pigz -k _allPS.tsv
mv _allPS.tsv.gz TCGA_LUAD_allPS.tsv.gz

pigz -k _intron_retention.tsv
mv _intron_retention.tsv.gz TCGA_LUAD_intron_retention.tsv.gz

# from mustard
scp ubuntu@10.50.100.128:/${this_base_dir}/analysis/TCGA_LUAD_allPS.tsv.gz .
scp ubuntu@10.50.100.128:/${this_base_dir}/analysis/TCGA_LUAD_intron_retention.tsv.gz .

```

