# Validation: ir-table-low-memory branch

Validates that the memory optimizations in `ir-table-low-memory` produce identical `_intron_retention.tsv` output to the known-good reference.

Changes validated:

- `getInclusionCounts`: sparse dict reader replacing pandas DataFrame pivot
- `getJunctions`: `imap_unordered` with incremental union replacing `pool.map`
- `getFilteredJunctions`: `imap_unordered` with incremental union replacing `pool.map`



# Server

ubuntu@hbeale-mesa



# Attempt 1

commit ebbfbd620d01632a8d62698e29781bd3c2d72de0

## Setup

### Define locations

```bash
this_base_dir=/mnt/sd/validation_ir_table_low_memory_2026.06.25_11.38.03/
known_good_dir=/mnt/sd/known_good/splicedice/validation/data/SUGP1_SRP286876
validation_source=/mnt/sd/validation_2b2c401_2026.06.23_15.07.59
genes=/mnt/ref/gencode.v47.primary_assembly.annotation.gtf

mkdir -p ${this_base_dir}/analysis/
```

### Get known-good reference data

(not necessary; it's already present)

```bash
mkdir -p /mnt/sd/known_good
cd /mnt/sd/known_good
git clone https://github.com/BrooksLabUCSC/splicedice.git
```

### Build Docker image from ir-table-low-memory branch

```bash
cd /mnt/gitCode/
git clone https://github.com/BrooksLabUCSC/splicedice.git splicedice
cd splicedice
git checkout ir-table-low-memory

docker build --build-arg CACHE_BUST=$(date +%s) \
    -t splicedice_analysis:ir_table_low_memory \
    -f /mnt/gitCode/splicedice/validation/processes/Validation_ir-table-low-memory_branch/Dockerfile_ir_table_low_memory.txt .
~/alert_msg.sh "docker build complete"

```

### prepare input

```bash
zcat ${known_good_dir}/_inclusionCounts.tsv.gz > ${this_base_dir}/analysis/_inclusionCounts.tsv
zcat ${known_good_dir}/_allClusters.tsv.gz > ${this_base_dir}/analysis/_allClusters.tsv
```



## Run ir_table

We reuse `_inclusionCounts.tsv`, `_allClusters.tsv`, and coverage files from the known-good validation run (same 6-sample SUGP1 dataset), so `quant` and `intron_coverage` do not need to be re-run.

```bash
date
time docker run --rm \
    -v /mnt/:/mnt \
    splicedice_analysis:ir_table_low_memory \
    splicedice ir_table \
    --annotation $genes \
    -i ${this_base_dir}/analysis/_inclusionCounts.tsv \
    -c ${this_base_dir}/analysis/_allClusters.tsv \
    -d ${validation_source}/coverage_output \
    -n 8 \
    -o ${this_base_dir}/analysis/
date
~/alert_msg.sh "ir_table validation run complete"

```

Expected output (from known-good run):

```
Starting ir_table with 6 samples
Loading annotation...
Annotation loaded: 528735 annotated junctions. ~60s
Gathering inclusion counts and clusters...
Loaded 6 samples and 333069 clusters. ~70s
Collecting junctions across all samples...
getJunctions complete: 210470 junctions. ~5s
RSD filtering complete: 96399 junctions retained. ~25s
Junction collection and RSD filtering complete: 96399 junctions retained. ~90s
Writing IR table...
IR calculated for 6/6 samples
IR table written. ~150s
Done. Total runtime: ~150s
```



Actual output

```bash
Thu Jun 25 20:19:00 UTC 2026
Starting ir_table with 6 samples
Loading annotation...
Annotation loaded: 528735 annotated junctions. 55.7s
Gathering inclusion counts and clusters...
Loaded 6 samples and 333069 clusters. 59.3s
Collecting junctions across all samples...
getJunctions complete: 210470 junctions. 5.1s
RSD filtering complete: 96399 junctions retained. 24.0s
Junction collection and RSD filtering complete: 96399 junctions retained. 83.3s
Writing IR table...
IR calculated for 6/6 samples
IR table written. 144.5s
Done. Total runtime: 144.5s
real    2m29.106s
user    0m0.049s
sys     0m0.072s
Thu Jun 25 20:21:29 UTC 2026
```



## Compare output to known good

```bash
diff ${this_base_dir}/analysis/_intron_retention.tsv \
    <(zcat ${known_good_dir}/_intron_retention.tsv.gz)
```

Expected: no output (files identical).



failed

```bash
ubuntu@hbeale-mesa:/mnt/gitCode/splicedice$ zcat ${known_good_dir}/_intron_retention.tsv.gz | head
Junction        SRR12801019     SRR12801023     SRR12801024     SRR12801028     SRR12801027     SRR12801020
GL000008.2:163999-164602:-      0.200   nan     nan     nan     0.000   0.000
GL000009.2:55162-78082:-        0.100   0.000   0.000   0.000   nan     0.083
GL000194.1:112851-114985:-      0.000   0.000   0.000   0.000   0.012   0.005
GL000194.1:11336-20237:-        1.000   0.000   nan     0.000   0.000   nan
GL000194.1:11336-28266:-        0.333   0.111   nan     0.000   0.333   nan
GL000194.1:20374-28266:-        0.333   0.250   nan     0.000   1.000   nan
GL000194.1:54833-55445:-        0.406   0.277   0.348   0.282   0.378   0.667
GL000194.1:55677-112791:-       0.014   0.019   0.011   0.039   0.038   0.008
GL000195.1:138141-139989:+      0.024   0.029   0.014   0.006   0.029   0.009
ubuntu@hbeale-mesa:/mnt/gitCode/splicedice$ head ${this_base_dir}/analysis/_intron_retention.tsv
Junction        SRR12801019     SRR12801023     SRR12801024     SRR12801028     SRR12801027     SRR12801020
GL000008.2:163999-164602:-      nan     nan     nan     nan     0.000   nan
GL000009.2:55162-78082:-        nan     nan     nan     nan     nan     0.083
GL000194.1:112851-114985:-      0.000   0.000   0.000   0.000   0.012   0.005
GL000194.1:11336-20237:-        nan     0.000   nan     nan     0.000   nan
GL000194.1:11336-28266:-        nan     0.111   nan     0.000   nan     nan
GL000194.1:20374-28266:-        0.333   nan     nan     nan     nan     nan
GL000194.1:54833-55445:-        0.406   0.277   0.348   0.282   0.378   0.667
GL000194.1:55677-112791:-       0.014   0.019   0.011   0.039   0.038   0.008
GL000195.1:138141-139989:+      0.024   nan     nan     nan     nan     nan
```



# Attempt 2

commit 15f0ebd1e724466930a66d04ca704d2f3d4a321b

## Setup

### Define locations

```bash
this_base_dir=/mnt/sd/validation_ir_table_low_memory_2026.06.25_13.31.46/
known_good_dir=/mnt/sd/known_good/splicedice/validation/data/SUGP1_SRP286876
validation_source=/mnt/sd/validation_2b2c401_2026.06.23_15.07.59
genes=/mnt/ref/gencode.v47.primary_assembly.annotation.gtf

mkdir -p ${this_base_dir}/analysis/
```

### Get known-good reference data

(not necessary; it's already present)

```bash
mkdir -p /mnt/sd/known_good
cd /mnt/sd/known_good
git clone https://github.com/BrooksLabUCSC/splicedice.git
```

### Build Docker image from ir-table-low-memory branch

```bash
cd /mnt/gitCode/
git clone https://github.com/BrooksLabUCSC/splicedice.git splicedice
cd splicedice
git checkout ir-table-low-memory
git pull

docker build --build-arg CACHE_BUST=$(date +%s) \
    -t splicedice_analysis:ir_table_low_memory \
    -f /mnt/gitCode/splicedice/validation/processes/Validation_ir-table-low-memory_branch/Dockerfile_ir_table_low_memory.txt .
~/alert_msg.sh "docker build complete"

```

### prepare input

```bash
zcat ${known_good_dir}/_inclusionCounts.tsv.gz > ${this_base_dir}/analysis/_inclusionCounts.tsv
zcat ${known_good_dir}/_allClusters.tsv.gz > ${this_base_dir}/analysis/_allClusters.tsv
```



## Run ir_table

We reuse `_inclusionCounts.tsv`, `_allClusters.tsv`, and coverage files from the known-good validation run (same 6-sample SUGP1 dataset), so `quant` and `intron_coverage` do not need to be re-run.

```bash
date
time docker run --rm \
    -v /mnt/:/mnt \
    splicedice_analysis:ir_table_low_memory \
    splicedice ir_table \
    --annotation $genes \
    -i ${this_base_dir}/analysis/_inclusionCounts.tsv \
    -c ${this_base_dir}/analysis/_allClusters.tsv \
    -d ${validation_source}/coverage_output \
    -n 8 \
    -o ${this_base_dir}/analysis/
date
~/alert_msg.sh "ir_table validation run complete"

```

Expected output (from known-good run):

```
Starting ir_table with 6 samples
Loading annotation...
Annotation loaded: 528735 annotated junctions. ~60s
Gathering inclusion counts and clusters...
Loaded 6 samples and 333069 clusters. ~70s
Collecting junctions across all samples...
getJunctions complete: 210470 junctions. ~5s
RSD filtering complete: 96399 junctions retained. ~25s
Junction collection and RSD filtering complete: 96399 junctions retained. ~90s
Writing IR table...
IR calculated for 6/6 samples
IR table written. ~150s
Done. Total runtime: ~150s
```



Actual output

```bash
Thu Jun 25 20:37:52 UTC 2026
Starting ir_table with 6 samples
Loading annotation...
Annotation loaded: 528735 annotated junctions. 55.5s
Gathering inclusion counts and clusters...
Loaded 6 samples and 333069 clusters. 59.2s
Collecting junctions across all samples...
getJunctions complete: 210470 junctions. 4.7s
RSD filtering complete: 96399 junctions retained. 24.7s
Junction collection and RSD filtering complete: 96399 junctions retained. 83.9s
Writing IR table...
IR calculated for 6/6 samples
IR table written. 145.5s
Done. Total runtime: 145.5s
real    2m30.017s
user    0m0.059s
sys     0m0.051s
Thu Jun 25 20:40:22 UTC 2026
```



## Compare output to known good

```bash
diff ${this_base_dir}/analysis/_intron_retention.tsv \
    <(zcat ${known_good_dir}/_intron_retention.tsv.gz) | head
```

Expected: no output (files identical).

Actual result:

no differences!



```bash
ubuntu@hbeale-mesa:/mnt/gitCode/splicedice$ diff ${this_base_dir}/analysis/_intron_retention.tsv \
    <(zcat ${known_good_dir}/_intron_retention.tsv.gz) | head
ubuntu@hbeale-mesa:/mnt/gitCode/splicedice$ head  ${this_base_dir}/analysis/_intron_retention.tsv
Junction        SRR12801019     SRR12801023     SRR12801024     SRR12801028     SRR12801027     SRR12801020
GL000008.2:163999-164602:-      0.200   nan     nan     nan     0.000   0.000
GL000009.2:55162-78082:-        0.100   0.000   0.000   0.000   nan     0.083
GL000194.1:112851-114985:-      0.000   0.000   0.000   0.000   0.012   0.005
GL000194.1:11336-20237:-        1.000   0.000   nan     0.000   0.000   nan
GL000194.1:11336-28266:-        0.333   0.111   nan     0.000   0.333   nan
GL000194.1:20374-28266:-        0.333   0.250   nan     0.000   1.000   nan
GL000194.1:54833-55445:-        0.406   0.277   0.348   0.282   0.378   0.667
GL000194.1:55677-112791:-       0.014   0.019   0.011   0.039   0.038   0.008
GL000195.1:138141-139989:+      0.024   0.029   0.014   0.006   0.029   0.009
ubuntu@hbeale-mesa:/mnt/gitCode/splicedice$ zcat ${known_good_dir}/_intron_retention.tsv.gz | head
Junction        SRR12801019     SRR12801023     SRR12801024     SRR12801028     SRR12801027     SRR12801020
GL000008.2:163999-164602:-      0.200   nan     nan     nan     0.000   0.000
GL000009.2:55162-78082:-        0.100   0.000   0.000   0.000   nan     0.083
GL000194.1:112851-114985:-      0.000   0.000   0.000   0.000   0.012   0.005
GL000194.1:11336-20237:-        1.000   0.000   nan     0.000   0.000   nan
GL000194.1:11336-28266:-        0.333   0.111   nan     0.000   0.333   nan
GL000194.1:20374-28266:-        0.333   0.250   nan     0.000   1.000   nan
GL000194.1:54833-55445:-        0.406   0.277   0.348   0.282   0.378   0.667
GL000194.1:55677-112791:-       0.014   0.019   0.011   0.039   0.038   0.008
GL000195.1:138141-139989:+      0.024   0.029   0.014   0.006   0.029   0.009
ubuntu@hbeale-mesa:/mnt/gitCode/splicedice$ 
```

