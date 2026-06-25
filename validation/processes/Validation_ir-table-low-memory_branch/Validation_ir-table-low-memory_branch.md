# Validation: ir-table-low-memory branch

Validates that the memory optimizations in `ir-table-low-memory` produce identical `_intron_retention.tsv` output to the known-good reference.

Changes validated:

- `getInclusionCounts`: sparse dict reader replacing pandas DataFrame pivot
- `getJunctions`: `imap_unordered` with incremental union replacing `pool.map`
- `getFilteredJunctions`: `imap_unordered` with incremental union replacing `pool.map`

# Server

ubuntu@hbeale-clin-validation

# Setup

## Define locations

```bash
this_base_dir=/mnt/sd/validation_ir_table_low_memory_2026.06.25_11.38.03/
known_good_dir=/mnt/sd/known_good/splicedice/validation/data/SUGP1_SRP286876
validation_source=/mnt/sd/validation_2b2c401_2026.06.23_15.07.59
genes=/mnt/ref/gencode.v47.primary_assembly.annotation.gtf

mkdir -p ${this_base_dir}/analysis/
```

## Get known-good reference data

(not necessary; it's already present)

```bash
mkdir -p /mnt/sd/known_good
cd /mnt/sd/known_good
git clone --branch validation \
    https://github.com/BrooksLabUCSC/splicedice.git
```

## Build Docker image from ir-table-low-memory branch

```bash
cd /mnt/git_code/splicedice
git checkout ir-table-low-memory

docker build --build-arg CACHE_BUST=$(date +%s) \
    -t splicedice_analysis:ir_table_low_memory \
    -f validation/processes/Validation_ir-table-low-memory_branch/Dockerfile_ir_table_low_memory.txt .
~/alert_msg.sh "docker build complete"
```

## prepare input

```bash
zcat ${known_good_dir}/_inclusionCounts.tsv.gz > ${this_base_dir}/analysis/_inclusionCounts.tsv
zcat ${known_good_dir}/_allClusters.tsv.gz > ${this_base_dir}/analysis/_allClusters.tsv
```



# Run ir_table

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

# Compare output to known good

```bash
diff ${this_base_dir}/analysis/_intron_retention.tsv \
    <(zcat ${known_good_dir}/_intron_retention.tsv.gz)
```

Expected: no output (files identical).

