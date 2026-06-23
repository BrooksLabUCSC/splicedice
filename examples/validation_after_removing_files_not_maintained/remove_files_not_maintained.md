# Removing files that are not used or not maintained

Note: if you are using this file as an example for running splicedice, it contains a lot of references to the splicedice_analysis repo. Consider instead using SUGP1_walkthrough_2026.06.23_11.24.14 .md as an example. 

# Summary

There are many files that are not used or not maintained in the repo, and I think this causes confusion. 

We can add a note to the readme that additional scripts were removed and can be viewed by going to commit X, and linking the commit. 

I propose the following, keeping only the functionality supporting quant, intron_coverage, and ir_table



## directories and top level files

### remove 



#### `scripts/` (entire directory)

It contains only `run-drim-seq.R` — the DRIMSeq integration script. You don't use DRIMSeq, so this is safe to remove. Note that `setup.py` references it via `script=["scripts/run-drim-seq.R"]`, but that line is a minor packaging artifact and removing the file won't break `pip install -e .` or any of the three subcommands you use.

#### `data/` (entire directory)

Example/test data files. Not needed at runtime.

#### `docs/` (entire directory)

Documentation only.

#### `tests/` (entire directory)

The only tests are under `tests/e2e/quant/`. Safe to remove since you're not running the test suite.

#### `Pipfile`

Only needed for `pipenv`. Since you use `pip`, irrelevant.

#### `MANIFEST.in`

Only needed for building a PyPI distribution package. Not needed for `pip install -e .`.



### keep

| Item                      | Why                                                          |
| ------------------------- | ------------------------------------------------------------ |
| `setup.py`                | Wires up the `splicedice` CLI entry point                    |
| `requirements.txt`        | Read by `setup.py` at install time                           |
| `README.md`               | Read by `setup.py` for `long_description` — removing it would cause `pip install` to crash |
| `.gitignore`              | Harmless, but irrelevant to function                         |
| `LICENSE`, `CHANGELOG.md` | Tiny, harmless                                               |

## in the splicedice code directory

### remove

```bash
splicedice/compareSampleSets.py   # compare_sample_sets
splicedice/pairwise_fisher.py     # pairwise
splicedice/findOutliers.py        # findOutliers
splicedice/subset.py              # subset
splicedice/similarity.py          # similarity
splicedice/select_samples.py      # select
splicedice/counts_to_ps.py        # counts_to_ps
```

### Keep

```bash
splicedice/__init__.py            # package init, always required
splicedice/__main__.py            # CLI entry point
splicedice/SPLICEDICE.py          # quant
splicedice/intron_coverage.py     # intron_coverage
splicedice/ir_table.py            # ir_table
```

I checked whether any of the "keep" modules import from the "remove" modules; they do not.

```bash
grep -l "import compareSampleSets\|import pairwise_fisher\|import findOutliers\|import subset\|import similarity\|import select_samples\|import counts_to_ps" \
  splicedice/SPLICEDICE.py splicedice/intron_coverage.py splicedice/ir_table.py
```



## Edits

change SPLICEDICE.py and `__main__.py` to remove references to other modules



# Test

(Note, this is the second test; the first one is documented in `remove_files_not_maintained v0 - obsolete.md`

server:

ubuntu@hbeale-mesa



# Setup per run

## define location

```bash
this_commit=4ba8c72
this_description=test_Remove-untested-code
this_datestamp=2026.06.23_11.14.13
this_branch=Remove-untested-code
this_dockerfile=Dockerfile_${this_description}_${this_commit}.txt
```

```bash

this_base_dir=/mnt/sd/${this_description}_${this_commit}_${this_datestamp}/
code_base=${this_base_dir}/git_code/splicedice_analysis/examples/test-Remove-untested-code
mkdir -p ${this_base_dir}/git_code/ ${this_base_dir}/analysis/

```



## get code

```bash
cd ${this_base_dir}/git_code/
git clone https://github.com/hbeale/splicedice_analysis.git
```

## build docker

```bash
cd ${code_base}
cat ${this_base_dir}/git_code/splicedice_analysis/2026_06_ps_ir_pipeline/Dockerfile_by_branch.txt | sed "s/replace_with_branch/$this_branch/" > ${code_base}/${this_dockerfile}

docker build --build-arg CACHE_BUST=$(date +%s) -t splicedice_analysis:latest -f ${code_base}/${this_dockerfile} .
bash ~/alert_msg.sh "docker build complete"

```

## identify manifests

```bash
bed_manifest=${this_base_dir}/git_code/splicedice_analysis/SUGP1/bed_manifest_2026-05-05_20-28-22.tsv
bam_manifest=${this_base_dir}/git_code/splicedice_analysis/SUGP1/bam_manifest.tsv 

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





# Run pipeline

## run quant

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

aside: ~/alert_msg.sh is a convenience script that notifies me when the run is complete

```bash
Tue Jun 23 18:23:08 UTC 2026
Parsing manifest...
        Done [0:00:0.00]
Getting all junctions from 6 files...
        Done [0:00:4.45]
Finding clusters from 333069 junctions...
        Done [0:00:3.43]
Writing cluster file...
        Done [0:00:3.66]
Writing junction bed file...
        Done [0:00:2.10]
Gathering junction counts...
        Done [0:00:5.65]
Writing inclusion counts...
        Done [0:00:5.34]
Calculating PS values...
        Done [0:00:14.01]
Writing PS values...
        Done [0:00:5.47]
All done [0:00:44.11]

real    0m48.101s
user    0m0.048s
sys     0m0.060s
Tue Jun 23 18:23:56 UTC 2026
```

## run intron coverage

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
Tue Jun 23 18:28:34 UTC 2026
getting paths for bam files
creating junction percentiles

```



## Create intron table

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



std out

```bash

Tue Jun 23 18:28:34 UTC 2026   
getting paths for bam files
creating junction percentiles  
SRR12801019 starting 4.838646411895752
SRR12801019 collected 608.4695711135864
SRR12801019 counted 1305.9411807060242
SRR12801019 done 1321.6828591823578
SRR12801023 starting 6.312715768814087
SRR12801023 collected 645.777664899826
SRR12801023 counted 1381.3058812618256
SRR12801023 done 1397.4367985725403
SRR12801024 starting 7.037007808685303
SRR12801024 collected 771.2980704307556
SRR12801024 counted 1529.4202527999878
SRR12801024 done 1544.9719219207764
SRR12801028 starting 8.533150672912598
SRR12801028 collected 800.1043570041656
SRR12801028 counted 1619.3926177024841
SRR12801028 done 1634.7489924430847
SRR12801027 starting 7.764169692993164
SRR12801027 collected 762.325615644455
SRR12801027 counted 1622.007576227188
SRR12801027 done 1637.4010829925537
Your runtime was 1654.5064775943756 seconds.

real    27m39.771s
user    0m0.164s
sys     0m0.056s
Tue Jun 23 18:56:14 UTC 2026

```



# Compare to previous output

```bash
prev_dir=/mnt/splicedice_ir_example_archives/2026.05.28_21.25.54/analysis/
current_dir=${this_base_dir}/analysis
for i in `ls $current_dir`; do echo $i; diff $current_dir/$i $prev_dir/${i}; done
```

```bash
ubuntu@hbeale-mesa:/mnt/sd/test_Remove-untested-code_4ba8c72_2026.06.23_11.14.13$ for i in `ls $current_dir`; do echo $i; diff $current_dir/$i $prev_dir/${i}; done
_allClusters.tsv
_allPS.tsv
_inclusionCounts.tsv
_junctions.bed
ubuntu@hbeale-mesa:/mnt/sd/test_Remove-untested-code_4ba8c72_2026.06.23_11.14.13$ 

```

make sure my code works:

```bash
echo 1 > $prev_dir/not_matching
echo 2 > $current_dir/not_matching
for i in `ls $current_dir`; do echo $i; diff $current_dir/$i $prev_dir/${i}; done
```

```bash
ubuntu@hbeale-mesa:/mnt/sd/test_Remove-untested-code_4ba8c72_2026.06.23_11.14.13$ echo 1 > $prev_dir/not_matching
echo 2 > $current_dir/not_matching
for i in `ls $current_dir`; do echo $i; diff $current_dir/$i $prev_dir/${i}; done
_allClusters.tsv
_allPS.tsv
_inclusionCounts.tsv
_junctions.bed
not_matching
1c1
< 2
---
> 1

```

