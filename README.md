# SpliceDICE
Splice Divergent Interval Co-Exclusion (Splice DICE) is a tool for detecting and quantifying splicing events
by mutually exclusive junctions. This tool is currently in development and
user discretion is advised.

## Table of Contents
  * [Dependencies](#dependencies)
  * [Installation](#installation)
  * [Usage](#usage)
    + [Aligned RNA sequencing reads](#aligned-rna-sequencing-reads)
    + [`splicedice bam_to_junc_bed`](#splicedice-bam_to_junc_bed)
    + [`splicedice quant`](#splicedice-quant)
      - [Output files](#output-files)
    + [`splicedice compare_sample_sets`](#splicedice-compare_sample_sets)
    + [`splicedice pairwise`](#splicedice-pairwise)
    + [Intron Retention](#intron-retention)

  * [Manifest Format](#manifest-format)
  * [Analyzing DRIMSeq output](#analyzing-drimseq-output)
  * [Contributing](#contributing)
  * [License](#license)

## Dependencies
- python=3.7+
- numpy
- samtools
- pysam
- scipy
- ...

## Installation

Use `git` to clone the package and install with the package manager [pip](https://pip.pypa.io/en/stable/).

```bash
$ git clone https://github.com/BrooksLabUCSC/splicedice.git
$ cd splicedice/
$ pip install --user .
```
### Development
If you are working on developing SpliceDICE it will likely be useful to install it in editable mode.
```bash
$ pip install --user -e .
```

## Usage

SpliceDICE uses counts of splice junctions from aligned RNA sequencing reads, to calculate a Percent-Spliced (PS) value for each junction. It performs best when the junction counts are gathered from a BAM file using `splicedice bam_to_junc_bed`, but can accept

### Aligned RNA sequencing reads
SpliceDICE requires RNA sequencing reads that are aligned to a reference genome. 

## Manifest files


### `splicedice bam_to_junc_bed`
Searches aligned RNA-seq reads (BAM files) for splice junctions, and outputs a bed file with junction counts. Takes a manifest with file paths, and outputs a .junc.bed file for each BAM.
```bash
$ splicedice bam_to_junc_bed -m bam_manifest.txt
```

### `splicedice quant`
Processes junction count files (bed files from `splicedice bam_to_junc_bed` or SJ.out.tab from STAR aligner) to calculate Percent-Spliced (PS) value for every splice junction in every sample in the manifest.
For information on the `bed_manifest.txt` format, see [Manifest Format](#manifest-format).
```bash
$ splicedice quant -m bed_manifest.txt -o output_prefix
```

#### Output files
Based on the `-o/--output_prefix` parameter, `splicedice quant` will output a number
of files for further processing.
- `{output_prefix}_allPS.tsv`: A tab-separated table of PS values, where each column is a sample, and each row is a splice junction.
- `{output_prefix}_allClusters.tsv`
- `{output_prefix}_inclusionCounts.tsv`
- `{output_prefix}_junctions.bed`

### `splicedice compare_sample_sets`
Compares the differences in splicing across two groups. Requires atleast 3
samples per condition, otherwise it will fail. If you have less than 3 samples
per condition, use `splicedice parwise`.
```bash
$ splicedice compare_sample_sets --psiSPLICEDICE my_output_allPSI.npz -m1 ctrl_manifest.txt -m2 mut_manifest.txt
```
This module will take two manifest files, that represents the two groups you
wish to compare. It also takes the allPS.tsv file that was previously outputted by
`splicedice quant`.

### `splicedice pairwise`
Example command:

```bash
$ splicedice pairwise --inclusionSPLICEDICE my_output_inclusionCounts.npz -c my_output_all_clusters2.tsv >pairwise_output.txt
```

This command performs a pairwise comparison of the junction usage for each
sample against each other. `pairwise_output.txt`, from the command above, will
output the p-value from running a Fisher's exact test, where the contingency
matrix is inclusion and exclusion read counts for each pair of samples for a
given junction. This command is recommend for datasets with less than three
samples per group where `splicedice compare_sample_sets` could not be used.

todo example output and explanation

### Intron Retention
The percent-spliced value does not quantify intron retention, so separate subprograms gives a table of IR values, in the same format as the PS table. The first subprogram, `splicedice intron_coverage`, measures the coverage across previously identified splice junctions, and outputs a table for each sample. The second subprogram, `splicedice ir_table`, takes those coverage values and calculates the IR value for each junction in each sample, outputting the final IR table.

```bash
$ splicedice intron_coverage -b bam_manifest.tsv -m project_allPS.tsv -j project_junctions.bed -n 4 -o coverage_output_dir
$ splicedice ir_table -i project_inclusionCounts.tsv -c project_allClusters.tsv -d coverage_output_dir -o project_output_prefix
```


## Manifest Format
The manifest is a tab-delimitted file used by `splicedice` provides information about
the samples and related files.
```bash
$ cat manifest.txt
sample1 /path/to/sample1/sj.tab.bed lung control
sample2 /path/to/sample2/sj.tab.bed lung control
sample3 /path/to/sample3/sj.tab.bed lung mutant
sample4 /path/to/sample4/sj.tab.bed lung mutant
```
- The first column is the sample identifier.
- The second column is the absolute path to the bed file version of the star junction output file produced by `splicedice star_junc_to_bed`
- The third column is additional metadata for the type of sample it is. This column is for convience for your own analyses and not used by `splicedice`.
- The fourth column is the condition. This is used to decide how the samples are
grouped and the statistical analysis uses the different groups to compare.

An example of the manifest format can be found [here](data/example_manifest.txt).

## Analyzing DRIMSeq output
`splicedice quant` can provide its output in a format for use with with the
alternative splicing quantifier tool DRIMSeq in the R programming language.
SpliceDICE provides a utility script [run-drim-seq.R](scripts/run-drim-seq.R). This
script depends on:
- R todo double check dependencies
- DRIMSeq
- ggplot2
- optparse

An example command to run our DRIMSeq script is shown below
```
$ Rscript run-drim-seq.R -m bed_manifest.txt -d drim_table.tsv -o drim_output -t 12
```

TODO explain various plots and data output

This script will automatically output a number of plots and an output table for
further analyzing splicing data. The `-t` parameter sets the number of threads
to be used, and we recommend setting it as high as you can because DRIMSeq is a
cpu-intensive tool.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[BSD-3](LICENSE)
