# SpliceDICE
Splice Divergent Interval Co-Exclusion (Splice DICE) is a tool for detecting and quantifying splicing events
by mutually exclusive junctions. This tool is currently in development and
user discretion is advised.

## Table of Contents
  * [Dependencies](#dependencies)
  * [Installation](#installation)
  * [Usage](#usage)
    + [Aligned RNA sequencing reads](#aligned-rna-sequencing-reads)
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

SpliceDICE uses splice junction counts derived from aligned RNA sequencing reads to calculate a Percent-Spliced (PS) value for each junction.

**Suggested tool for generating junction counts:**

- [intronProspector](https://github.com/diekhans/intronProspector)

### Aligned RNA sequencing reads
Calculating intron coverage with SpliceDICE requires RNA sequencing reads that are aligned to a reference genome. 

## Manifest files

### `splicedice quant`
Processes junction count files (bed6 files) to calculate Percent-Spliced (PS) value for every splice junction in every sample in the manifest.
For information on the `bed_manifest.txt` format, see [Manifest Format](#manifest-format).

BED input format (per sample file):
- Tab-delimited BED6 with columns: `chrom`, `start`, `end`, `name`, `score`, `strand`.
- Genomic coordinates must use 0-based, half-open BED convention (UCSC standard).
- `score` is used as the junction read count.
- `strand` must be `+` or `-`.

Input parameters:
- `-m, --manifest` (required): tab-separated manifest file with sample names, bed file paths, and metadata.
- `-o, --output_prefix` (required): prefix used for all generated output files.
- `--drim` (optional): also write `{output_prefix}_drimTable.tsv` for DRIMSeq.

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
- `{output_prefix}_drimTable.tsv` (optional): created only when `--drim` is used.

### Intron Retention
The percent-spliced value does not quantify intron retention, so separate subprograms gives a table of IR values, in the same format as the PS table. The first subprogram, `splicedice intron_coverage`, measures the coverage across previously identified splice junctions, and outputs a table for each sample. The second subprogram, `splicedice ir_table`, takes those coverage values and calculates the IR value for each junction in each sample, outputting the final IR table.

```bash
$ splicedice intron_coverage -b bam_manifest.tsv -m project_allPS.tsv -j project_junctions.bed -n 4 -o coverage_output_dir
$ splicedice ir_table -i project_inclusionCounts.tsv -c project_allClusters.tsv -d coverage_output_dir -n 4 -o project_output_prefix
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

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please include tests as appropriate.

## License
[BSD-3](LICENSE)
