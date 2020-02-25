# MESA
MESA (**M**utually **E**xclusive **S**plicing **A**nalysis) is a tool for detecting and quantifying splicing events
by mutually exclusive junctions.  This tool is currently in development and
user discretion is advised.

## Dependencies
- STAR aligner
- bedtools

## Installation

Use `git` to clone the package and install with the package manager [pip](https://pip.pypa.io/en/stable/).

```bash
$ git clone https://github.com/BrooksLabUCSC/mesa.git
$ cd mesa/
$ pip install --user .
```
### Development
If you are working on developing MESA it will likely be useful to install it in editable mode.
```bash
$ pip install --user -e .
```

## Usage

### Alignment with star
MESA uses output files from the [STAR](https://github.com/alexdobin/STAR)
RNA-seq aligner. An example of using star aligner is shown below.
```bash
todo example if needed
```

### `mesa star_junc_to_bed`
Converts star aligner tab output files into bed files
```bash
$ mesa star_junc_to_bed -s sj.tab
```

### `mesa quant`
Processes STAR junction bed files and genome to produce output for later steps.
For information on the `bed_manifest.txt` format, see [Manifest Format](#manifest-format).
```bash
$ mesa quant -m bed_manifest.txt -g genome.fa -o my_output
```

#### Additional options
todo example sub arguments

#### Output files
Based on the `-o/--output_prefix` parameter, `mesa quant` will output a number
of files for further processing.
- `{output_prefix}_allPSI.npz`: a numpy output file containing various
  information from `mesa quant`
- `{output_prefix}_all_clusters2.tsv`
- `{output_prefix}_inclusionCounts.tsv`
- `{output_prefix}_junctions.bed`
- `order.npz`
- `{output_prefix}_drimTable.tsv`: If the `--drim-table` option was passed, it
  will output its data in a format for use with DRIMSeq. See [Analyzing DRIMSeq
  output](#analyzing-drimseq-output) for how to analyze this data with a
  provided script.

### `mesa compare_sample_sets`
Compares the differences in splicing across two groups. Requires atleast 3
samples per condition, otherwise it will fail. If you have less than 3 samples
per condition, use `mesa parwise`.
```bash
$ mesa compare_sample_sets --psiMESA my_output_allPSI.npz -m1 ctrl_manifest.txt -m2 mut_manifest.txt
```
This module will take two manifest files, that represents the two groups you
wish to compare as well as the allPSI.npz file that was previously outputted by
`mesa quant`. TODO what is the output

### `mesa pairwise`

### `mesa cluster`

## Manifest Format
The manifest is a tab-delimitted file used by `mesa` provides information about
the samples and related files.
```bash
$ cat manifest.tsv
sample1 /path/to/sample1/sj.tab.bed lung control
sample2 /path/to/sample2/sj.tab.bed lung control
sample3 /path/to/sample3/sj.tab.bed lung mutant
sample4 /path/to/sample4/sj.tab.bed lung mutant
```
- The first column is the sample identifier.
- The second column is the absolute path to the bed file version of the star junction output file produced by `mesa star_junc_to_bed`
- The third column is additional metadata for the type of sample it is. This column is for convience for your own analyses and not used by `mesa`.
- The fourth column is the condition. This is used to decide how the samples are
grouped and the statistical analysis uses the different groups to compare.

An example of the manifest format can be found [here](data/example_manifest.tsv).

### Analyzing DRIMSeq output
`mesa quant` can provide its output in a format for use with with the
alternative splicing quantifier tool DRIMSeq in the R programming language.
MESA provides a utility script [run-drim-seq.R](scripts/run-drim-seq.R). This
script depends on:
- R todo double check dependencies
- DRIMSeq
- ggplot2
- optparse

An example command to run our DRIMSeq script is shown below
```
$ Rscript run-drim-seq.R -m bed_manifest.tsv -d drim_table.tsv -o drim_output -t 12
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
[MIT](LICENSE)
