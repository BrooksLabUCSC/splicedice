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
```bash
$ mesa quant -m bed_manifest.tsv -g genome.fa -o my_output
```

#### Additional options
todo example sub arguments

### `mesa compare_sample_sets`

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


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](LICENSE)
