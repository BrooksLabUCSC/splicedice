
# Changelog

## [1.1.0]

### Changed
- `quant` input uses standard 6 column bed file
- Updated dependencies
    - Build uses `requirements.txt`
    - Added `pysam` to dependencies. Was previously missing
    - Pinned all tool versions

### Removed
- `bam_to_junc_bed.py`. Use diekhans/intronProspector instead.
- `quant` input support and cl parameters
    - Remove all junction filtering parameters
    - removes support for STAR SJ.out.tab files
    - removes support for `splicedice bam_to_junc_bed` bed files

## [1.0.0] - Mesa

Initial release.
