# End-to-End Tests

This folder contains end-to-end tests for project components. Each subdirectory provides a test
suite for a specific component.

Available test suites

- **quant**: tests for the quantification pipeline (located in `quant/`).

## quant

### exon_skip.sh

Runs validation using the toy dataset for the quant pipeline.

What this test verifies

- **Inputs accepted:** expected input files are handled correctly.
- **Output format & content:** output tables are produced with the expected format and values.
- **Indexing:** output uses 1-based indexing where applicable.

Run the test from the project root:

```bash
bash tests/e2e/quant/exon_skip.sh
```

About the dataset

The dataset is a small toy example representing a single exon-skipping event. Three samples
represent different percent-spliced (PS) values for the event; a fourth sample contains an
unrelated junction used as a control.

Inputs and expected outputs for the test are under `data/example_data/exon_skip/`. See
`data/example_data/exon_skip/README.md` for dataset details and expected output files.