#!/usr/bin/env bash
# run_intron_coverage.sh
# Runs splicedice intron_coverage on a batch manifest, then optionally removes
# BAMs for samples that produced output successfully.
#
# Usage: run_intron_coverage.sh --manifest PATH --coverage-dir PATH --analysis-base PATH
#                               --disk-constraint [yes|no] [--max-threads N]
set -euo pipefail

# ── defaults ──────────────────────────────────────────────────────────────────
manifest=""
coverage_dir=""
analysis_base=""
disk_constraint=""
max_threads=8

# ── argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case $1 in
        --manifest)        manifest="$2";        shift 2 ;;
        --coverage-dir)    coverage_dir="$2";    shift 2 ;;
        --analysis-base)   analysis_base="$2";   shift 2 ;;
        --disk-constraint) disk_constraint="$2"; shift 2 ;;
        --max-threads)     max_threads="$2";     shift 2 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "$manifest" || -z "$coverage_dir" || -z "$analysis_base" || -z "$disk_constraint" ]]; then
    echo "ERROR: --manifest, --coverage-dir, --analysis-base, and --disk-constraint are required"
    echo "Usage: $0 --manifest PATH --coverage-dir PATH --analysis-base PATH --disk-constraint [yes|no] [--max-threads N]"
    exit 1
fi

if [[ "$disk_constraint" != "yes" && "$disk_constraint" != "no" ]]; then
    echo "ERROR: --disk-constraint must be 'yes' or 'no'"
    exit 1
fi

mkdir -p "$coverage_dir"

# ── check BAMs are actually present before running ────────────────────────────
# Primary manifest columns: dataset_id  bam_location  bed_location  phenotype  download_id
available_manifest="${manifest%.tsv}_available.tsv"
> "$available_manifest"

while IFS=$'\t' read -r id bam_location rest; do
    if [[ -f "$bam_location" ]]; then
        printf '%s\t%s\t%s\n' "$id" "$bam_location" "$rest" >> "$available_manifest"    
    else
        echo "WARNING: $id: BAM not found at $bam_location — skipping"
    fi
done < <(tail -n +2 "$manifest")

n_samples=$(wc -l < "$available_manifest")

if [[ $n_samples -eq 0 ]]; then
    echo "No BAMs available; skipping intron_coverage."
    exit 0
fi

n_threads=$(( n_samples < max_threads ? n_samples : max_threads ))
echo "Running intron_coverage on ${n_samples} samples with ${n_threads} threads..."
date

time docker run --rm \
    -v /mnt/:/mnt \
    splicedice_analysis:latest \
    splicedice intron_coverage \
    -b "$available_manifest" \
    -j "${analysis_base}/_junctions.bed" \
    -n "${n_threads}" \
    -o "${coverage_dir}"

date
echo "intron_coverage complete."

# ── optionally clean up BAMs for samples with confirmed output ────────────────
if [[ "$disk_constraint" == "yes" ]]; then
    echo ""
    echo "--- disk_constraint=yes: cleaning up BAMs ---"
    while IFS=$'\t' read -r id bam_location _rest; do
        expected_output="${coverage_dir}/${id}_intron_coverage.txt"
        if [[ -f "$bam_location" && -f "$expected_output" ]]; then
            echo "$id: output confirmed, deleting BAM"
            rm "$bam_location"
        elif [[ -f "$bam_location" && ! -f "$expected_output" ]]; then
            echo "WARNING: $id: BAM present but no output found — keeping BAM for investigation"
        fi
    done < <(tail -n +2 "$available_manifest")
else
    echo "--- disk_constraint=no: keeping BAMs ---"
fi