splicedice quant \
    --manifest=data/example_data/exon_skip/manifest.tsv \
    --drim \
    --output_prefix=data/output/

failures=0

run_diff_test() {
    test_name="$1"
    actual_file="$2"
    expected_file="$3"

    if diff -u "$actual_file" "$expected_file" >/dev/null; then
        echo "PASS: $test_name"
    else
        echo "FAIL: $test_name"
        diff -u "$actual_file" "$expected_file"
        failures=$((failures + 1))
    fi
}

run_diff_test "allClusters" data/output/_allClusters.tsv data/example_data/exon_skip/expected_out/_allClusters.tsv
run_diff_test "allPS" data/output/_allPS.tsv data/example_data/exon_skip/expected_out/_allPS.tsv
run_diff_test "inclusionCounts" data/output/_inclusionCounts.tsv data/example_data/exon_skip/expected_out/_inclusionCounts.tsv
run_diff_test "junctions.bed" data/output/_junctions.bed data/example_data/exon_skip/expected_out/_junctions.bed 
run_diff_test "drimTable" data/output/_drimTable.tsv data/example_data/exon_skip/expected_out/_drimTable.tsv

if [ "$failures" -eq 0 ]; then
    echo "PASS: all quant e2e checks"
else
    echo "FAIL: $failures quant e2e check(s) failed"
    exit 1
fi