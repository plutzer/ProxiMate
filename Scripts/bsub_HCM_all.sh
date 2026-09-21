#!/bin/bash
set -euo pipefail

# Submit HCM parameter-sweep jobs against proximate:v0.1.10.
# Runs 7 combinations of quant type x imputation mode.
# Run this from the directory that contains ProxiMate_test_datasets/.

IMAGE="plutzer/proximate:v0.1.10"
HCM_DATA="$(pwd)/ProxiMate_test_datasets/HCM"
OUTPUT_BASE="$(pwd)/HCM_All"
LOG_DIR="${OUTPUT_BASE}/logs"
ERR_DIR="${OUTPUT_BASE}/errors"
JOB_GROUP="/plutzer/proximate_test"
N_ITER=1000

mkdir -p "$LOG_DIR" "$ERR_DIR"

submit_run() {
    local name="$1"
    local quant="$2"
    local imputation="$3"

    local run_dir="${OUTPUT_BASE}/${name}"
    mkdir -p "$run_dir"

    LSF_DOCKER_VOLUMES="${HCM_DATA}:/testdata ${run_dir}:/output" \
        bsub \
        -G compute-bmajor \
        -g "$JOB_GROUP" \
        -q general \
        -J "HCM_${name}" \
        -oo "${LOG_DIR}/${name}.txt" \
        -eo "${ERR_DIR}/${name}.txt" \
        -a "docker(${IMAGE})" \
        bash /run_pipeline.sh --format maxquant \
            /testdata/HCM_ED_corrected_replicates.csv \
            /testdata/proteinGroups.txt \
            "$quant" /output "$N_ITER" "$imputation"

    echo "  Submitted: HCM_${name}  (quant=${quant}, imputation=${imputation})"
}

echo "Submitting HCM parameter sweep..."
echo "  Image:  $IMAGE"
echo "  Output: $OUTPUT_BASE"
echo ""

# Imputation codes (from run_pipeline.sh):
#   0 = none, 2 = two-component (refactored) AFT, 3 = one-component AFT
submit_run "SpectralCounts_noimp"   "Spectral Counts"  0
submit_run "LFQ_noimp"              "LFQ"              0
submit_run "LFQ_aft_one_component"  "LFQ"              3
submit_run "LFQ_aft_two_component"  "LFQ"              2
submit_run "Intensity_noimp"        "Intensity"        0
submit_run "Intensity_aft_one_component"  "Intensity"  3
submit_run "Intensity_aft_two_component"  "Intensity"  2

echo ""
echo "All jobs submitted."
