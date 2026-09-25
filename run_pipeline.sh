#!/bin/bash

# A stage's exit status must survive the `| tee` below, or a failed parse would
# be reported as success and scoring would run on absent inputs.
set -o pipefail

# Kept for the run header: the option parsing below shifts these away.
invocation="$0 $*"

# Usage function
usage() {
    echo "Usage:"
    echo "  $0 [OPTIONS] --format maxquant  ED_file PG_file quant_type output_dir n_iterations imputation"
    echo "  $0 [OPTIONS] --format diann     ED_file matrix_file output_dir n_iterations imputation"
    echo "  $0 [OPTIONS] --format pioneer   ED_file protein_groups_wide.tsv output_dir n_iterations imputation"
    echo "  $0 [OPTIONS] --format fragpipe  ED_file FP_file quant_type output_dir n_iterations imputation"
    echo "  $0 [OPTIONS] --format msstats   ED_file ProteinLevelData.csv output_dir n_iterations imputation"
    echo "  $0 [OPTIONS] --format saint     bait_file prey_file interaction_file quant_type output_dir n_iterations imputation"
    echo ""
    echo "Formats:"
    echo "  maxquant  - MaxQuant proteinGroups.txt + experimental design CSV"
    echo "  diann     - DIA-NN report.pg_matrix.tsv + experimental design CSV"
    echo "  pioneer   - Pioneer protein_groups_wide.tsv + experimental design CSV"
    echo "  fragpipe  - FragPipe combined_protein.tsv + experimental design CSV"
    echo "  msstats   - MSstats ProteinLevelData.csv (already log2/normalized/imputed) + experimental design CSV"
    echo "  saint     - SAINT bait.txt, prey.txt, interaction.txt"
    echo ""
    echo "Options:"
    echo "  --organism    - human (default), mouse, or yeast"
    echo "  --exclude-hcm - annotate against BioGRID with Human Cell Map (Go et al. 2021) evidence removed; human only"
    echo "  --pi-method   - weighted_average (default) or single_bait; applies when imputation=2"
    echo "  --pi-bait     - required when --pi-method=single_bait: control Bait name"
    echo "  --aft-min-obs - preys observed in fewer than N runs take the dataset median per-prey SD"
    echo "                  as their sigma lower bound (imputation 2 or 3); default 0 (off)"
    echo "  --seed        - CompPASS permutation seed, so WD p-values reproduce"
    echo ""
    echo "Arguments:"
    echo "  quant_type    - Intensity, LFQ, or 'Spectral Counts'"
    echo "  n_iterations  - Number of CompPASS resampling iterations"
    echo "  imputation    - 0 (none), 1 (prey-specific AFT), 2 (refactored AFT), or 3 (one-component AFT)"
    echo ""
    echo "Output written to output_dir:"
    echo "  proximate.log - timestamped log from each stage, with the run ID"
    echo "  run.json      - parameters, input checksums and row counts for the run"
    echo "  log.txt       - raw console transcript, including output no logger sees"
    echo ""
    echo "The pipeline stops at the first stage that fails."
    exit 1
}

# Parse optional top-level flags in any order
organism="human"
pi_method="weighted_average"
pi_bait=""
aft_min_obs="0"
seed=""
hcm_args=()
while true; do
    case "$1" in
        --organism)   organism="$2"; shift 2 ;;
        --pi-method)  pi_method="$2"; shift 2 ;;
        --pi-bait)    pi_bait="$2"; shift 2 ;;
        --aft-min-obs) aft_min_obs="$2"; shift 2 ;;
        --seed)       seed="$2"; shift 2 ;;
        --exclude-hcm) hcm_args=(--excludeHCM); shift ;;
        *) break ;;
    esac
done

# Reusable arg arrays threaded into each score.py invocation
pi_args=(--pi-method "$pi_method")
[ -n "$pi_bait" ] && pi_args+=(--pi-bait "$pi_bait")
aft_args=(--aft-min-obs "$aft_min_obs")
seed_args=()
[ -n "$seed" ] && seed_args=(--seed "$seed")

# Check for --format flag
if [ "$1" != "--format" ] || [ -z "$2" ]; then
    usage
fi

format=$2
shift 2

# Every stage of one invocation shares a run ID, so their log lines and their
# run.json entries can be tied together.  Inherited if the caller set one.
if [ -z "${PROXIMATE_RUN_ID:-}" ]; then
    PROXIMATE_RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)-$(head -c4 /dev/urandom | od -An -tx1 | tr -d ' \n')"
fi
export PROXIMATE_RUN_ID

# Collect the per-format arguments; the stages themselves are shared below.
case "$format" in
    maxquant)
        [ "$#" -ne 6 ] && { echo "Error: maxquant format requires 6 arguments"; usage; }
        ed_file=$1; pg_file=$2; quant=$3; output_dir=$4; niters=$5; imp=$6
        stage_label="MaxQuant"
        parse_args=(--experimentalDesign "$ed_file" --proteinGroups "$pg_file" --quantType "$quant")
        score_ed="$ed_file"
        ;;

    fragpipe)
        [ "$#" -ne 6 ] && { echo "Error: fragpipe format requires 6 arguments"; usage; }
        ed_file=$1; fp_file=$2; quant=$3; output_dir=$4; niters=$5; imp=$6
        stage_label="FragPipe"
        parse_args=(--experimentalDesign "$ed_file" --fragpipeFile "$fp_file" --quantType "$quant")
        score_ed="$ed_file"
        ;;

    diann)
        [ "$#" -ne 5 ] && { echo "Error: diann format requires 5 arguments"; usage; }
        ed_file=$1; matrix_file=$2; output_dir=$3; niters=$4; imp=$5; quant="Intensity"
        stage_label="DIA-NN"
        parse_args=(--experimentalDesign "$ed_file" --diannMatrix "$matrix_file" --quantType "$quant")
        score_ed="$ed_file"
        ;;

    pioneer)
        [ "$#" -ne 5 ] && { echo "Error: pioneer format requires 5 arguments"; usage; }
        ed_file=$1; matrix_file=$2; output_dir=$3; niters=$4; imp=$5; quant="Intensity"
        stage_label="Pioneer"
        parse_args=(--experimentalDesign "$ed_file" --pioneerMatrix "$matrix_file" --quantType "$quant")
        score_ed="$ed_file"
        ;;

    msstats)
        [ "$#" -ne 5 ] && { echo "Error: msstats format requires 5 arguments"; usage; }
        ed_file=$1; msstats_file=$2; output_dir=$3; niters=$4; imp=$5; quant="Intensity"
        stage_label="MSstats"
        parse_args=(--experimentalDesign "$ed_file" --msstatsFile "$msstats_file" --quantType "$quant")
        score_ed="$ed_file"
        ;;

    saint)
        [ "$#" -ne 7 ] && { echo "Error: saint format requires 7 arguments"; usage; }
        bait_file=$1; prey_file=$2; interaction_file=$3
        quant=$4; output_dir=$5; niters=$6; imp=$7
        stage_label="SAINT"
        parse_args=(--bait "$bait_file" --prey "$prey_file" \
                    --interaction "$interaction_file" --quantType "$quant")
        # parse.py reconstructs the experimental design in the output directory.
        score_ed="$output_dir/ED.csv"
        ;;

    *)
        echo "Error: Unknown format '$format'. Must be maxquant, diann, pioneer, fragpipe, msstats, or saint."
        usage
        ;;
esac

mkdir -p "$output_dir"
log_file="$output_dir/log.txt"
export PROXIMATE_LOG_DIR="${PROXIMATE_LOG_DIR:-$output_dir}"

# run_stage <label> <command...>
# Tees the banner as well as the output, so log.txt shows which stage produced
# what, and aborts the pipeline when the stage fails.
run_stage() {
    local label="$1"; shift
    echo "=== $label ===" | tee -a "$log_file"
    "$@" 2>&1 | tee -a "$log_file"
    local status=${PIPESTATUS[0]}
    if [ "$status" -ne 0 ]; then
        echo "ERROR: $label failed with exit code $status; stopping." | tee -a "$log_file"
        exit "$status"
    fi
}

# log.txt is appended across runs, so mark where this one starts.
{
    echo ""
    echo "############################################################"
    echo "# ProxiMate run $PROXIMATE_RUN_ID"
    echo "# Started:  $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "# Command:  $invocation"
    echo "# Format:   $format    Organism: $organism"
    echo "############################################################"
} | tee -a "$log_file"

if [ "$format" = "msstats" ] && [ "$imp" != "0" ]; then
    echo "WARNING: MSstats data is already imputed by MBimpute; running AFT (imputation=$imp) will re-fit on imputed values. Pass imputation=0 to skip." | tee -a "$log_file"
fi

run_stage "Parsing ($stage_label)" \
    python3 /Scripts/parse.py "${parse_args[@]}" --outputPath "$output_dir"

run_stage "Scoring" \
    python3 /Scripts/score.py \
        --experimentalDesign "$score_ed" \
        --scoreInputs "$output_dir" \
        --outputPath "$output_dir" \
        --n-iterations "$niters" \
        --imputation "$imp" \
        --quantType "$quant" "${pi_args[@]}" "${aft_args[@]}" "${seed_args[@]}"

run_stage "Annotating" \
    python3 /Scripts/annotator.py \
        --organism "$organism" \
        --scoreFile "$output_dir/merged.csv" \
        --outputDir "$output_dir" "${hcm_args[@]}"

echo "=== Done (run $PROXIMATE_RUN_ID) ===" | tee -a "$log_file"
