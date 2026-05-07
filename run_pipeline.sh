#!/bin/bash

# Usage function
usage() {
    echo "Usage:"
    echo "  $0 [--organism ORG] --format maxquant  ED_file PG_file quant_type output_dir n_iterations imputation"
    echo "  $0 [--organism ORG] --format diann     ED_file matrix_file output_dir n_iterations imputation"
    echo "  $0 [--organism ORG] --format fragpipe  ED_file FP_file quant_type output_dir n_iterations imputation"
    echo "  $0 [--organism ORG] --format msstats   ED_file ProteinLevelData.csv output_dir n_iterations imputation"
    echo "  $0 [--organism ORG] --format saint     bait_file prey_file interaction_file quant_type output_dir n_iterations imputation"
    echo ""
    echo "Formats:"
    echo "  maxquant  - MaxQuant proteinGroups.txt + experimental design CSV"
    echo "  diann     - DIA-NN report.pg_matrix.tsv + experimental design CSV"
    echo "  fragpipe  - FragPipe combined_protein.tsv + experimental design CSV"
    echo "  msstats   - MSstats ProteinLevelData.csv (already log2/normalized/imputed) + experimental design CSV"
    echo "  saint     - SAINT bait.txt, prey.txt, interaction.txt"
    echo ""
    echo "Arguments:"
    echo "  --organism    - human (default), mouse, or yeast"
    echo "  --pi-method   - weighted_average (default) or single_bait; applies when imputation=2"
    echo "  --pi-bait     - required when --pi-method=single_bait: control Bait name"
    echo "  quant_type    - Intensity, LFQ, or 'Spectral Counts'"
    echo "  n_iterations  - Number of CompPASS resampling iterations"
    echo "  imputation    - 0 (none), 1 (prey-specific AFT), or 2 (refactored AFT)"
    exit 1
}

# Parse optional top-level flags (--organism, --pi-method, --pi-bait) in any order
organism="human"
pi_method="weighted_average"
pi_bait=""
while true; do
    case "$1" in
        --organism)   organism="$2"; shift 2 ;;
        --pi-method)  pi_method="$2"; shift 2 ;;
        --pi-bait)    pi_bait="$2"; shift 2 ;;
        *) break ;;
    esac
done

# Reusable arg array threaded into each score.py invocation
pi_args=(--pi-method "$pi_method")
[ -n "$pi_bait" ] && pi_args+=(--pi-bait "$pi_bait")

# Check for --format flag
if [ "$1" != "--format" ] || [ -z "$2" ]; then
    usage
fi

format=$2
shift 2

case "$format" in
    maxquant)
        if [ "$#" -ne 6 ]; then
            echo "Error: maxquant format requires 6 arguments"
            usage
        fi
        ed_file=$1
        pg_file=$2
        quant=$3
        output_dir=$4
        niters=$5
        imp=$6

        mkdir -p "$output_dir"

        echo "=== Parsing (MaxQuant) ==="
        python3 /Scripts/parse.py \
            --experimentalDesign "$ed_file" \
            --proteinGroups "$pg_file" \
            --quantType "$quant" \
            --outputPath "$output_dir" 2>&1 | tee -a "$output_dir/log.txt"

        echo "=== Scoring ==="
        python3 /Scripts/score.py \
            --experimentalDesign "$ed_file" \
            --scoreInputs "$output_dir" \
            --outputPath "$output_dir" \
            --n-iterations "$niters" \
            --imputation "$imp" \
            --quantType "$quant" "${pi_args[@]}" 2>&1 | tee -a "$output_dir/log.txt"

        echo "=== Annotating ==="
        python3 /Scripts/annotator.py \
            --organism "$organism" \
            --scoreFile "$output_dir/merged.csv" \
            --outputDir "$output_dir" 2>&1 | tee -a "$output_dir/log.txt"
        ;;

    fragpipe)
        if [ "$#" -ne 6 ]; then
            echo "Error: fragpipe format requires 6 arguments"
            usage
        fi
        ed_file=$1
        fp_file=$2
        quant=$3
        output_dir=$4
        niters=$5
        imp=$6

        mkdir -p "$output_dir"

        echo "=== Parsing (FragPipe) ==="
        python3 /Scripts/parse.py \
            --experimentalDesign "$ed_file" \
            --fragpipeFile "$fp_file" \
            --quantType "$quant" \
            --outputPath "$output_dir" 2>&1 | tee -a "$output_dir/log.txt"

        echo "=== Scoring ==="
        python3 /Scripts/score.py \
            --experimentalDesign "$ed_file" \
            --scoreInputs "$output_dir" \
            --outputPath "$output_dir" \
            --n-iterations "$niters" \
            --imputation "$imp" \
            --quantType "$quant" "${pi_args[@]}" 2>&1 | tee -a "$output_dir/log.txt"

        echo "=== Annotating ==="
        python3 /Scripts/annotator.py \
            --organism "$organism" \
            --scoreFile "$output_dir/merged.csv" \
            --outputDir "$output_dir" 2>&1 | tee -a "$output_dir/log.txt"
        ;;

    diann)
        if [ "$#" -ne 5 ]; then
            echo "Error: diann format requires 5 arguments"
            usage
        fi
        ed_file=$1
        matrix_file=$2
        output_dir=$3
        niters=$4
        imp=$5
        quant="Intensity"

        mkdir -p "$output_dir"

        echo "=== Parsing (DIA-NN) ==="
        python3 /Scripts/parse.py \
            --experimentalDesign "$ed_file" \
            --diannMatrix "$matrix_file" \
            --quantType "$quant" \
            --outputPath "$output_dir" 2>&1 | tee -a "$output_dir/log.txt"

        echo "=== Scoring ==="
        python3 /Scripts/score.py \
            --experimentalDesign "$ed_file" \
            --scoreInputs "$output_dir" \
            --outputPath "$output_dir" \
            --n-iterations "$niters" \
            --imputation "$imp" \
            --quantType "$quant" "${pi_args[@]}" 2>&1 | tee -a "$output_dir/log.txt"

        echo "=== Annotating ==="
        python3 /Scripts/annotator.py \
            --organism "$organism" \
            --scoreFile "$output_dir/merged.csv" \
            --outputDir "$output_dir" 2>&1 | tee -a "$output_dir/log.txt"
        ;;

    msstats)
        if [ "$#" -ne 5 ]; then
            echo "Error: msstats format requires 5 arguments"
            usage
        fi
        ed_file=$1
        msstats_file=$2
        output_dir=$3
        niters=$4
        imp=$5
        quant="Intensity"

        mkdir -p "$output_dir"

        if [ "$imp" != "0" ]; then
            echo "WARNING: MSstats data is already imputed by MBimpute; running AFT (imputation=$imp) will re-fit on imputed values. Pass imputation=0 to skip." | tee -a "$output_dir/log.txt"
        fi

        echo "=== Parsing (MSstats) ==="
        python3 /Scripts/parse.py \
            --experimentalDesign "$ed_file" \
            --msstatsFile "$msstats_file" \
            --quantType "$quant" \
            --outputPath "$output_dir" 2>&1 | tee -a "$output_dir/log.txt"

        echo "=== Scoring ==="
        python3 /Scripts/score.py \
            --experimentalDesign "$ed_file" \
            --scoreInputs "$output_dir" \
            --outputPath "$output_dir" \
            --n-iterations "$niters" \
            --imputation "$imp" \
            --quantType "$quant" "${pi_args[@]}" 2>&1 | tee -a "$output_dir/log.txt"

        echo "=== Annotating ==="
        python3 /Scripts/annotator.py \
            --organism "$organism" \
            --scoreFile "$output_dir/merged.csv" \
            --outputDir "$output_dir" 2>&1 | tee -a "$output_dir/log.txt"
        ;;

    saint)
        if [ "$#" -ne 7 ]; then
            echo "Error: saint format requires 7 arguments"
            usage
        fi
        bait_file=$1
        prey_file=$2
        interaction_file=$3
        quant=$4
        output_dir=$5
        niters=$6
        imp=$7

        mkdir -p "$output_dir"

        echo "=== Parsing (SAINT) ==="
        python3 /Scripts/parse.py \
            --bait "$bait_file" \
            --prey "$prey_file" \
            --interaction "$interaction_file" \
            --quantType "$quant" \
            --outputPath "$output_dir" 2>&1 | tee -a "$output_dir/log.txt"

        # For SAINT, the ED is reconstructed by parse.py in the output directory
        echo "=== Scoring ==="
        python3 /Scripts/score.py \
            --experimentalDesign "$output_dir/ED.csv" \
            --scoreInputs "$output_dir" \
            --outputPath "$output_dir" \
            --n-iterations "$niters" \
            --imputation "$imp" \
            --quantType "$quant" "${pi_args[@]}" 2>&1 | tee -a "$output_dir/log.txt"

        echo "=== Annotating ==="
        python3 /Scripts/annotator.py \
            --organism "$organism" \
            --scoreFile "$output_dir/merged.csv" \
            --outputDir "$output_dir" 2>&1 | tee -a "$output_dir/log.txt"
        ;;

    *)
        echo "Error: Unknown format '$format'. Must be maxquant, diann, fragpipe, msstats, or saint."
        usage
        ;;
esac

echo "=== Done ==="
