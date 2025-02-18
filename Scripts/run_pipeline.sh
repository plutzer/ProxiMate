#!/bin/bash

# Check if the correct number of arguments is passed
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 ED_file PG_file quant_type output_directory n_iterations imputation"
    exit 1
fi

# Assign variables to arguments
file1=$1
file2=$2
quant=$3
output_dir=$4
niters=$5
imp=$6

# Check if the output directory exists, if not create it
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

# Call the parse script
python3 /parse.py --experimentalDesign "$file1" --proteinGroups "$file2" --quantType "$quant" --outputPath "$output_dir" >> "$output_dir/log.txt"

# Call the score script
python3 /score.py --experimentalDesign "$file1" --scoreInputs "$output_dir" --outputPath "$output_dir" --n-iterations "$niters" --imputation "$imp" --quantType "$quant" >> "$output_dir/log.txt"

# Call the annotate script
python3 /annotator.py --scoreFile "$output_dir/merged.csv" --outputDir "$output_dir" >> "$output_dir/log.txt"

# Fixed