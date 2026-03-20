import pandas as pd
import numpy as np
import argparse

description = "This is the entry point to the program. It will execute the requested tasks."

# initialize the parser
parser = argparse.ArgumentParser(description=description)

# Add arguments for the biogrid_all file
parser.add_argument("--biogrid_all",
                    help="path to the biogrid_all file",
                    required=True)

# Add arguments for the biogrid_mv file
parser.add_argument("--biogrid_mv",
                    help="path to the biogrid_mv file",
                    required=True)

parser.add_argument("--output_dir",
                    help="path to the output directory",
                    default="/Datasets")

# NCBI Taxonomy ID for organism filtering. Common IDs:
#   Human: 9606, Mouse: 10090, Yeast (S. cerevisiae S288C): 559292
# See ORGANISMS config in setup_datasets.py and Scripts/annotator.py for the full list.
parser.add_argument("--organism_id", type=int, default=9606,
                    help="NCBI Taxonomy ID for organism filtering (default: 9606 for human)")

def preprocess_biogrid(biogrid_all_path, biogrid_mv_path, output_dir, organism_id=9606):
    # Read the biogrid_all file
    all_biogrid = pd.read_csv(biogrid_all_path, sep='\t')
    all_biogrid = all_biogrid[
        (all_biogrid['Organism ID Interactor A'] == organism_id) &
        (all_biogrid['Organism ID Interactor B'] == organism_id) &
        (all_biogrid['Experimental System Type'] == 'physical')
    ]

    # Read the biogrid_mv file
    mv_biogrid = pd.read_csv(biogrid_mv_path, sep='\t')
    mv_biogrid = mv_biogrid[
        (mv_biogrid['Organism ID Interactor A'] == organism_id) &
        (mv_biogrid['Organism ID Interactor B'] == organism_id)
    ]
    
    # Tag mv_biogrid for multivalidation
    mv_biogrid['multivalidated'] = True
    mv_biogrid = mv_biogrid.groupby(['SWISS-PROT Accessions Interactor A', 'SWISS-PROT Accessions Interactor B']).agg(
        Multivalidated=('multivalidated', 'any')
    ).reset_index()
    
    # Merge all_biogrid with mv_biogrid
    all_biogrid = pd.merge(
        all_biogrid, mv_biogrid, 
        how='left', 
        on=['SWISS-PROT Accessions Interactor A', 'SWISS-PROT Accessions Interactor B']
    )
    
    # Group by and summarize data
    summ_biogrid = all_biogrid.groupby(['SWISS-PROT Accessions Interactor A', 'SWISS-PROT Accessions Interactor B']).agg({
        'Experimental System': lambda x: '; '.join(np.array(x).astype(str)),
        'Author': lambda x: '; '.join(np.array(x).astype(str)),
        'Publication Source': lambda x: '; '.join(np.array(x).astype(str)),
        'Multivalidated': 'any'
    }).reset_index()

    summ_biogrid['In.BioGRID'] = True

    # Save the summary to CSV
    output_path = f"{output_dir}/biogrid_summary.csv"
    summ_biogrid.to_csv(output_path, index=False)

# Script entry point

args = parser.parse_args()
preprocess_biogrid(args.biogrid_all, args.biogrid_mv,
                   output_dir=args.output_dir, organism_id=args.organism_id)

