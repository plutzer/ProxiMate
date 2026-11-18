import pandas as pd
import numpy as np
import argparse

A = 'SWISS-PROT Accessions Interactor A'
B = 'SWISS-PROT Accessions Interactor B'


def mirror(df):
    """Add each A-B row again as B-A.

    BioGRID records a pair once, bait in A.  Listing both orientations lets a lookup
    find the pair whichever partner was the bait.  Self-pairs are not repeated.
    """
    swapped = df[df[A] != df[B]].rename(columns={A: B, B: A})
    return pd.concat([df, swapped], ignore_index=True)


def preprocess_biogrid(biogrid_all_path, biogrid_mv_path, output_dir, organism_id=9606,
                       exclude_publication=None, output_filename="biogrid_summary.csv"):
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

    # Excluding a publication removes its evidence rows from both files; a pair with no
    # rows left drops out of the summary below, one with other evidence keeps it.
    if exclude_publication is not None:
        all_biogrid = all_biogrid[all_biogrid['Publication Source'] != exclude_publication]
        mv_biogrid = mv_biogrid[mv_biogrid['Publication Source'] != exclude_publication]

    all_biogrid = mirror(all_biogrid[[A, B, 'Experimental System', 'Author', 'Publication Source']])
    mv_biogrid = mirror(mv_biogrid[[A, B, 'Experimental System', 'Publication Source']])

    # Every pair in the multivalidated file has two publications or two experimental
    # systems.  With a publication's rows excluded, a pair has to meet that again on the
    # evidence that is left.
    mv_biogrid = mv_biogrid.groupby([A, B])[['Experimental System', 'Publication Source']].nunique()
    if exclude_publication is not None:
        mv_biogrid = mv_biogrid[mv_biogrid.max(axis=1) >= 2]
    mv_biogrid = mv_biogrid.assign(Multivalidated=True)[['Multivalidated']].reset_index()

    # Merge all_biogrid with mv_biogrid
    all_biogrid = pd.merge(all_biogrid, mv_biogrid, how='left', on=[A, B])

    # Group by and summarize data
    summ_biogrid = all_biogrid.groupby([A, B]).agg({
        'Experimental System': lambda x: '; '.join(np.array(x).astype(str)),
        'Author': lambda x: '; '.join(np.array(x).astype(str)),
        'Publication Source': lambda x: '; '.join(np.array(x).astype(str)),
        'Multivalidated': 'any'
    }).reset_index()

    summ_biogrid['In.BioGRID'] = True

    # Save the summary to CSV
    output_path = f"{output_dir}/{output_filename}"
    summ_biogrid.to_csv(output_path, index=False)

def main():
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

    parser.add_argument("--output_filename",
                        help="name of the summary file written into output_dir",
                        default="biogrid_summary.csv")

    # NCBI Taxonomy ID for organism filtering. Common IDs:
    #   Human: 9606, Mouse: 10090, Yeast (S. cerevisiae S288C): 559292
    # See ORGANISMS in setup_datasets.py for the full list.
    parser.add_argument("--organism_id", type=int, default=9606,
                        help="NCBI Taxonomy ID for organism filtering (default: 9606 for human)")

    parser.add_argument("--exclude_publication", default=None,
                        help="drop every evidence row with this 'Publication Source' value "
                             "(e.g. PUBMED:34079125) before summarizing; interactions that keep "
                             "other evidence remain, multivalidated only if that evidence has "
                             "two publications or two experimental systems")

    args = parser.parse_args()
    preprocess_biogrid(args.biogrid_all, args.biogrid_mv,
                       output_dir=args.output_dir, organism_id=args.organism_id,
                       exclude_publication=args.exclude_publication,
                       output_filename=args.output_filename)


if __name__ == "__main__":
    main()
