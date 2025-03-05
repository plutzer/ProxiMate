import numpy as np
import pandas as pd
import argparse
import re
import csv
import time
import subprocess

def get_first_SCL(item):
    if pd.isnull(item):
        return np.nan
    else:
        # Check if there are commas in the item
        annotation = item.split('SUBCELLULAR LOCATION: ')[1].split(';')[0].split(' {')[0]
        if ']: ' in annotation:
            annotation = annotation.split(']: ')[1]
        if ',' in annotation:
            annotation = annotation.split(',')[0]
        if '.' in annotation:
            annotation = annotation.split('.')[0]
        return annotation

# GO CC Annotations:
def clean_gocc(s):
    # Use a regular expression to remove text within brackets
    cleaned = re.sub(r'\s*\[.*?\]\s*', '', s)
    # Clean up any resulting multiple semicolons and trim spaces around semicolons
    cleaned = re.sub(r'\s*;\s*', '; ', cleaned).strip()
    # If the cleaned string ends with a semicolon, remove it
    if cleaned.endswith(';'):
        cleaned = cleaned[:-1]
    return cleaned

def trim_GO_CC(item):
    if pd.isnull(item):
        return np.nan
    else:
        return clean_gocc(item)
    
# Motif Annotations:
def trim_motifs(item):
    if pd.isnull(item):
        return np.nan
    else:
        return clean_motif(item)

def clean_motif(s):
    # Regular expression to find /note="..." patterns
    motif_names = re.findall(r'/note="([^"]*)"', s)
    # Join the extracted motif names with a semicolon and space
    return '; '.join(motif_names)

# TODO: Add an annotation for if it is a self-interaction
def self_inter(prey_id, bait_id):
    prey_id=prey_id.split(';')
    for prey in prey_id:
            if prey == bait_id:
                return True
    return False

# TODO: Add an annotation for if the prey is a bait
def prey_is_bait(prey_id, bait_values):
    prey_id=prey_id.split(';')
    for prey in prey_id:
        if prey in bait_values:
            return True
    return False

# TODO: Add an annotation for the main location
#Get first prey-gene -> new column
def get_first_pg(item):
    return item.split(';')[0]

#make sure each First_Prey_Gene name is in subcellular
#if name cannot be found, search uniprot['Gene Names'] for it
#take list in cell that includes First_Prey_Gene name and split it
#search subcellular for each name until one hits, and return that name
#works, total runtime after adding is 3 minutes in docker
def get_match(item, subcellular, uniprot):
    if item in subcellular:
        return item

    for cells in uniprot:
        # print(cells)
        names = str(cells).split(' ')
        if item in names:
            for name in names:
                if name in subcellular:
                    return name

# Function that can add GOGO dictionary scores back to the proximity data
def get_cco_score(bait_gene,prey_gene,cc_dict):
    if str(bait_gene) in cc_dict:
        if ';' in str(prey_gene):
            for item in prey_gene.split(';'):
                if item in cc_dict[str(bait_gene)]:
                    return cc_dict[str(bait_gene)][item]
        else:
            if str(prey_gene) in cc_dict[str(bait_gene)]:
                return cc_dict[str(bait_gene)][str(prey_gene)]
    return np.nan # return nan if the gene pair is not in the dictionary

#make sure each First_Prey_Gene name is in subcellular
#if name cannot be found, search uniprot['Gene Names'] for it
#take list in cell that includes First_Prey_Gene name and split it
#search subcellular for each name until one hits, and return that name
#works, total runtime after adding is 3 minutes in docker
def get_match(item, subcellular, uniprot):
    if item in subcellular:
        return item
 
    for cells in uniprot:
        names = str(cells).split(' ')
        if item in names:
            for name in names:
                if name in subcellular:
                    return name

def complex_id(prey_id, complex_dict):
    prey_id = prey_id.split(';')
    for prey in prey_id:
        for key,value in complex_dict.items():
            ids = value.split(';')
            if prey in ids:
                return key
    return None

def main():
    description = "This is the entry point to the program. It will execute the requested tasks."

    # initialize the parser
    parser = argparse.ArgumentParser(description=description)

    # Add arguments for the annotator
    parser.add_argument("--scoreFile",
                        help="path to raw scores file",
                        default="/srv/shiny-server/myapp/score_outputs/merged.csv")

    parser.add_argument("--uniprotFile",
                        help="path to uniprot annotation file",
                        default="/Datasets/uniprot_anns.tsv")

    parser.add_argument("--biogridFile",
                        help="path to biogrid annotation file",
                        default="/Datasets/biogrid_summary.csv")

    parser.add_argument("--locationFile", 
                        help="path to subcellular locations file", 
                        default="/Datasets/subcellular_location.tsv")

    parser.add_argument("--complexFile",
                        help="path to human complexes file",
                        default="/Datasets/humanComplexes.txt")

    # Add arguments for the prey and bait columns
    parser.add_argument("--preyColumn",
                        help="name of the prey column",
                        default="Prey.ID")

    parser.add_argument("--baitColumn",
                        help="name of the bait column",
                        default="Bait.ID")

    # Add argument for PreyGene column
    parser.add_argument("--preyGeneColumn",
                        help="name of the prey gene column",
                        default="PreyGene")

    parser.add_argument("--outputDir",
                        help="path to the output directory",
                        default="/srv/shiny-server/myapp/output")

    args = parser.parse_args()

    prey_col = args.preyColumn
    bait_col = args.baitColumn
    gene_col = args.preyGeneColumn

    # Import the uniprot annotations
    uniprot = pd.read_csv(args.uniprotFile, sep='\t')

    # SCL Annotations:
    # Take the first item from the SCL column and move to a new column
    uniprot['first_SCL'] = uniprot['Subcellular location [CC]'].apply(get_first_SCL)
    # For any item in the first_SCL column that only occurs once, change it to NaN - this is a bit lazy but will only affect ~16 proteins with shitty annotations
    uniprot['first_SCL'] = uniprot['first_SCL'].where(uniprot['first_SCL'].map(uniprot['first_SCL'].value_counts()) > 1, np.nan)

    # Creating cleaner versions of other useful annotations:
    # Now some of the columns to produce lists:
    uniprot['GO_CC'] = uniprot['Gene Ontology (cellular component)'].apply(trim_GO_CC)
    uniprot['Motifs'] = uniprot['Motif'].apply(trim_motifs)
    uniprot['Regions'] = uniprot['Region'].apply(trim_motifs)
    uniprot['Repeats'] = uniprot['Repeat'].apply(trim_motifs)
    uniprot['Compositions'] = uniprot['Compositional bias'].apply(trim_motifs)
    uniprot['Domains'] = uniprot['Domain [FT]'].apply(trim_motifs)

    # Now Merge the uniprot annotations with the raw scores
    raw_scores = pd.read_csv(args.scoreFile)

    raw_scores['First_ID'] = raw_scores[prey_col].apply(lambda x: x.split(';')[0])
    first_prey_col = 'First_ID'

    # Merge the raw scores with the uniprot annotations
    annotated_scores = raw_scores.merge(uniprot, left_on=first_prey_col, right_on='Entry', how='left')

    annotated_scores['First_Prey_Gene'] = annotated_scores[gene_col].apply(get_first_pg)

    #subset human protein atlas to just two columns
    hpa = pd.read_csv(args.locationFile, sep='\t')
    name_loc = pd.DataFrame(hpa.iloc[:, [1,3]]) # TODO: Make this not hardcoded

    annotated_scores['Matched_Gene_Name'] = annotated_scores['First_Prey_Gene'].apply(get_match, subcellular=name_loc['Gene name'].to_numpy(), uniprot=uniprot['Gene Names'].to_numpy())

    annotated_scores = annotated_scores.merge(name_loc, left_on=['Matched_Gene_Name'], right_on=['Gene name'], how='left')

    bait_values = set(annotated_scores[bait_col])
    annotated_scores['Prey_Is_Bait'] = annotated_scores[prey_col].apply(prey_is_bait, bait_values=bait_values)

    for bait_id in bait_values:
        annotated_scores['Self-Interaction'] = annotated_scores[prey_col].apply(self_inter, bait_id=bait_id)

    # TODO: add human complex annotations
    # subunits(UniProt IDs) match with Prey ID(?), split over semicolon
    # if match, add annotation for ComplexName
    human_complex = pd.read_table(parser.parse_args().complexFile, encoding='latin-1')
    complex_cols = human_complex[['complex_name','subunits_uniprot_id']]
    complex_dict = complex_cols.set_index('complex_name').to_dict()['subunits_uniprot_id']

    annotated_scores['Human_Complex'] = annotated_scores[prey_col].apply(complex_id, complex_dict=complex_dict)

    # Now for BioGrid
    # Import the BioGrid annotations
    biogrid = pd.read_csv(parser.parse_args().biogridFile)
    biogrid = pd.read_csv(args.biogridFile)

    # Convert the integer column to string in both dataframes before merging
    annotated_scores[first_prey_col] = annotated_scores[first_prey_col].astype(str)
    annotated_scores[bait_col] = annotated_scores[bait_col].astype(str)
    biogrid['SWISS-PROT Accessions Interactor A'] = biogrid['SWISS-PROT Accessions Interactor A'].astype(str)
    biogrid['SWISS-PROT Accessions Interactor B'] = biogrid['SWISS-PROT Accessions Interactor B'].astype(str)
    # Merge the annotated scores with the BioGrid annotations
    annotated_scores = annotated_scores.merge(biogrid, left_on=[first_prey_col, bait_col], right_on=['SWISS-PROT Accessions Interactor A', 'SWISS-PROT Accessions Interactor B'], how='left')

    # Fill in the In.BioGRID column with False for rows that have nan
    annotated_scores['In.BioGRID'] = annotated_scores['In.BioGRID'].fillna(False)
    # Do the same with Multivalidated
    annotated_scores['Multivalidated'] = annotated_scores['Multivalidated'].fillna(False)

    # First need to write an input file for GOGO
    # Input file should be written to the same directory that the merged scores file came from. So from the scoreFile path remove everything after the last / and add gogo_input.txt
    gogo_input_path = args.scoreFile.rsplit('/', 1)[0] + '/gogo_input.txt'
    print("GOGO input: " + gogo_input_path)

    start_time = time.time()
    bait_anns = {}
    with open(gogo_input_path, "w") as f:
        for index, row in annotated_scores.iterrows():
            bait_id = row[bait_col]
            prey_id = row[first_prey_col]
            go_anns_raw = row['Gene Ontology (cellular component)'] # Go anns are in this format: cytosolic small ribosomal subunit [GO:0022627]; nucleus [GO:0005634]; plasma membrane [GO:0005886]
            if pd.isnull(go_anns_raw):
                # print(f"GO anns are null for {prey_id}")
                # Skip this row
                continue
            else:
                go_anns_split = go_anns_raw.split('; ')
                # Now for each element get just the GO ID
                go_ids = [ann.split(' [')[1].split(']')[0] for ann in go_anns_split]

                # Now need the GO anns for the bait - use the uniprot dataframe for this
                if bait_id in bait_anns:
                    if bait_anns[bait_id] == []:
                        # print(f"GO anns are empty for {bait_id}")
                        continue
                    else:
                        bait_go_ids = bait_anns[bait_id]
                else:
                    bait_go_anns = uniprot[uniprot['Entry'] == bait_id]['Gene Ontology (cellular component)'].values
                    if bait_go_anns.size == 0:
                        bait_anns[bait_id] = []
                        # print(f"GO anns are null for {bait_id}")
                        continue
                    elif pd.isnull(bait_go_anns):
                        bait_anns[bait_id] = []
                        # print(f"GO anns are null for {bait_id}")
                        continue
                    else:
                        bait_go_anns_split = bait_go_anns[0].split('; ')
                        bait_go_ids = [ann.split(' [')[1].split(']')[0] for ann in bait_go_anns_split]
                        bait_anns[bait_id] = bait_go_ids
            # Now append the line to the file
            f.write(f"{bait_id} {' '.join(bait_go_ids)};{prey_id} {' '.join(go_ids)}\n")

    print(f"Time taken to write GOGO input file: {time.time() - start_time}")

    # TODO: GOGO subprocess call goes here
    p = subprocess.run(["perl",
                        "gene_pair_comb.pl",
                        str(gogo_input_path),
                        str(args.outputDir) + "/gogo_output.txt"],
                        cwd="/GOGO")

    # Now read in the GOGO output file
    output_filename = str(args.outputDir) + "/gogo_output.txt" # This will be the output file from GOGO

    # Initialize the dictionary
    cc_dict = {}

    counter = 0

    # Read in the output file line by line
    with open(output_filename) as f:
        content = f.readlines()
        for item in content:
            line = item.strip()

            baitgene = line.split(' ') [0]
            preygene = line.split(';')[1].split(' ')[0]

            # Now I need to get the CCO score
            # The score comes after 'CCO'
            score = float(line.split('CCO')[1].strip().split(' ')[0])

            # Now I need to add this to the dictionary
            if baitgene in cc_dict:
                cc_dict[baitgene][preygene] = score
            else:
                cc_dict[baitgene] = {preygene: score}


    # Now I can add the CCO scores to the annotated scores
    annotated_scores['CCO'] = annotated_scores.apply(lambda x: get_cco_score(x[bait_col], x[first_prey_col], cc_dict), axis=1)

    # TODO: This is where we will filter the annotations based on user preference.
    # TODO: Also need to re-order the columns in a way that makes the most sense.

    # Save the annotated scores
    output_path = f"{args.outputDir}/annotated_scores.csv"
    annotated_scores.to_csv(output_path, index=False)

if __name__ == "__main__":
    main()
