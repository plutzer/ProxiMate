import numpy as np
import pandas as pd
import argparse
import re
import csv
import os
import sys
import shutil
import time
import subprocess
import provenance
from log_config import get_logger, add_file_handler

logger = get_logger(__name__)

# Supported organisms for ProxiMate annotation.
# To add a new organism:
#   1. Add an entry here with its NCBI Taxonomy ID and feature flags
#   2. Sync this config in both setup_datasets.py and Scripts/annotator.py
#   3. Add the organism to the GUI dropdown in GUI/app.py (scoring panel)
#   4. If the organism has a species-specific database (like HPA for human),
#      add a download function in setup_datasets.py and conditional logic below
ORGANISMS = {
    "human": {"organism_id": 9606, "has_hpa": True, "has_corum": True, "has_hcm": True},
    "mouse": {"organism_id": 10090, "has_hpa": False, "has_corum": False, "has_hcm": False},
    "yeast": {"organism_id": 559292, "has_hpa": False, "has_corum": False, "has_hcm": False},
}

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

def complex_id(prey_id, complex_dict):
    prey_id = prey_id.split(';')
    for prey in prey_id:
        for key,value in complex_dict.items():
            ids = value.split(';')
            if prey in ids:
                return key
    return None

# A UniProt accession, with an optional isoform suffix.
ACCESSION_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?$")


def symbol_accession_map(uniprot):
    """Gene symbol -> accession, from the space-separated 'Gene Names' column.

    A symbol listed under more than one entry is left out: picking one would annotate
    the wrong protein without any sign of it.
    """
    seen = {}
    for accession, names in zip(uniprot['Entry'], uniprot['Gene Names']):
        if pd.isnull(names):
            continue
        for symbol in str(names).split():
            seen.setdefault(symbol, set()).add(accession)
    return {symbol: next(iter(accs)) for symbol, accs in seen.items() if len(accs) == 1}


def resolve_accessions(ids, symbol_map):
    """Resolve one ';'-joined identifier string to accessions.

    Accessions pass through, a known symbol becomes its accession, and anything else is
    kept as written so the caller can count what did not resolve.
    """
    resolved = []
    for item in str(ids).split(';'):
        item = item.strip()
        if ACCESSION_RE.match(item):
            resolved.append(item)
        else:
            resolved.append(symbol_map.get(item, item))
    return ';'.join(resolved)


def unresolved_ids(resolved_series):
    """The distinct identifiers in a resolved column that are still not accessions."""
    return sorted({item for ids in resolved_series for item in str(ids).split(';')
                   if not ACCESSION_RE.match(item)})


def check_annotation_coverage(matched, total, source):
    """Report an annotation source that matched nothing at all.

    Across a whole dataset, zero matches is far more often a mismatch between
    the identifiers being joined — the wrong column, or placeholder IDs from an
    input format that has none — than a real absence of annotation.  Either way
    the output is a column of False, which reads as a negative result.
    """
    if total and not matched:
        logger.warning(
            "%s matched none of the %d interactions; check that the identifiers "
            "being joined are of the same kind.", source, total)
    return matched


def collapse_hpa_locations(name_loc):
    """Reduce the HPA gene/location table to one row per gene.

    ``subcellular_location.tsv`` repeats some gene names.  Merged as-is on gene
    name, each repeat multiplies every interaction row whose prey is that gene,
    inflating network sizes and enrichment counts with rows that look real.

    Most repeats carry an identical location and collapse cleanly.  A few carry
    genuinely different ones, which are joined rather than resolved by row order
    — dropping one would discard a real annotation — and reported, so a future
    release that introduces a new conflict is visible rather than silent.
    """
    deduped = name_loc.drop_duplicates()

    conflicting = deduped[deduped.duplicated(subset=['Gene name'], keep=False)]
    if not conflicting.empty:
        genes = sorted(conflicting['Gene name'].dropna().unique())
        logger.warning(
            "HPA lists differing main locations for %d gene(s); joining them: %s",
            len(genes), ", ".join(genes))

    def join_locations(values):
        distinct = sorted(values.dropna().unique())
        return "; ".join(distinct) if distinct else None

    collapsed = (deduped.groupby('Gene name', as_index=False, sort=False)
                        ['Main location'].agg(join_locations))
    logger.info("HPA subcellular locations: %d genes from %d rows",
                len(collapsed), len(name_loc))
    return collapsed


def main():
    description = "This is the entry point to the program. It will execute the requested tasks."

    # initialize the parser
    parser = argparse.ArgumentParser(description=description)

    # Add arguments for the annotator
    parser.add_argument("--organism",
                        choices=list(ORGANISMS.keys()), default="human",
                        help="Organism for annotation databases (default: human)")

    parser.add_argument("--scoreFile",
                        help="path to raw scores file",
                        default="/srv/shiny-server/myapp/score_outputs/merged.csv")

    parser.add_argument("--uniprotFile", default=None,
                        help="path to uniprot annotation file (default: /Datasets/{organism}/uniprot_anns.tsv)")

    parser.add_argument("--biogridFile", default=None,
                        help="path to biogrid annotation file (default: /Datasets/{organism}/biogrid_summary.csv)")

    parser.add_argument("--excludeHCM", action="store_true",
                        help="annotate against the BioGRID summary with Human Cell Map "
                             "(Go et al. 2021) evidence removed; human only")

    parser.add_argument("--locationFile", default=None,
                        help="path to subcellular locations file (HPA, human only)")

    parser.add_argument("--complexFile", default=None,
                        help="path to protein complexes file (CORUM, human only)")

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

    # Attach the dataset log before resolving anything, so a failure while
    # resolving organism defaults is recorded in the output directory too.
    os.makedirs(args.outputDir, exist_ok=True)
    add_file_handler(os.path.join(args.outputDir, "proximate.log"))

    with provenance.stage(args.outputDir, "annotate", entrypoint="annotator.main",
                          cli_args=vars(args)) as record:
        _annotate(args, record)


def _annotate(args, record):
    """Annotate the scored interactions, recording provenance into `record`."""
    # Resolve organism-specific default paths
    datasets_dir = "/Datasets"
    organism_dir = f"{datasets_dir}/{args.organism}"
    org_config = ORGANISMS[args.organism]

    if args.excludeHCM and not org_config["has_hcm"]:
        raise ValueError(f"--excludeHCM: no Human Cell Map variant exists for {args.organism}")

    if args.uniprotFile is None:
        args.uniprotFile = f"{organism_dir}/uniprot_anns.tsv"
    if args.biogridFile is None:
        args.biogridFile = provenance.biogrid_summary_path(
            args.organism, datasets_dir, exclude_hcm=args.excludeHCM)
    if args.locationFile is None and org_config["has_hpa"]:
        args.locationFile = f"{organism_dir}/subcellular_location.tsv"
    if args.complexFile is None and org_config["has_corum"]:
        args.complexFile = f"{datasets_dir}/corum_humanComplexes.txt"

    prey_col = args.preyColumn
    bait_col = args.baitColumn
    gene_col = args.preyGeneColumn

    record.extra(organism=args.organism, exclude_hcm=args.excludeHCM, prey_column=prey_col,
                 bait_column=bait_col, gene_column=gene_col)
    for role in ("scoreFile", "uniprotFile", "biogridFile", "locationFile", "complexFile"):
        path = getattr(args, role)
        if path:
            record.add_input(path, role=role)

    # Validate that the score file exists
    if not os.path.exists(args.scoreFile):
        logger.error("Score file not found: %s", args.scoreFile)
        sys.exit(1)

    # Import the uniprot annotations
    logger.info("Loading UniProt annotations from %s", args.uniprotFile)
    try:
        uniprot = pd.read_csv(args.uniprotFile, sep='\t')
    except Exception:
        logger.exception("Failed to load UniProt annotations")
        sys.exit(1)

    # SCL Annotations:
    # Take the first item from the SCL column and move to a new column
    uniprot['first_SCL'] = uniprot['Subcellular location [CC]'].apply(get_first_SCL)
    # For any item in the first_SCL column that only occurs once, change it to NaN - this is a bit lazy but will only affect ~16 proteins with unreliable annotations
    uniprot['first_SCL'] = uniprot['first_SCL'].where(uniprot['first_SCL'].map(uniprot['first_SCL'].value_counts()) > 1, np.nan)

    # Creating cleaner versions of other useful annotations:
    # Now some of the columns to produce lists:
    uniprot['GO_CC'] = uniprot['Gene Ontology (cellular component)'].apply(trim_GO_CC)
    uniprot['GO_BP'] = uniprot['Gene Ontology (biological process)'].apply(trim_GO_CC)
    uniprot['GO_MF'] = uniprot['Gene Ontology (molecular function)'].apply(trim_GO_CC)
    uniprot['Motifs'] = uniprot['Motif'].apply(trim_motifs)
    uniprot['Regions'] = uniprot['Region'].apply(trim_motifs)
    uniprot['Repeats'] = uniprot['Repeat'].apply(trim_motifs)
    uniprot['Compositions'] = uniprot['Compositional bias'].apply(trim_motifs)
    uniprot['Domains'] = uniprot['Domain [FT]'].apply(trim_motifs)

    # Now Merge the uniprot annotations with the raw scores
    logger.info("Loading scored data from %s", args.scoreFile)
    try:
        raw_scores = pd.read_csv(args.scoreFile)
        logger.info("Scored data: %d rows, columns: %s", len(raw_scores), list(raw_scores.columns))
        record.metric("scored_rows_in", len(raw_scores))
    except Exception:
        logger.exception("Failed to load score file")
        sys.exit(1)

    # Every lookup below is keyed on accessions.  Inputs that carry gene symbols instead
    # are resolved here, once; the supplied Prey.ID and Bait.ID columns stay as given
    # because the GUI keys on them.
    symbol_map = symbol_accession_map(uniprot)
    raw_scores['Prey_Accessions'] = raw_scores[prey_col].apply(
        resolve_accessions, symbol_map=symbol_map)
    raw_scores['First_ID'] = raw_scores['Prey_Accessions'].str.split(';').str[0]
    # SAINT bait files carry no protein ID; the bait name is then tried as a symbol.
    bait_source = raw_scores[bait_col].where(raw_scores[bait_col].notna(), raw_scores['Experiment.ID'])
    raw_scores['Bait_Accession'] = bait_source.astype(str).apply(
        resolve_accessions, symbol_map=symbol_map)
    first_prey_col = 'First_ID'

    n_symbols = int((raw_scores['First_ID'] != raw_scores[prey_col].str.split(';').str[0]).sum())
    record.metric("prey_symbols_resolved", n_symbols)
    if n_symbols:
        logger.info("Resolved gene symbols to accessions for %d prey rows", n_symbols)
    for label, column in (("prey", 'First_ID'), ("bait", 'Bait_Accession')):
        unresolved = unresolved_ids(raw_scores[column])
        record.metric(f"{label}_ids_unresolved", len(unresolved))
        if unresolved:
            logger.warning(
                "%d distinct %s identifier(s) are neither UniProt accessions nor known "
                "gene symbols and will not be annotated, e.g. %s",
                len(unresolved), label, ", ".join(unresolved[:5]))

    # Merge the raw scores with the uniprot annotations
    annotated_scores = raw_scores.merge(uniprot, left_on=first_prey_col, right_on='Entry', how='left')

    annotated_scores['First_Prey_Gene'] = annotated_scores[gene_col].apply(get_first_pg)

    # HPA subcellular location annotations (human only)
    if args.locationFile:
        logger.info("Loading HPA annotations from %s", args.locationFile)
        try:
            hpa = pd.read_csv(args.locationFile, sep='\t')
            name_loc = collapse_hpa_locations(hpa[['Gene name', 'Main location']])

            annotated_scores['Matched_Gene_Name'] = annotated_scores['First_Prey_Gene'].apply(get_match, subcellular=name_loc['Gene name'].to_numpy(), uniprot=uniprot['Gene Names'].to_numpy())

            annotated_scores = annotated_scores.merge(name_loc, left_on=['Matched_Gene_Name'], right_on=['Gene name'], how='left')
        except Exception:
            logger.exception("HPA annotation failed (non-fatal, continuing)")

    bait_values = set(annotated_scores[bait_col])
    annotated_scores['Prey_Is_Bait'] = annotated_scores[prey_col].apply(prey_is_bait, bait_values=bait_values)

    for bait_id in bait_values:
        annotated_scores['Self-Interaction'] = annotated_scores[prey_col].apply(self_inter, bait_id=bait_id)

    # CORUM protein complex annotations (human only)
    if args.complexFile:
        logger.info("Loading CORUM annotations from %s", args.complexFile)
        try:
            human_complex = pd.read_table(args.complexFile, encoding='latin-1')
            complex_cols = human_complex[['complex_name','subunits_uniprot_id']]
            complex_dict = complex_cols.set_index('complex_name').to_dict()['subunits_uniprot_id']
            annotated_scores['Human_Complex'] = annotated_scores['Prey_Accessions'].apply(complex_id, complex_dict=complex_dict)
        except Exception:
            logger.exception("CORUM annotation failed (non-fatal, continuing)")

    # Now for BioGrid
    logger.info("Loading BioGRID annotations from %s", args.biogridFile)
    try:
        biogrid = pd.read_csv(args.biogridFile)

        # Convert the integer column to string in both dataframes before merging
        annotated_scores[first_prey_col] = annotated_scores[first_prey_col].astype(str)
        biogrid['SWISS-PROT Accessions Interactor A'] = biogrid['SWISS-PROT Accessions Interactor A'].astype(str)
        biogrid['SWISS-PROT Accessions Interactor B'] = biogrid['SWISS-PROT Accessions Interactor B'].astype(str)
        # Merge the annotated scores with the BioGrid annotations
        annotated_scores = annotated_scores.merge(biogrid, left_on=[first_prey_col, 'Bait_Accession'], right_on=['SWISS-PROT Accessions Interactor A', 'SWISS-PROT Accessions Interactor B'], how='left')

        # Fill in the In.BioGRID column with False for rows that have nan
        annotated_scores['In.BioGRID'] = annotated_scores['In.BioGRID'].fillna(False)
        # Do the same with Multivalidated
        annotated_scores['Multivalidated'] = annotated_scores['Multivalidated'].fillna(False)
    except Exception:
        logger.exception("BioGRID annotation failed (non-fatal, continuing)")

    # First need to write an input file for GOGO
    # Input file should be written to the same directory that the merged scores file came from. So from the scoreFile path remove everything after the last / and add gogo_input.txt
    gogo_input_path = args.scoreFile.rsplit('/', 1)[0] + '/gogo_input.txt'
    logger.info("Writing GOGO input to %s", gogo_input_path)

    start_time = time.time()
    bait_anns = {}
    with open(gogo_input_path, "w") as f:
        for index, row in annotated_scores.iterrows():
            bait_id = row['Bait_Accession']
            prey_id = row[first_prey_col]
            go_anns_raw = row['Gene Ontology (cellular component)'] # Go anns are in this format: cytosolic small ribosomal subunit [GO:0022627]; nucleus [GO:0005634]; plasma membrane [GO:0005886]
            if pd.isnull(go_anns_raw):
                continue
            else:
                go_anns_split = go_anns_raw.split('; ')
                # Now for each element get just the GO ID
                go_ids = [ann.split(' [')[1].split(']')[0] for ann in go_anns_split]

                # Now need the GO anns for the bait - use the uniprot dataframe for this
                if bait_id in bait_anns:
                    if bait_anns[bait_id] == []:
                        continue
                    else:
                        bait_go_ids = bait_anns[bait_id]
                else:
                    bait_go_anns = uniprot[uniprot['Entry'] == bait_id]['Gene Ontology (cellular component)'].values
                    if bait_go_anns.size == 0:
                        bait_anns[bait_id] = []
                        continue
                    elif pd.isnull(bait_go_anns):
                        bait_anns[bait_id] = []
                        continue
                    else:
                        bait_go_anns_split = bait_go_anns[0].split('; ')
                        bait_go_ids = [ann.split(' [')[1].split(']')[0] for ann in bait_go_anns_split]
                        bait_anns[bait_id] = bait_go_ids
            # Now append the line to the file
            f.write(f"{bait_id} {' '.join(bait_go_ids)};{prey_id} {' '.join(go_ids)}\n")

    logger.info("GOGO input written in %.1f seconds", time.time() - start_time)

    # Run GOGO subprocess
    gogo_output_path = str(args.outputDir) + "/gogo_output.txt"
    logger.info("Running GOGO (gene_pair_comb.pl)...")
    p = subprocess.run(["perl",
                        "/Scripts/GOGO/gene_pair_comb.pl",
                        str(gogo_input_path),
                        gogo_output_path],
                        cwd="/Scripts/GOGO",
                        capture_output=True, text=True)

    if p.returncode != 0:
        logger.error("GOGO subprocess failed (exit code %d)", p.returncode)
        if p.stderr:
            logger.error("GOGO stderr:\n%s", p.stderr)
        sys.exit(1)
    else:
        logger.info("GOGO completed successfully")

    # Now read in the GOGO output file
    logger.info("Reading GOGO output from %s", gogo_output_path)
    if not os.path.exists(gogo_output_path):
        logger.error("GOGO output file not found: %s", gogo_output_path)
        sys.exit(1)

    # Initialize the dictionary
    cc_dict = {}

    try:
        with open(gogo_output_path) as f:
            content = f.readlines()
            for line_num, item in enumerate(content, 1):
                line = item.strip()
                if not line:
                    continue

                baitgene = line.split(' ') [0]
                preygene = line.split(';')[1].split(' ')[0]

                # Now I need to get the CCO score
                # The score comes after 'CCO'
                score_str = line.split('CCO')[1].strip().split(' ')[0]
                if score_str == 'NA':
                    continue
                score = float(score_str)

                # Now I need to add this to the dictionary
                if baitgene in cc_dict:
                    cc_dict[baitgene][preygene] = score
                else:
                    cc_dict[baitgene] = {preygene: score}
    except Exception:
        logger.exception("Failed to parse GOGO output")
        sys.exit(1)

    # Now I can add the CCO scores to the annotated scores
    annotated_scores['CCO'] = annotated_scores.apply(lambda x: get_cco_score(x['Bait_Accession'], x[first_prey_col], cc_dict), axis=1)

    # Save the annotated scores
    output_path = f"{args.outputDir}/annotated_scores.csv"
    annotated_scores.to_csv(output_path, index=False)
    logger.info("Annotated scores written to %s (%d rows)", output_path, len(annotated_scores))

    record.metric("annotated_rows", len(annotated_scores))
    if "In.BioGRID" in annotated_scores.columns:
        in_biogrid = int(annotated_scores["In.BioGRID"].sum())
        check_annotation_coverage(in_biogrid, len(annotated_scores), "BioGRID")
        record.metric("in_biogrid", in_biogrid)
    record.add_output(output_path, rows=len(annotated_scores))

    # Copy build info to output directory so users know dataset versions
    build_info = f"{datasets_dir}/build_info.txt"
    if os.path.isfile(build_info):
        shutil.copy2(build_info, f"{args.outputDir}/build_info.txt")

if __name__ == "__main__":
    main()
