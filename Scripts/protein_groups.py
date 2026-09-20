import pandas as pd
import os
import csv
from ed_exceptions import EDNoTestExperimentsError, EDPGMismatchError
from log_config import get_logger

logger = get_logger(__name__)

# MaxQuant 2.4 and later name the decoy flag column "Decoy"; earlier versions "Reverse".
# The table is renamed to the older name on load so the rest of the pipeline sees one.
MAXQUANT_COLUMN_ALIASES = {"Decoy": "Reverse"}

# Column prefix under which proteinGroups.txt carries each quantification.
QUANT_COL_PREFIXES = {
    "Intensity": "Intensity ",
    "LFQ": "LFQ intensity ",
    "Spectral Counts": "MS/MS count ",
}

# At most this many identifiers are kept in the short protein/gene names SAINT sees.
MAX_SHORT_IDS = 3


def apply_column_aliases(df):
    """Rename alias columns to their canonical names where the canonical one is absent."""
    renames = {alias: name for alias, name in MAXQUANT_COLUMN_ALIASES.items()
               if alias in df.columns and name not in df.columns}
    return df.rename(columns=renames)


def get_quant_col_prefix(quantification):
    """The proteinGroups column prefix for a quantification type; raises on an unknown one."""
    try:
        return QUANT_COL_PREFIXES[quantification]
    except KeyError:
        raise ValueError(f"Unknown quantification type {quantification!r}; "
                         f"expected one of {sorted(QUANT_COL_PREFIXES)}")


def _shorten(ids):
    """Truncate a ';'-joined identifier list to MAX_SHORT_IDS with a dropped count, and
    replace spaces with hyphens so SAINT's whitespace-delimited files stay intact."""
    parts = ids.split(";")
    if len(parts) > MAX_SHORT_IDS:
        ids = ";".join(parts[:MAX_SHORT_IDS]) + ";+" + str(len(parts) - MAX_SHORT_IDS)
    return ids.replace(" ", "-")


class ProteinGroups:
    def __init__(self, experimental_design, file_path, quantification):
        logger.info("Parsing MaxQuant proteinGroups: %s", file_path)

        self.experimental_design = experimental_design
        self.quantification = quantification
        self.quant_col_prefix = get_quant_col_prefix(quantification)
        self.data = apply_column_aliases(
            pd.read_csv(file_path, sep="\t", quotechar="'", low_memory=False))

        design = experimental_design.name2experiment
        if not any(e.attributes["Type"] == "T" for e in design.values()):
            raise EDNoTestExperimentsError()

        # Quantification columns whose experiment the design names; the rest are ignored.
        value_type = int if quantification == "Spectral Counts" else float
        self.quant_cols = []
        self.bait_cols = []
        for col in self.data.columns:
            if not col.startswith(self.quant_col_prefix):
                continue
            exp_name = col[len(self.quant_col_prefix):]
            if exp_name not in design:
                logger.warning("Experiment found in proteinGroups but missing from experimental design: %s", exp_name)
                continue
            self.data[col] = self.data[col].astype(value_type)
            self.quant_cols.append(col)
            if design[exp_name].attributes["Type"] == "T":
                self.bait_cols.append(col)

        ed_only = sorted(set(design) - set(self._experiment_of(c) for c in self.quant_cols))
        if ed_only:
            raise EDPGMismatchError(ed_only, [])

        for column, label in (("Reverse", "reverses"),
                              ("Only identified by site", "only identified by site"),
                              ("Potential contaminant", "potential contaminants")):
            flagged = self.data[column] == "+"
            logger.info("Removed %d %s", flagged.sum(), label)
            self.data = self.data[~flagged]

        # remove proteins only found in controls
        in_a_bait = self.data[self.bait_cols].sum(axis=1) > 0
        logger.info("Removed %d proteins found only in controls", (~in_a_bait).sum())
        self.data = self.data[in_a_bait].copy()

        logger.info("Kept %d proteins", len(self.data.index))

        # SAINT-facing identifiers: gene names backfilled from the protein ID, both
        # truncated.  The original columns are kept for annotation joins.
        genes = self.data["Gene names"].where(self.data["Gene names"].notnull(),
                                              self.data["Majority protein IDs"])
        self.data["Short protein IDs"] = self.data["Majority protein IDs"].map(_shorten)
        self.data["Short Gene names"] = genes.map(_shorten)

    def _experiment_of(self, quant_col):
        return quant_col[len(self.quant_col_prefix):]

    def write_prey_file(self, out_path):
        """SAINTexpress-spc normalizes by protein length and reads it from prey.txt;
        the intensity binaries read a two-column file."""
        if self.quantification == "Spectral Counts":
            columns = ["Short protein IDs", "Sequence length", "Short Gene names"]
        else:
            columns = ["Short protein IDs", "Short Gene names"]
        self.data[columns].to_csv(out_path, index=False, sep="\t", header=False)

    def write_bait_file(self, out_path):
        with open(out_path, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            for e in self.experimental_design.name2experiment.values():
                writer.writerow([e.attributes["Experiment Name"], e.attributes["Bait"], e.attributes["Type"]])

    def write_interaction_file(self, out_path):
        design = self.experimental_design.name2experiment
        with open(out_path, 'w') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            for index, row in self.data.iterrows():
                for quant_col in self.quant_cols:
                    exp_name = self._experiment_of(quant_col)
                    writer.writerow([exp_name, design[exp_name].attributes["Bait"],
                                     row["Short protein IDs"], row[quant_col]])

    def write_CompPASS(self, out_path):
        """CompPASS input: one row per protein per experiment with a positive value.

        The header names are CompPASS's; "Experiment.ID" holds the bait name and
        "Bait" the bait's protein ID, which is what score_compPass compares to Prey.
        """
        design = self.experimental_design.name2experiment
        with open(out_path, 'w', newline='\n', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            writer.writerow(["Experiment.ID", "Replicate", "Bait", "Prey", "Prey.Name", "Spectral.Count"])
            for index, row in self.data.iterrows():
                for quant_col in self.quant_cols:
                    if row[quant_col] > 0:
                        attributes = design[self._experiment_of(quant_col)].attributes
                        writer.writerow([attributes["Bait"], attributes["Replicate"],
                                         attributes["Bait ID"], row["Short protein IDs"],
                                         row["Gene names"], row[quant_col]])

    def to_CompPASS(self, out_path):
        self.write_CompPASS(os.path.join(out_path, "to_CompPASS.csv"))

    def to_SAINT(self, out_path):
        self.write_prey_file(os.path.join(out_path, "prey.txt"))
        self.write_bait_file(os.path.join(out_path, "bait.txt"))
        self.write_interaction_file(os.path.join(out_path, "interaction.txt"))
