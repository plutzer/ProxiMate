"""Reads and writes the interaction file SAINTexpress consumes.

Every scoring run produces ``filtered_interaction.txt``, whatever imputation method was
selected: score.py names it on the SAINTexpress command line and _build_group_saint_inputs
filters it per group.  SAINTexpress parses it positionally, so the four columns and the
absent header are part of the contract, not formatting.
"""

import pandas as pd

SAINT_COLUMNS = ['ExperimentID', 'Bait', 'Prey', 'Intensity']


def read_saint_inputs(interaction_path, ed_path):
    """Read the SAINT interaction file and the experimental design.

    Returns (interaction, ed, bait_dict): the interaction frame with SAINT_COLUMNS as
    header, the design frame, and a bait name -> bait protein ID map the imputation
    paths use to drop a prey's own-bait rows before fitting.
    """
    interaction = pd.read_csv(interaction_path, sep='	', header=None)
    interaction.columns = SAINT_COLUMNS
    ed = pd.read_csv(ed_path)
    bait_dict = dict(zip(ed['Bait'], ed['Bait ID']))
    return interaction, ed, bait_dict


def write_filtered_interaction(interaction, output_dir):
    """Drop non-positive intensities and write the 4-column SAINT interaction file.

    `interaction` may carry extra working columns — the imputation paths attach a BaitID
    helper — so the SAINT columns are selected explicitly rather than written wholesale.
    `output_dir` is concatenated directly, and must therefore end in a separator.

    Returns (kept, total) row counts so callers can report what the filter removed.
    """
    positive = interaction[interaction['Intensity'] > 0]
    positive[SAINT_COLUMNS].to_csv(
        output_dir + 'filtered_interaction.txt', sep='\t', index=False, header=False)
    return len(positive), len(interaction)
