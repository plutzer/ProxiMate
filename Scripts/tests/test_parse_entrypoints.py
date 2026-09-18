"""End-to-end tests for the MaxQuant, DIA-NN, Pioneer and FragPipe parse entry points.

Each writes the five files the scoring stage reads, and each reports the experiment counts
that reach run.json and the GUI.  The MSstats entry point is covered by
test_parse_msstats.py and the SAINT one by test_parse_saint.py, so between the four files
every supported input format is exercised through to the SAINT and CompPASS inputs.
"""

import numpy as np
import pandas as pd
import pytest

import parse
from experimental_design import ExperimentalDesign


EXPERIMENTS = ["t1_1", "t1_2", "c_1", "c_2"]
N_EXPERIMENTS, N_CONTROLS = 2, 2

# P3 is seen only in the controls and is dropped; P1 and P2 survive.
PROTEINS = ["P1", "P2", "P3"]
QUANT = {
    "P1": {"t1_1": 100.0, "t1_2": 120.0, "c_1": 0.0, "c_2": 5.0},
    "P2": {"t1_1": 40.0, "t1_2": 0.0, "c_1": 0.0, "c_2": 0.0},
    "P3": {"t1_1": 0.0, "t1_2": 0.0, "c_1": 60.0, "c_2": 55.0},
}
SURVIVING_PROTEINS = 2


@pytest.fixture
def ed_file(tmp_path):
    path = tmp_path / "ED.csv"
    pd.DataFrame([
        {"Experiment Name": "t1_1", "Type": "T", "Bait": "BaitA", "Replicate": 1,
         "Bait ID": "P_A"},
        {"Experiment Name": "t1_2", "Type": "T", "Bait": "BaitA", "Replicate": 2,
         "Bait ID": "P_A"},
        {"Experiment Name": "c_1", "Type": "C", "Bait": "Ctrl", "Replicate": 1,
         "Bait ID": "P_C"},
        {"Experiment Name": "c_2", "Type": "C", "Bait": "Ctrl", "Replicate": 2,
         "Bait ID": "P_C"},
    ]).to_csv(path, index=False)
    return str(path)


@pytest.fixture
def maxquant_file(tmp_path):
    frame = pd.DataFrame({
        "Majority protein IDs": PROTEINS,
        "Gene names": ["G_{}".format(p) for p in PROTEINS],
        "Reverse": ["-"] * 3,
        "Only identified by site": ["-"] * 3,
        "Potential contaminant": ["-"] * 3,
        "Sequence length": [100, 200, 300],
    })
    for experiment in EXPERIMENTS:
        frame["Intensity {}".format(experiment)] = [QUANT[p][experiment]
                                                    for p in PROTEINS]
    path = tmp_path / "proteinGroups.txt"
    frame.to_csv(path, sep="\t", index=False)
    return str(path)


@pytest.fixture
def diann_file(tmp_path):
    frame = pd.DataFrame({
        "Protein.Group": PROTEINS,
        "Protein.Names": ["{}_HUMAN".format(p) for p in PROTEINS],
        "Genes": ["G_{}".format(p) for p in PROTEINS],
    })
    for experiment in EXPERIMENTS:
        # A zero reaches DIA-NN's matrix as a blank, which the converter fills.
        frame[experiment] = [QUANT[p][experiment] or np.nan for p in PROTEINS]
    path = tmp_path / "report.pg_matrix.tsv"
    frame.to_csv(path, sep="\t", index=False)
    return str(path)


@pytest.fixture
def pioneer_file(tmp_path):
    frame = pd.DataFrame({
        "gene_names": ["G_{}".format(p) for p in PROTEINS],
        "protein": PROTEINS,
        "target": [True] * 3,
    })
    for experiment in EXPERIMENTS:
        # A zero reaches Pioneer's wide table as an empty cell, which the converter fills.
        frame[experiment] = [QUANT[p][experiment] or np.nan for p in PROTEINS]
    path = tmp_path / "protein_groups_wide.tsv"
    frame.to_csv(path, sep="	", index=False)
    return str(path)


@pytest.fixture
def fragpipe_file(tmp_path):
    frame = pd.DataFrame({
        "Protein": ["sp|{}|X_HUMAN".format(p) for p in PROTEINS],
        "Protein ID": PROTEINS,
        "Gene": ["G_{}".format(p) for p in PROTEINS],
        "Protein Length": [100, 200, 300],
        "Description": ["desc"] * 3,
    })
    for experiment in EXPERIMENTS:
        frame["{} Intensity".format(experiment)] = [QUANT[p][experiment]
                                                    for p in PROTEINS]
    path = tmp_path / "combined_protein.tsv"
    frame.to_csv(path, sep="\t", index=False)
    return str(path)


@pytest.fixture
def run(tmp_path, ed_file, maxquant_file, diann_file, pioneer_file, fragpipe_file):
    """Run one entry point into its own output directory and return that directory."""
    def _run(fmt):
        out = tmp_path / "out_{}".format(fmt)
        if fmt == "maxquant":
            counts = parse.parse_ed_pg(maxquant_file, ed_file, "Intensity", str(out))
        elif fmt == "diann":
            counts = parse.parse_diann(diann_file, ed_file, "Intensity", str(out))
        elif fmt == "pioneer":
            counts = parse.parse_pioneer(pioneer_file, ed_file, "Intensity", str(out))
        elif fmt == "fragpipe":
            counts = parse.parse_fragpipe(fragpipe_file, ed_file, "Intensity", str(out))
        else:
            raise AssertionError("unknown format: {}".format(fmt))
        return out, counts

    return _run


FORMATS = ["maxquant", "diann", "pioneer", "fragpipe"]


# --- outputs -------------------------------------------------------------------

@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("filename", list(parse.PARSE_OUTPUTS) + ["proteinGroups.txt"])
def test_the_scoring_inputs_are_written(run, fmt, filename):
    out, _ = run(fmt)

    assert (out / filename).exists()


@pytest.mark.parametrize("fmt", FORMATS)
def test_the_reported_counts_exclude_controls(run, fmt):
    _, counts = run(fmt)

    assert counts == (N_EXPERIMENTS, N_CONTROLS)


@pytest.mark.parametrize("fmt", FORMATS)
def test_the_counts_match_the_design_that_was_written(run, fmt):
    """Scoring reads the copied ED.csv, so the counts reported at parse time have to be
    the counts that file yields."""
    out, counts = run(fmt)
    ed = ExperimentalDesign(str(out / "ED.csv"))

    assert counts == (ed.num_experiments, ed.num_controls)


@pytest.mark.parametrize("fmt", FORMATS)
def test_control_only_proteins_are_dropped(run, fmt):
    out, _ = run(fmt)
    prey = pd.read_csv(out / "prey.txt", sep="\t", header=None, names=["Prey", "Gene"])

    assert set(prey["Prey"]) == {"P1", "P2"}


@pytest.mark.parametrize("fmt", FORMATS)
def test_the_interaction_file_is_dense_across_every_experiment(run, fmt):
    out, _ = run(fmt)
    interaction = pd.read_csv(out / "interaction.txt", sep="\t", header=None,
                              names=["Experiment", "Bait", "Prey", "Intensity"])

    assert len(interaction) == SURVIVING_PROTEINS * len(EXPERIMENTS)
    assert set(interaction["Experiment"]) == set(EXPERIMENTS)


@pytest.mark.parametrize("fmt", FORMATS)
def test_the_bait_file_covers_every_experiment(run, fmt):
    out, _ = run(fmt)
    bait = pd.read_csv(out / "bait.txt", sep="\t", header=None,
                       names=["Experiment", "Bait", "Type"])

    assert list(bait["Experiment"]) == EXPERIMENTS
    assert list(bait["Type"]) == ["T", "T", "C", "C"]


@pytest.mark.parametrize("fmt", FORMATS)
def test_comppass_input_follows_the_bait_name_and_id_convention(run, fmt):
    out, _ = run(fmt)
    written = pd.read_csv(out / "to_CompPASS.csv")

    assert set(written["Experiment.ID"]) <= {"BaitA", "Ctrl"}
    assert set(written["Bait"]) <= {"P_A", "P_C"}


def test_every_format_produces_the_same_scoring_inputs(run):
    """The converters differ only in what they read.  A prey universe that varied
    by input format would make results incomparable between them."""
    preys = {}
    for fmt in FORMATS:
        out, _ = run(fmt)
        prey = pd.read_csv(out / "prey.txt", sep="\t", header=None, names=["Prey", "G"])
        preys[fmt] = set(prey["Prey"])

    assert preys["maxquant"] == preys["diann"] == preys["pioneer"] == preys["fragpipe"]


# --- copies of the inputs ------------------------------------------------------

def test_maxquant_copies_the_protein_groups_file_verbatim(run, maxquant_file):
    out, _ = run("maxquant")

    assert (out / "proteinGroups.txt").read_bytes() == open(maxquant_file, "rb").read()


@pytest.mark.parametrize("fmt", ["diann", "fragpipe"])
def test_the_converted_frame_is_saved_as_protein_groups(run, fmt):
    """The vendor file itself is not copied; what is saved is the MaxQuant-shaped frame
    that was actually parsed, which is the one a later run would have to reproduce."""
    out, _ = run(fmt)
    saved = pd.read_csv(out / "proteinGroups.txt", sep="\t")

    assert "Majority protein IDs" in saved.columns
    assert "Intensity t1_1" in saved.columns


@pytest.mark.parametrize("fmt", FORMATS)
def test_the_design_is_copied_unchanged(run, fmt, ed_file):
    out, _ = run(fmt)

    assert (out / "ED.csv").read_bytes() == open(ed_file, "rb").read()


@pytest.mark.parametrize("fmt", FORMATS)
def test_a_dataset_log_is_written(run, fmt):
    out, _ = run(fmt)

    assert (out / "proximate.log").exists()


@pytest.mark.parametrize("fmt", FORMATS)
def test_the_run_is_recorded_in_the_manifest(run, fmt):
    out, _ = run(fmt)

    assert (out / "run.json").exists()


# --- format-specific behavior --------------------------------------------------

def test_diann_ignores_the_requested_quantification(tmp_path, ed_file, diann_file):
    """Documented, not fixed: parse_diann accepts quantType and always builds
    ProteinGroups with "Intensity".  Asking for LFQ therefore yields an intensity parse
    rather than an error, though the request is still recorded in the manifest."""
    as_lfq = tmp_path / "out_lfq"
    as_intensity = tmp_path / "out_int"

    parse.parse_diann(diann_file, ed_file, "LFQ", str(as_lfq))
    parse.parse_diann(diann_file, ed_file, "Intensity", str(as_intensity))

    assert ((as_lfq / "interaction.txt").read_bytes()
            == (as_intensity / "interaction.txt").read_bytes())


def test_pioneer_ignores_the_requested_quantification(tmp_path, ed_file, pioneer_file):
    """As for DIA-NN: quantType is recorded but the parse is always by intensity."""
    as_lfq = tmp_path / "out_lfq"
    as_intensity = tmp_path / "out_int"

    parse.parse_pioneer(pioneer_file, ed_file, "LFQ", str(as_lfq))
    parse.parse_pioneer(pioneer_file, ed_file, "Intensity", str(as_intensity))

    assert ((as_lfq / "interaction.txt").read_bytes()
            == (as_intensity / "interaction.txt").read_bytes())


def test_fragpipe_honors_the_requested_quantification(tmp_path, ed_file):
    """Unlike DIA-NN, the quantification reaches both the converter and ProteinGroups."""
    frame = pd.DataFrame({
        "Protein": ["sp|P1|X_HUMAN"],
        "Protein ID": ["P1"],
        "Gene": ["G_P1"],
        "Protein Length": [100],
    })
    for experiment in EXPERIMENTS:
        frame["{} Intensity".format(experiment)] = [10.0]
        frame["{} Total Spectral Count".format(experiment)] = [7]
    fp = tmp_path / "combined_protein.tsv"
    frame.to_csv(fp, sep="\t", index=False)

    out = tmp_path / "out_spc"
    parse.parse_fragpipe(str(fp), ed_file, "Spectral Counts", str(out))
    interaction = pd.read_csv(out / "interaction.txt", sep="\t", header=None,
                              names=["Experiment", "Bait", "Prey", "Value"])

    assert set(interaction["Value"]) == {7}


def test_a_design_experiment_missing_from_the_data_is_rejected(tmp_path, ed_file,
                                                              maxquant_file):
    """Validation runs before anything is written, so the output directory is left with
    no half-built scoring inputs."""
    from ed_exceptions import EDPGMismatchError

    frame = pd.read_csv(maxquant_file, sep="\t").drop(columns=["Intensity c_2"])
    incomplete = tmp_path / "incomplete.txt"
    frame.to_csv(incomplete, sep="\t", index=False)
    out = tmp_path / "out_bad"

    with pytest.raises(EDPGMismatchError):
        parse.parse_ed_pg(str(incomplete), ed_file, "Intensity", str(out))

    assert not (out / "prey.txt").exists()
