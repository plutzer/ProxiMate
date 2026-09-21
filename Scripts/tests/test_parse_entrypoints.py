"""End-to-end tests for the parse entry points.

Each writes the files the scoring stage reads and reports the experiment counts that
reach run.json and the GUI.  Every supported input format is exercised through to the
SAINT and CompPASS inputs, and the formats are checked against each other.
"""

import json

import numpy as np
import pandas as pd
import pytest

import parse
from ed_exceptions import EDPGMismatchError


EXPERIMENTS = ["t1_1", "t1_2", "c_1", "c_2"]
N_EXPERIMENTS, N_CONTROLS = 2, 2

# P3 is seen only in the controls and is dropped; P1 and P2 survive.
PROTEINS = ["P1", "P2", "P3"]
QUANT = {
    "P1": {"t1_1": 128.0, "t1_2": 64.0, "c_1": 0.0, "c_2": 8.0},
    "P2": {"t1_1": 32.0, "t1_2": 0.0, "c_1": 0.0, "c_2": 0.0},
    "P3": {"t1_1": 0.0, "t1_2": 0.0, "c_1": 16.0, "c_2": 4.0},
}
SURVIVING_PROTEINS = 2

INTERACTION_COLUMNS = ["Experiment", "Bait", "Prey", "Value"]


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
        frame[experiment] = [QUANT[p][experiment] or np.nan for p in PROTEINS]
    path = tmp_path / "protein_groups_wide.tsv"
    frame.to_csv(path, sep="\t", index=False)
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
def msstats_file(tmp_path):
    """The same quantification in MSstats long form: log2 values, unobserved cells absent."""
    rows = [
        {"Protein": p, "originalRUN": e, "LABEL": "L", "GROUP": "g", "SUBJECT": "s",
         "LogIntensities": np.log2(QUANT[p][e])}
        for p in PROTEINS for e in EXPERIMENTS if QUANT[p][e] > 0
    ]
    path = tmp_path / "ProteinLevelData.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)


@pytest.fixture
def run(tmp_path, ed_file, maxquant_file, diann_file, pioneer_file, fragpipe_file,
        msstats_file):
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
        elif fmt == "msstats":
            counts = parse.parse_msstats(msstats_file, ed_file, str(out))
        else:
            raise AssertionError("unknown format: {}".format(fmt))
        return out, counts

    return _run


FORMATS = ["maxquant", "diann", "pioneer", "fragpipe", "msstats"]


def _interaction(out):
    return pd.read_csv(out / "interaction.txt", sep="\t", header=None,
                       names=INTERACTION_COLUMNS)


# --- outputs ---------------------------------------------------------------------

@pytest.mark.parametrize("fmt", FORMATS)
def test_the_scoring_inputs_are_written(run, fmt):
    out, _ = run(fmt)

    for filename in list(parse.PARSE_OUTPUTS) + ["proteinGroups.txt"]:
        assert (out / filename).exists(), filename


@pytest.mark.parametrize("fmt", FORMATS)
def test_the_reported_counts_exclude_controls(run, fmt):
    _, counts = run(fmt)

    assert counts == (N_EXPERIMENTS, N_CONTROLS)


@pytest.mark.parametrize("fmt", FORMATS)
def test_control_only_proteins_are_dropped(run, fmt):
    out, _ = run(fmt)
    prey = pd.read_csv(out / "prey.txt", sep="\t", header=None, names=["Prey", "Gene"])

    assert set(prey["Prey"]) == {"P1", "P2"}


def test_every_format_produces_the_same_scoring_inputs(run):
    """The converters differ only in what they read; a value that varied by input
    format would make results incomparable between them.  The interaction file is
    dense, so this also checks every surviving prey against every experiment."""
    tables = {}
    for fmt in FORMATS:
        out, _ = run(fmt)
        table = _interaction(out).sort_values(["Experiment", "Prey"]).reset_index(drop=True)
        table["Value"] = table["Value"].round(6)
        tables[fmt] = table

    reference = tables["maxquant"]
    assert len(reference) == SURVIVING_PROTEINS * len(EXPERIMENTS)
    for fmt in FORMATS[1:]:
        pd.testing.assert_frame_equal(tables[fmt], reference, check_dtype=False)


# --- format-specific behavior ----------------------------------------------------

def test_fragpipe_honors_the_requested_quantification(tmp_path, ed_file):
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

    assert set(_interaction(out)["Value"]) == {7}


def test_a_design_experiment_missing_from_the_data_is_rejected(tmp_path, ed_file,
                                                              maxquant_file):
    """Validation runs before anything is written, so the output directory is left with
    no half-built scoring inputs."""
    frame = pd.read_csv(maxquant_file, sep="\t").drop(columns=["Intensity c_2"])
    incomplete = tmp_path / "incomplete.txt"
    frame.to_csv(incomplete, sep="\t", index=False)
    out = tmp_path / "out_bad"

    with pytest.raises(EDPGMismatchError):
        parse.parse_ed_pg(str(incomplete), ed_file, "Intensity", str(out))

    assert not (out / "prey.txt").exists()


# --- SAINT-format input ------------------------------------------------------------

@pytest.fixture
def saint_inputs(tmp_path):
    """Two test baits and one control, two replicates each, two preys."""
    bait_df = pd.DataFrame({
        "Experiment Name": ["t1_1", "t1_2", "t2_1", "t2_2", "c_1", "c_2"],
        "Bait": ["BaitA", "BaitA", "BaitB", "BaitB", "Ctrl", "Ctrl"],
        "Type": ["T", "T", "T", "T", "C", "C"],
    })
    bait_df["Bait ID"] = "None"

    prey = tmp_path / "prey.txt"
    prey.write_text("P1\tGene1\nP2\tGene2\n")

    interaction = tmp_path / "interaction.txt"
    interaction.write_text("".join(
        "{}\t{}\t{}\t{}\n".format(e, b, p, v)
        for e, b in zip(bait_df["Experiment Name"], bait_df["Bait"])
        for p, v in (("P1", 10), ("P2", 5))))

    return bait_df, str(prey), str(interaction)


def test_saint_experiment_counts_exclude_controls(saint_inputs, tmp_path):
    bait_df, prey, interaction = saint_inputs

    counts = parse.parse_from_saint(bait_df, prey, interaction, str(tmp_path / "out"))

    assert counts == (4, 2)


def test_saint_input_is_rebuilt_into_a_comppass_table(saint_inputs, tmp_path):
    """Replicates are numbered by order within each bait, and every interaction row is
    joined to its prey name and its bait's protein ID under the CompPASS convention."""
    bait_df, prey, interaction = saint_inputs
    out = tmp_path / "out"

    parse.parse_from_saint(bait_df, prey, interaction, str(out))
    written = pd.read_csv(out / "to_CompPASS.csv", keep_default_na=False)

    assert list(written.columns) == [
        "Experiment.ID", "Replicate", "Bait", "Prey", "Prey.Name", "Spectral.Count"]
    assert len(written) == 12
    by_bait = written.groupby("Experiment.ID")["Replicate"].apply(sorted).to_dict()
    assert by_bait == {"BaitA": [1, 1, 2, 2], "BaitB": [1, 1, 2, 2], "Ctrl": [1, 1, 2, 2]}
    assert set(written["Bait"]) == {"None"}
    assert dict(zip(written["Prey"], written["Prey.Name"])) == {"P1": "Gene1", "P2": "Gene2"}
    assert set(written[written["Prey"] == "P1"]["Spectral.Count"]) == {10}


# --- helpers ---------------------------------------------------------------------

@pytest.mark.parametrize("name, taken, valid", [
    ("dataset_01", [], True),
    ("", [], False),
    ("has space", [], False),
    ("has-hyphen", [], False),
    ("existing", ["existing"], False),
])
def test_dataset_names_are_alphanumeric_underscore_and_unused(name, taken, valid):
    """Success is the int 0; anything else is a message describing the problem."""
    result = parse.validate_name(name, taken)

    assert (result == 0) is valid


@parse._parse_stage
def _records_two_counts(proteinGroups, quantType, outputPath):
    return 3, 2


@parse._parse_stage
def _fails(proteinGroups, outputPath):
    raise RuntimeError("parsing blew up")


def _stage_entry(out_dir):
    """The single stage recorded in the manifest at `out_dir`."""
    document = json.loads((out_dir / "run.json").read_text())
    runs = list(document["runs"].values())
    assert len(runs) == 1 and len(runs[0]["stages"]) == 1
    return runs[0]["stages"][0]


def test_a_parse_stage_records_its_inputs_parameters_and_counts(tmp_path):
    """A checksum of each file is what makes a run reproducible; a value parameter
    has nothing to check and is recorded as is."""
    input_file = tmp_path / "proteinGroups.txt"
    input_file.write_text("Majority protein IDs\nP1\n")
    out = tmp_path / "out"

    assert _records_two_counts(str(input_file), "LFQ", str(out)) == (3, 2)
    entry = _stage_entry(out)

    assert entry["entrypoint"] == "parse._records_two_counts"
    assert entry["metrics"] == {"n_experiments": 3, "n_controls": 2}
    assert [i["role"] for i in entry["inputs"]] == ["proteinGroups"]
    assert entry["params"] == {"proteinGroups": str(input_file), "quantType": "LFQ"}


def test_a_failing_parse_still_leaves_a_manifest_and_reraises(tmp_path):
    out = tmp_path / "out"

    with pytest.raises(RuntimeError, match="parsing blew up"):
        _fails(str(tmp_path / "proteinGroups.txt"), str(out))

    entry = _stage_entry(out)
    assert entry["status"] == "error"
    assert entry["error"]["type"] == "RuntimeError"
