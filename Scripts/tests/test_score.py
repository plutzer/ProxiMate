"""Tests for the scoring stage's orchestration.

SAINTexpress is a compiled binary that exists only in the container, so these tests stub
the one function that invokes it and exercise everything around it: which of the three
installed builds a run selects, which prey file it is handed, how a grouped run's inputs
are cut down, and the column schema the merged output has to keep for annotation to read
it.
"""

import argparse
import os

import numpy as np
import pandas as pd
import pytest

import score


# --- binary and prey file selection --------------------------------------------

@pytest.mark.parametrize("quant_type, imputation, binary, prey_file", [
    ("Intensity", "0", score.SAINT_EXPRESS_INT_DEFAULT_DIR, "prey.txt"),
    ("Intensity", "1", score.SAINT_EXPRESS_INT_DIR, "imputed_prey.txt"),
    ("Intensity", "2", score.SAINT_EXPRESS_INT_DIR, "imputed_prey.txt"),
    ("Intensity", "3", score.SAINT_EXPRESS_INT_DIR, "imputed_prey.txt"),
    ("LFQ", "0", score.SAINT_EXPRESS_INT_DEFAULT_DIR, "prey.txt"),
    ("LFQ", "1", score.SAINT_EXPRESS_INT_DIR, "imputed_prey.txt"),
    ("Spectral Counts", "0", score.SAINT_EXPRESS_SPC_DIR, "prey.txt"),
])
def test_quant_type_and_imputation_select_the_build_and_prey_file(
        quant_type, imputation, binary, prey_file):
    """Which build ran is not recoverable from list.txt, so the selection is pinned."""
    assert score._select_saint(quant_type, imputation) == (binary, prey_file)


def test_the_command_names_its_inputs_relative_to_the_working_directory():
    """SAINTexpress is invoked with cwd set to the group or score directory, so absolute
    paths here would read the wrong group's files."""
    cmd = score._build_saint_cmd(1000, "Intensity", "2")

    assert cmd == [score.SAINT_EXPRESS_INT_DIR, "-L", "1000", "filtered_interaction.txt",
                   "imputed_prey.txt", "bait.txt"]


# --- _run_saint ----------------------------------------------------------------

class _FakeCompleted:
    def __init__(self, returncode=0, stdout="", stderr=""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


@pytest.fixture
def fake_saint(monkeypatch):
    """Replace the SAINTexpress subprocess, recording the call and writing list.txt."""
    calls = []

    def _install(returncode=0, write_list=True):
        def _run(cmd, cwd=None, capture_output=False, text=False):
            calls.append({"cmd": cmd, "cwd": cwd})
            if write_list and returncode == 0:
                _write_saint_list(cwd)
            return _FakeCompleted(returncode, "out", "err")

        monkeypatch.setattr(score.subprocess, "run", _run)
        return calls

    return _install


SAINT_COLUMNS = ["Bait", "Prey", "PreyGene", "AvgSpec", "AvgP", "SaintScore", "BFDR",
                 "boosted_by"]

TEST_BAITS = ["BaitA", "BaitB", "BaitC"]
PREYS = ["P1", "P2", "P3", "P4"]

# AvgP per bait and prey in the fake SAINT output.  No two values tie, so the pooled
# BFDR of a row depends on which other baits were scored alongside it.
AVGP = {
    "BaitA": [0.9, 0.7, 0.5, 0.3],
    "BaitB": [0.8, 0.6, 0.4, 0.2],
    "BaitC": [0.95, 0.85, 0.75, 0.65],
}

# BFDR the fake writes; any real value would be recomputed from AvgP.
FAKE_BFDR = 0.5


def _write_saint_list(cwd):
    """One row per test bait and prey, which is the shape SAINTexpress reports.

    Only baits present in the directory's bait.txt are reported, so a grouped run's
    per-group list.txt covers that group's baits alone.
    """
    baits = TEST_BAITS
    bait_path = os.path.join(cwd, "bait.txt")
    if os.path.exists(bait_path):
        bait_file = pd.read_csv(bait_path, sep="\t", header=None, dtype=str)
        baits = [b for b in TEST_BAITS if b in set(bait_file[1])]
    rows = [(bait, prey, "G_" + prey, 10.0 * (i + 1), AVGP[bait][i], AVGP[bait][i],
             FAKE_BFDR, None)
            for bait in baits
            for i, prey in enumerate(PREYS)]
    pd.DataFrame(rows, columns=SAINT_COLUMNS).to_csv(
        os.path.join(cwd, "list.txt"), sep="\t", index=False)


def test_a_failing_saint_run_stops_the_pipeline(tmp_path, fake_saint):
    fake_saint(returncode=1)

    with pytest.raises(SystemExit) as excinfo:
        score._run_saint(str(tmp_path), 1000, "Intensity", "0")

    assert excinfo.value.code == 1


def test_a_silent_failure_to_produce_results_is_caught(tmp_path, fake_saint):
    """SAINTexpress can exit zero having written nothing; continuing would merge an
    absent list.txt and report no interactions rather than an error."""
    fake_saint(write_list=False)

    with pytest.raises(SystemExit) as excinfo:
        score._run_saint(str(tmp_path), 1000, "Intensity", "0")

    assert excinfo.value.code == 1


# --- _build_group_saint_inputs -------------------------------------------------

class _FakeExperiment:
    def __init__(self, name):
        self.attributes = {"Experiment Name": name}


class _FakeED:
    def __init__(self, tests, controls):
        self._tests = [_FakeExperiment(n) for n in tests]
        self._controls = [_FakeExperiment(n) for n in controls]

    def get_experiments_for_group(self, group):
        return self._tests, self._controls


@pytest.fixture
def group_inputs(tmp_path):
    """A source directory holding the global SAINT inputs, and an empty group directory."""
    src = tmp_path / "src"
    src.mkdir()
    group = tmp_path / "group_1"
    group.mkdir()

    def _write(names, prey_lines="P1\tG1\nP2\tG2\n"):
        with open(src / "bait.txt", "w", newline="") as handle:
            for name in names:
                handle.write("{}\tBait_{}\tT\n".format(name, name))
        with open(src / "filtered_interaction.txt", "w", newline="") as handle:
            for name in names:
                for prey in ("P1", "P2"):
                    handle.write("{}\tBait_{}\t{}\t100.0\n".format(name, name, prey))
        (src / "prey.txt").write_text(prey_lines)
        return src, group

    return _write


def _read(path, names):
    return pd.read_csv(path, sep="\t", header=None, names=names, dtype=str)


def test_only_the_groups_experiments_are_kept(group_inputs):
    src, group = group_inputs(["t1_1", "t1_2", "t2_1", "c_1"])
    ed = _FakeED(["t1_1", "t1_2"], ["c_1"])

    score._build_group_saint_inputs(str(src), str(group), ed, 1, False)
    bait = _read(group / "bait.txt", ["Experiment", "Bait", "Type"])

    assert set(bait["Experiment"]) == {"t1_1", "t1_2", "c_1"}


def test_the_interaction_file_is_filtered_to_the_same_experiments(group_inputs):
    src, group = group_inputs(["t1_1", "t2_1", "c_1"])
    ed = _FakeED(["t1_1"], ["c_1"])

    score._build_group_saint_inputs(str(src), str(group), ed, 1, False)
    interaction = _read(group / "filtered_interaction.txt",
                        ["Experiment", "Bait", "Prey", "Intensity"])

    assert set(interaction["Experiment"]) == {"t1_1", "c_1"}
    assert len(interaction) == 4


@pytest.mark.parametrize("design_names", [(["1"], ["3"]), ([1], [3])],
                         ids=["design-strings", "design-ints"])
def test_numeric_experiment_names_are_matched_as_text(group_inputs, design_names):
    """pandas infers int64 for a column of digits.  Without the explicit dtype on the
    files, the filter compares integers against the design's strings, matches nothing,
    and every group is scored on an empty file.  A design supplying names as numbers
    must match the files just the same."""
    src, group = group_inputs(["1", "2", "3"])
    ed = _FakeED(*design_names)

    score._build_group_saint_inputs(str(src), str(group), ed, 1, False)
    bait = _read(group / "bait.txt", ["Experiment", "Bait", "Type"])

    assert set(bait["Experiment"]) == {"1", "3"}


def test_the_prey_universe_is_copied_unchanged(group_inputs):
    """Every group must see the same preys, or the pooled BFDR recomputation is not
    comparing like with like."""
    src, group = group_inputs(["t1_1", "t2_1", "c_1"])
    ed = _FakeED(["t1_1"], ["c_1"])

    score._build_group_saint_inputs(str(src), str(group), ed, 1, False)

    assert (group / "prey.txt").read_bytes() == (src / "prey.txt").read_bytes()


# --- _score --------------------------------------------------------------------

# Three baits in duplicate plus a duplicate control.  CompPASS scores a prey against the
# spread of its values across baits, so a one-bait design gives it nothing to compare.
ED_ROWS = [
    {"Experiment Name": "{}_{}".format(bait.lower(), replicate), "Type": "T",
     "Bait": bait, "Replicate": replicate, "Bait ID": "P_{}".format(bait[-1])}
    for bait in TEST_BAITS for replicate in (1, 2)
] + [
    {"Experiment Name": "c_{}".format(replicate), "Type": "C", "Bait": "Ctrl",
     "Replicate": replicate, "Bait ID": "P_Ctrl"}
    for replicate in (1, 2)
]

# Prey counts per bait, chosen so P1 is specific to BaitA and P4 is seen everywhere.
COUNTS = {
    "BaitA": {"P1": 40, "P2": 20, "P3": 4, "P4": 12},
    "BaitB": {"P1": 3, "P2": 30, "P3": 6, "P4": 14},
    "BaitC": {"P1": 2, "P2": 5, "P3": 35, "P4": 13},
    "Ctrl": {"P1": 1, "P2": 2, "P3": 1, "P4": 11},
}


@pytest.fixture
def score_inputs(tmp_path):
    """Lay out the four files _score requires, plus the design, and return an args
    namespace pointing at them."""
    def _build(ed_rows=ED_ROWS, group_column=None):
        work = tmp_path / "score"
        work.mkdir(exist_ok=True)

        rows = [dict(row) for row in ed_rows]
        if group_column is not None:
            for row, group in zip(rows, group_column):
                row["Group"] = group
        ed_path = work / "ED.csv"
        pd.DataFrame(rows).to_csv(ed_path, index=False)

        (work / "prey.txt").write_text("P1\tG1\nP2\tG2\n")
        with open(work / "bait.txt", "w", newline="") as handle:
            for row in rows:
                handle.write("{}\t{}\t{}\n".format(
                    row["Experiment Name"], row["Bait"], row["Type"]))
        with open(work / "interaction.txt", "w", newline="") as handle:
            for row in rows:
                for prey in PREYS:
                    handle.write("{}\t{}\t{}\t{}\n".format(
                        row["Experiment Name"], row["Bait"], prey,
                        float(COUNTS[row["Bait"]][prey])))

        comppass = []
        for row in rows:
            for prey in PREYS:
                comppass.append({
                    "Experiment.ID": row["Bait"], "Replicate": row["Replicate"],
                    "Bait": row["Bait ID"], "Prey": prey, "Prey.Name": prey,
                    "Spectral.Count": COUNTS[row["Bait"]][prey]})
        pd.DataFrame(comppass).to_csv(work / "to_CompPASS.csv", index=False)

        return argparse.Namespace(
            scoreInputs=str(work), outputPath=str(work),
            experimentalDesign=str(ed_path), imputation="0", quantType="Intensity",
            compress_n_rep=1000, n_iterations=2, seed=0,
            pi_method="weighted_average", pi_bait=None)

    return _build


class _Record:
    """Stands in for a provenance record; every method only accumulates."""

    def __init__(self):
        self.inputs, self.outputs, self.metrics, self.extras = [], [], {}, {}

    def add_input(self, path, role=None):
        self.inputs.append((path, role))

    def add_output(self, path, rows=None, role=None):
        self.outputs.append((path, rows))

    def metric(self, key, value):
        self.metrics[key] = value

    def extra(self, **fields):
        self.extras.update(fields)


@pytest.fixture
def run_score(fake_saint):
    def _run(args):
        calls = fake_saint()
        record = _Record()
        score._score(args, record)
        return record, calls

    return _run


def test_the_merged_output_keeps_the_columns_annotation_reads(score_inputs, run_score):
    """Annotation joins on these names, so a rename here surfaces as missing annotations
    rather than as an error.  The schema is the same whether or not the run was
    grouped: a legacy run carries an all-NA source_group column."""
    args = score_inputs()
    record, _ = run_score(args)

    merged = pd.read_csv(os.path.join(args.outputPath, "merged.csv"))

    assert record.extras["run_mode"] == "legacy"
    assert {"Experiment.ID", "Prey.ID", "Bait.ID", "source_group"} <= set(merged.columns)
    assert not {"Bait_x", "Bait_y"} & set(merged.columns)
    assert merged["source_group"].isna().all()


def test_the_bait_name_and_its_protein_id_end_up_in_separate_columns(score_inputs,
                                                                    run_score):
    """SAINT's own Bait column is dropped as redundant: the bait name survives as
    Experiment.ID, which is what the GUI groups on, and Bait.ID carries the protein ID
    CompPASS supplied."""
    args = score_inputs()
    run_score(args)

    merged = pd.read_csv(os.path.join(args.outputPath, "merged.csv"))

    assert set(merged["Experiment.ID"]) == set(TEST_BAITS)
    assert set(merged["Bait.ID"].dropna()) == {"P_A", "P_B", "P_C"}
    assert "Bait" not in merged.columns


def test_a_missing_required_input_stops_the_run(score_inputs, run_score):
    args = score_inputs()
    os.remove(os.path.join(args.scoreInputs, "bait.txt"))

    with pytest.raises(SystemExit) as excinfo:
        run_score(args)

    assert excinfo.value.code == 1


def _pooled_bfdr(avgp):
    """SAINT's BFDR definition over one pooled set of AvgP values."""
    avgp = np.asarray(avgp)
    return np.array([1.0 - avgp[avgp > p].mean() if (avgp > p).any() else 0.0
                     for p in avgp])


def test_a_grouped_run_scores_each_group_and_repools_bfdr(score_inputs, run_score):
    """One SAINTexpress run per group, rows labelled with their group, and BFDR
    recomputed over the merged table rather than kept from each group's own run."""
    args = score_inputs(group_column=["1"] * 4 + ["2"] * 2 + ["*", "*"])

    record, calls = run_score(args)
    merged = pd.read_csv(os.path.join(args.outputPath, "merged.csv"))

    assert record.extras["run_mode"] == "grouped"
    assert [os.path.basename(c["cwd"]) for c in calls] == ["1", "2"]
    assert merged.set_index("Experiment.ID")["source_group"].to_dict() == {
        "BaitA": 1, "BaitB": 1, "BaitC": 2}
    assert not (merged["BFDR"] == FAKE_BFDR).any()
    assert np.allclose(merged["BFDR"], _pooled_bfdr(merged["AvgP"]))
