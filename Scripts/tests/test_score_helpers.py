"""Tests for the scoring stage's orchestration.

SAINTexpress is a compiled binary that exists only in the container, so these tests stub
the one function that invokes it and exercise everything around it: which of the three
installed builds a run selects, which prey file it is handed, how a grouped run's inputs
are cut down, and the column schema the merged output has to keep for annotation to read
it.
"""

import argparse
import os

import pandas as pd
import pytest

import score


# --- _choose_prey_filename -----------------------------------------------------

@pytest.mark.parametrize("imputation", ["1", "2", "3"])
def test_an_imputed_intensity_run_reads_the_imputed_prey_file(imputation):
    assert score._choose_prey_filename("Intensity", imputation) == "imputed_prey.txt"


@pytest.mark.parametrize("imputation", ["0", "", None])
def test_an_unimputed_run_reads_the_plain_prey_file(imputation):
    assert score._choose_prey_filename("Intensity", imputation) == "prey.txt"


@pytest.mark.parametrize("imputation", ["0", "1", "2", "3"])
def test_spectral_counts_never_reads_an_imputed_prey_file(imputation):
    """Imputation is not implemented for spectral counts, so there is no imputed prey
    file to read whatever was asked for."""
    assert score._choose_prey_filename("Spectral Counts", imputation) == "prey.txt"


def test_the_imputation_flag_is_compared_as_a_string():
    """Documented, not fixed: argparse leaves --imputation a string, and an int 1 built
    by hand takes the unimputed branch instead of raising."""
    assert score._choose_prey_filename("Intensity", 1) == "prey.txt"
    assert score._choose_prey_filename("Intensity", "1") == "imputed_prey.txt"


# --- _build_saint_cmd ----------------------------------------------------------

@pytest.mark.parametrize("imputation", ["1", "2", "3"])
def test_an_imputed_intensity_run_uses_the_intensity_build(imputation):
    cmd = score._build_saint_cmd(1000, "Intensity", imputation, "imputed_prey.txt")

    assert cmd[0] == score.SAINT_EXPRESS_INT_DIR


@pytest.mark.parametrize("imputation", ["0", "", None])
def test_an_unimputed_intensity_run_uses_the_default_build(imputation):
    """Three SAINTexpress builds are installed and they do not agree; which one ran is
    not recoverable from list.txt, which is why the command is recorded."""
    cmd = score._build_saint_cmd(1000, "Intensity", imputation, "prey.txt")

    assert cmd[0] == score.SAINT_EXPRESS_INT_DEFAULT_DIR


def test_a_spectral_counts_run_uses_the_spectral_counts_build():
    cmd = score._build_saint_cmd(1000, "Spectral Counts", "0", "prey.txt")

    assert cmd[0] == score.SAINT_EXPRESS_SPC_DIR


def test_the_three_builds_are_distinct():
    assert len({score.SAINT_EXPRESS_INT_DIR, score.SAINT_EXPRESS_INT_DEFAULT_DIR,
                score.SAINT_EXPRESS_SPC_DIR}) == 3


def test_the_command_names_its_inputs_relative_to_the_working_directory():
    """SAINTexpress is invoked with cwd set to the group or score directory, so absolute
    paths here would read the wrong group's files."""
    cmd = score._build_saint_cmd(1000, "Intensity", "2", "imputed_prey.txt")

    assert cmd[1:] == ["-L", "1000", "filtered_interaction.txt", "imputed_prey.txt",
                       "bait.txt"]


def test_the_replicate_compression_is_stringified():
    assert score._build_saint_cmd(4, "Intensity", "0", "prey.txt")[2] == "4"


def test_a_spectral_counts_run_ignores_the_chosen_prey_file():
    """The branch hardcodes prey.txt, so an imputed prey file passed in is discarded
    rather than used."""
    cmd = score._build_saint_cmd(1000, "Spectral Counts", "2", "imputed_prey.txt")

    assert "imputed_prey.txt" not in cmd
    assert "prey.txt" in cmd


def test_imputation_requested_for_spectral_counts_is_reported(caplog):
    """The run proceeds unimputed, which would otherwise look like the imputation had
    been applied and made no difference."""
    with caplog.at_level("WARNING", logger="proximate.score"):
        score._build_saint_cmd(1000, "Spectral Counts", "1", "prey.txt")

    assert "not yet implemented" in caplog.text


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
                _write_saint_list(os.path.join(cwd, "list.txt"))
            return _FakeCompleted(returncode, "out", "err")

        monkeypatch.setattr(score.subprocess, "run", _run)
        return calls

    return _install


SAINT_COLUMNS = ["Bait", "Prey", "PreyGene", "AvgSpec", "AvgP", "SaintScore", "BFDR",
                 "boosted_by"]


TEST_BAITS = ["BaitA", "BaitB", "BaitC"]
PREYS = ["P1", "P2", "P3", "P4"]


def _write_saint_list(path):
    """One row per test bait and prey, which is the shape SAINTexpress reports."""
    rows = [(bait, prey, "G_" + prey, 10.0 * (i + 1), 0.9, 0.9, 0.01, None)
            for bait in TEST_BAITS
            for i, prey in enumerate(PREYS)]
    pd.DataFrame(rows, columns=SAINT_COLUMNS).to_csv(path, sep="\t", index=False)


def test_the_built_command_is_what_runs(tmp_path, fake_saint):
    calls = fake_saint()

    returned = score._run_saint(str(tmp_path), 1000, "Intensity", "0", "prey.txt")

    assert calls[0]["cmd"] == returned
    assert calls[0]["cwd"] == str(tmp_path)


def test_a_failing_saint_run_stops_the_pipeline(tmp_path, fake_saint):
    fake_saint(returncode=1)

    with pytest.raises(SystemExit) as excinfo:
        score._run_saint(str(tmp_path), 1000, "Intensity", "0", "prey.txt")

    assert excinfo.value.code == 1


def test_saint_output_is_reported_when_the_run_fails(tmp_path, fake_saint, caplog):
    fake_saint(returncode=1)

    with caplog.at_level("ERROR", logger="proximate.score"):
        with pytest.raises(SystemExit):
            score._run_saint(str(tmp_path), 1000, "Intensity", "0", "prey.txt")

    assert "err" in caplog.text


def test_a_silent_failure_to_produce_results_is_caught(tmp_path, fake_saint):
    """SAINTexpress can exit zero having written nothing; continuing would merge an
    absent list.txt and report no interactions rather than an error."""
    fake_saint(write_list=False)

    with pytest.raises(SystemExit) as excinfo:
        score._run_saint(str(tmp_path), 1000, "Intensity", "0", "prey.txt")

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

    def _write(names, prey_lines="P1\tG1\nP2\tG2\n", imputed=False):
        with open(src / "bait.txt", "w", newline="") as handle:
            for name in names:
                handle.write("{}\tBait_{}\tT\n".format(name, name))
        with open(src / "filtered_interaction.txt", "w", newline="") as handle:
            for name in names:
                for prey in ("P1", "P2"):
                    handle.write("{}\tBait_{}\t{}\t100.0\n".format(name, name, prey))
        (src / "prey.txt").write_text(prey_lines)
        if imputed:
            (src / "imputed_prey.txt").write_text(prey_lines)
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


def test_numeric_experiment_names_are_read_as_text(group_inputs):
    """pandas infers int64 for a column of digits.  Without the explicit dtype on the
    files, the filter compares integers against the design's strings, matches nothing,
    and every group is scored on an empty file."""
    src, group = group_inputs(["1", "2", "3"])
    ed = _FakeED(["1"], ["3"])

    score._build_group_saint_inputs(str(src), str(group), ed, 1, False)
    bait = _read(group / "bait.txt", ["Experiment", "Bait", "Type"])

    assert set(bait["Experiment"]) == {"1", "3"}


def test_numeric_experiment_names_from_the_design_are_coerced(group_inputs):
    """The other half of the same guard: a design supplying names as numbers rather than
    as the strings a CSV round-trip yields must still match the files."""
    src, group = group_inputs(["1", "2", "3"])
    ed = _FakeED([1], [3])

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


def test_the_imputed_prey_file_is_copied_only_when_it_is_used(group_inputs):
    src, group = group_inputs(["t1_1", "c_1"], imputed=True)
    ed = _FakeED(["t1_1"], ["c_1"])

    score._build_group_saint_inputs(str(src), str(group), ed, 1, False)
    assert not (group / "imputed_prey.txt").exists()

    score._build_group_saint_inputs(str(src), str(group), ed, 1, True)
    assert (group / "imputed_prey.txt").exists()


def test_the_group_files_carry_no_header(group_inputs):
    """SAINTexpress reads them positionally."""
    src, group = group_inputs(["t1_1", "c_1"])
    ed = _FakeED(["t1_1"], ["c_1"])

    score._build_group_saint_inputs(str(src), str(group), ed, 1, False)
    first = (group / "bait.txt").read_text().splitlines()[0]

    assert first.split("\t")[0] == "t1_1"


def test_a_group_matching_nothing_is_reported(group_inputs, caplog):
    """Documented, not fixed: the group is still written, empty, and the run continues to
    SAINTexpress rather than stopping here."""
    src, group = group_inputs(["t1_1", "c_1"])
    ed = _FakeED(["absent"], [])

    with caplog.at_level("ERROR", logger="proximate.score"):
        score._build_group_saint_inputs(str(src), str(group), ed, 1, False)

    assert "zero rows" in caplog.text
    assert (group / "bait.txt").read_text() == ""


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
        fake_saint()
        record = _Record()
        score._score(args, record)
        return record

    return _run


def test_a_legacy_run_writes_the_merged_output(score_inputs, run_score):
    args = score_inputs()

    record = run_score(args)

    assert os.path.exists(os.path.join(args.outputPath, "merged.csv"))
    assert record.extras["run_mode"] == "legacy"


def test_the_merged_output_keeps_the_columns_annotation_reads(score_inputs, run_score):
    """Annotation joins on these names, so a rename here surfaces as missing annotations
    rather than as an error."""
    args = score_inputs()
    run_score(args)

    merged = pd.read_csv(os.path.join(args.outputPath, "merged.csv"))

    assert {"Experiment.ID", "Prey.ID", "Bait.ID", "source_group"} <= set(merged.columns)
    assert "Bait_x" not in merged.columns
    assert "Bait_y" not in merged.columns


def test_a_legacy_run_still_carries_a_source_group_column(score_inputs, run_score):
    """The schema is the same whether or not the run was grouped, so anything reading
    merged.csv does not have to branch on it."""
    args = score_inputs()
    run_score(args)

    merged = pd.read_csv(os.path.join(args.outputPath, "merged.csv"))

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


def test_comppass_results_are_written_alongside_the_merge(score_inputs, run_score):
    args = score_inputs()
    record = run_score(args)

    assert os.path.exists(os.path.join(args.outputPath, "compPASS.csv"))
    assert record.metrics["comppass_rows_out"] > 0


def test_the_filter_writes_the_interaction_file_saint_reads(score_inputs, run_score):
    args = score_inputs()
    run_score(args)

    assert os.path.exists(os.path.join(args.scoreInputs, "filtered_interaction.txt"))


@pytest.mark.parametrize("missing", ["prey.txt", "interaction.txt", "bait.txt",
                                     "to_CompPASS.csv"])
def test_every_required_input_is_checked(score_inputs, run_score, missing):
    args = score_inputs()
    os.remove(os.path.join(args.scoreInputs, missing))

    with pytest.raises(SystemExit) as excinfo:
        run_score(args)

    assert excinfo.value.code == 1


def test_a_missing_input_is_still_recorded(score_inputs, fake_saint):
    """The manifest names the file that was looked for, which is what makes a failed run
    diagnosable after the fact."""
    args = score_inputs()
    os.remove(os.path.join(args.scoreInputs, "bait.txt"))
    fake_saint()
    record = _Record()

    with pytest.raises(SystemExit):
        score._score(args, record)

    assert any(role == "bait.txt" for _, role in record.inputs)


def test_single_bait_imputation_requires_a_bait(score_inputs, run_score):
    args = score_inputs()
    args.imputation = "2"
    args.pi_method = "single_bait"

    with pytest.raises(SystemExit) as excinfo:
        run_score(args)

    assert excinfo.value.code == 1


def test_a_non_numeric_imputation_stops_the_run(score_inputs, run_score, caplog):
    """Documented, not fixed: the unrecognized value reaches int(), and the ValueError is
    caught by the imputation handler, so the log says imputation failed rather than that
    the option was invalid."""
    args = score_inputs()
    args.imputation = "yes"

    with caplog.at_level("ERROR", logger="proximate.score"):
        with pytest.raises(SystemExit):
            run_score(args)

    assert "Imputation failed" in caplog.text


def test_a_grouped_run_scores_each_group_separately(score_inputs, fake_saint):
    args = score_inputs(group_column=["1"] * 6 + ["*", "*"])
    calls = fake_saint()
    record = _Record()

    score._score(args, record)

    assert record.extras["run_mode"] == "grouped"
    assert os.path.isdir(os.path.join(args.outputPath, "groups", "1"))
    assert len(calls) == 1


def test_a_grouped_run_labels_rows_with_their_group(score_inputs, fake_saint):
    args = score_inputs(group_column=["1"] * 6 + ["*", "*"])
    fake_saint()

    score._score(args, _Record())
    merged = pd.read_csv(os.path.join(args.outputPath, "merged.csv"))

    assert set(merged["source_group"]) == {1}


def test_a_grouped_run_records_a_manifest(score_inputs, fake_saint):
    args = score_inputs(group_column=["1"] * 6 + ["*", "*"])
    fake_saint()

    score._score(args, _Record())

    assert os.path.exists(os.path.join(args.outputPath, "groups", "manifest.json"))
