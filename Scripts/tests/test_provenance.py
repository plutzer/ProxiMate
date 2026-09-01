"""Tests for the run manifest written alongside each dataset's results."""

import hashlib
import json
import os

import pytest

import provenance
import setup_datasets


def _manifest(output_dir):
    with open(os.path.join(output_dir, provenance.RUN_JSON_FILENAME), encoding="utf-8") as handle:
        return json.load(handle)


def _only_run(output_dir):
    runs = _manifest(output_dir)["runs"]
    assert len(runs) == 1
    return next(iter(runs.values()))


@pytest.fixture
def run_id(monkeypatch):
    """Pin the run ID so stages in a test share one run without minting."""
    value = "20260831T120000Z-0badcafe"
    monkeypatch.setenv("PROXIMATE_RUN_ID", value)
    return value


# ---------------------------------------------------------------------------
# Stage lifecycle
# ---------------------------------------------------------------------------

def test_stage_writes_a_manifest_with_the_schema_version(tmp_path, run_id):
    with provenance.stage(tmp_path, "parse"):
        pass

    assert _manifest(tmp_path)["schema_version"] == provenance.SCHEMA_VERSION


def test_successful_stage_is_recorded_as_ok(tmp_path, run_id):
    with provenance.stage(tmp_path, "parse"):
        pass

    stage = _only_run(tmp_path)["stages"][0]
    assert stage["stage"] == "parse"
    assert stage["status"] == "ok"
    assert stage["wall_seconds"] >= 0


def test_failing_stage_records_the_error_and_reraises(tmp_path, run_id):
    with pytest.raises(ValueError, match="boom"):
        with provenance.stage(tmp_path, "score"):
            raise ValueError("boom")

    stage = _only_run(tmp_path)["stages"][0]
    assert stage["status"] == "error"
    assert stage["error"]["type"] == "ValueError"
    assert "boom" in stage["error"]["message"]
    assert stage["error"]["traceback"]


def test_sys_exit_is_recorded_as_a_failure(tmp_path, run_id):
    """score.py and annotator.py signal failure with sys.exit(1); catching only
    Exception would file those runs as successful."""
    with pytest.raises(SystemExit):
        with provenance.stage(tmp_path, "score"):
            raise SystemExit(1)

    stage = _only_run(tmp_path)["stages"][0]
    assert stage["status"] == "error"
    assert stage["exit_code"] == 1


def test_sys_exit_zero_is_recorded_as_success(tmp_path, run_id):
    with pytest.raises(SystemExit):
        with provenance.stage(tmp_path, "parse"):
            raise SystemExit(0)

    assert _only_run(tmp_path)["stages"][0]["status"] == "ok"


def test_stage_creates_the_output_directory(tmp_path, run_id):
    target = tmp_path / "not_yet_created"

    with provenance.stage(target, "parse"):
        pass

    assert (target / provenance.RUN_JSON_FILENAME).exists()


# ---------------------------------------------------------------------------
# Accumulating stages and runs
# ---------------------------------------------------------------------------

def test_three_stages_accumulate_under_one_run(tmp_path, run_id):
    for name in ("parse", "score", "annotate"):
        with provenance.stage(tmp_path, name):
            pass

    run = _only_run(tmp_path)
    assert [s["stage"] for s in run["stages"]] == ["parse", "score", "annotate"]
    assert run["run_id"] == run_id


def test_rescoring_a_dataset_keeps_the_earlier_run(tmp_path, monkeypatch):
    """Re-running must not erase the record of how the previous result was made."""
    monkeypatch.setenv("PROXIMATE_RUN_ID", "20260831T120000Z-aaaaaaaa")
    with provenance.stage(tmp_path, "score"):
        pass
    monkeypatch.setenv("PROXIMATE_RUN_ID", "20260831T130000Z-bbbbbbbb")
    with provenance.stage(tmp_path, "score"):
        pass

    assert set(_manifest(tmp_path)["runs"]) == {
        "20260831T120000Z-aaaaaaaa", "20260831T130000Z-bbbbbbbb"}


def test_stage_mints_a_run_id_when_none_is_set(tmp_path, clean_logging):
    with provenance.stage(tmp_path, "parse"):
        pass

    assert _only_run(tmp_path)["run_id"] == os.environ["PROXIMATE_RUN_ID"]


# ---------------------------------------------------------------------------
# Recorded content
# ---------------------------------------------------------------------------

def test_stage_records_inputs_outputs_and_metrics(tmp_path, run_id):
    source = tmp_path / "ED.csv"
    source.write_text("Experiment Name,Type\ne1,T\n")

    with provenance.stage(tmp_path, "parse") as record:
        record.add_input(source, role="experimentalDesign")
        record.add_output(source, rows=1)
        record.metric("n_experiments", 12)
        record.extra(quant_type="Intensity")

    stage = _only_run(tmp_path)["stages"][0]
    assert stage["inputs"][0]["role"] == "experimentalDesign"
    assert stage["inputs"][0]["bytes"] == source.stat().st_size
    assert stage["outputs"][0]["rows"] == 1
    assert stage["metrics"]["n_experiments"] == 12
    assert stage["extra"]["quant_type"] == "Intensity"


def test_recorded_inputs_carry_a_sha256(tmp_path, run_id):
    source = tmp_path / "input.txt"
    source.write_bytes(b"proteinGroups")

    with provenance.stage(tmp_path, "parse") as record:
        record.add_input(source)

    expected = hashlib.sha256(b"proteinGroups").hexdigest()
    assert _only_run(tmp_path)["stages"][0]["inputs"][0]["sha256"] == expected


def test_a_missing_input_is_recorded_rather_than_raising(tmp_path, run_id):
    with provenance.stage(tmp_path, "parse") as record:
        record.add_input(tmp_path / "absent.txt", role="proteinGroups")

    entry = _only_run(tmp_path)["stages"][0]["inputs"][0]
    assert entry["role"] == "proteinGroups"
    assert entry["missing"] is True


def test_cli_args_are_recorded(tmp_path, run_id):
    import argparse
    args = argparse.Namespace(quantType="LFQ", n_iterations=5, seed=1234)

    with provenance.stage(tmp_path, "score", cli_args=vars(args)):
        pass

    recorded = _only_run(tmp_path)["stages"][0]["cli_args"]
    assert recorded == {"quantType": "LFQ", "n_iterations": 5, "seed": 1234}


def test_cli_args_are_copied_not_captured_by_reference(tmp_path, run_id):
    """`vars(args)` hands over argparse's live namespace dict."""
    live = {"quantType": "LFQ"}

    with provenance.stage(tmp_path, "score", cli_args=live):
        live["quantType"] = "MUTATED"

    assert _only_run(tmp_path)["stages"][0]["cli_args"]["quantType"] == "LFQ"


def test_unserializable_values_do_not_break_the_manifest(tmp_path, run_id):
    """argparse namespaces and metrics can hold numpy scalars and Paths."""
    with provenance.stage(tmp_path, "score", cli_args={"path": tmp_path}) as record:
        record.metric("rows", 10)

    assert isinstance(_only_run(tmp_path)["stages"][0]["cli_args"]["path"], str)


def test_environment_is_recorded_once_per_run(tmp_path, run_id):
    with provenance.stage(tmp_path, "parse"):
        pass

    environment = _only_run(tmp_path)["environment"]
    assert environment["python"]
    assert environment["packages"]["pandas"]


def test_version_prefers_the_baked_in_environment_variable(tmp_path, monkeypatch):
    monkeypatch.setenv("PROXIMATE_VERSION", "9f2c1ab")

    version = provenance.proximate_version()

    assert version["version"] == "9f2c1ab"
    assert version["source"] == "env"


def test_a_baked_version_is_labelled_bare(monkeypatch):
    """What a built image shows: the stamp and nothing else."""
    monkeypatch.setattr(provenance, "proximate_version",
                        lambda: {"version": "9f2c1ab", "source": "env"})

    assert provenance.version_label() == "version 9f2c1ab"


def test_a_source_checkout_is_labelled_as_one(monkeypatch):
    """A commit read from a working tree may include uncommitted edits, so it must not
    be mistaken for the build that commit produced."""
    monkeypatch.setattr(provenance, "proximate_version",
                        lambda: {"version": "9f2c1ab", "source": "git"})

    assert provenance.version_label() == "version 9f2c1ab (source checkout)"


def test_an_unidentifiable_build_is_labelled_unknown(monkeypatch):
    monkeypatch.setattr(provenance, "proximate_version",
                        lambda: {"version": None, "source": "unknown"})

    assert provenance.version_label() == "version unknown"


def test_an_unstamped_image_is_also_labelled_unknown(monkeypatch):
    """The Dockerfile defaults PROXIMATE_VERSION to the literal "unknown", so an image
    built without the build argument reports that string rather than nothing.  It says
    as little as an absent version and must read the same way."""
    monkeypatch.setenv("PROXIMATE_VERSION", "unknown")

    assert provenance.version_label() == "version unknown"


# ---------------------------------------------------------------------------
# Robustness — provenance must never fail a scientific run
# ---------------------------------------------------------------------------

def test_a_corrupt_manifest_is_replaced_rather_than_raising(tmp_path, run_id):
    target = tmp_path / provenance.RUN_JSON_FILENAME
    target.write_text("{ this is not json")

    with provenance.stage(tmp_path, "parse"):
        pass

    assert _only_run(tmp_path)["stages"][0]["stage"] == "parse"


def test_a_failed_manifest_write_does_not_mask_the_real_work(tmp_path, run_id, monkeypatch):
    """A provenance problem must never turn a successful run into a failure."""
    monkeypatch.setattr(provenance.json, "dump",
                        lambda *a, **k: (_ for _ in ()).throw(OSError("disk full")))

    with provenance.stage(tmp_path, "parse"):
        pass  # must not raise


def test_the_manifest_write_leaves_no_temporary_files(tmp_path, run_id):
    with provenance.stage(tmp_path, "parse"):
        pass

    assert [p.name for p in tmp_path.iterdir()] == [provenance.RUN_JSON_FILENAME]


# ---------------------------------------------------------------------------
# Dataset build info
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# parse.py's stage decorator
# ---------------------------------------------------------------------------

def test_scalar_arguments_are_not_recorded_as_input_files(tmp_path, run_id):
    """`quantType` is a value like "LFQ", not a path.  Fingerprinting it files a
    bogus missing-input entry beside the real ones."""
    import parse

    source = tmp_path / "proteinGroups.txt"
    source.write_text("data\n")

    @parse._parse_stage
    def fake_entry(proteinGroups, quantType, outputPath):
        return 3, 1

    fake_entry(str(source), "LFQ", str(tmp_path))

    inputs = _only_run(tmp_path)["stages"][0]["inputs"]
    assert [i["role"] for i in inputs] == ["proteinGroups"]
    assert inputs[0]["sha256"]


def test_a_parse_that_fails_validation_still_writes_a_dataset_log(tmp_path, run_id):
    """Validation runs before each entry point attaches its own log handler, so
    without the stage attaching one first, the failure people most need to read
    about leaves no proximate.log at all."""
    import parse

    @parse._parse_stage
    def fake_entry(proteinGroups, outputPath):
        parse.logger.error("inputs rejected: no such file")
        raise ValueError("bad input")

    with pytest.raises(ValueError):
        fake_entry(str(tmp_path / "absent.txt"), str(tmp_path))

    log = tmp_path / "proximate.log"
    assert log.exists()
    assert "inputs rejected" in log.read_text(encoding="utf-8")


def test_scalar_arguments_are_still_recorded_as_parameters(tmp_path, run_id):
    import parse

    @parse._parse_stage
    def fake_entry(quantType, outputPath):
        return 1, 0

    fake_entry("Spectral Counts", str(tmp_path))

    assert _only_run(tmp_path)["stages"][0]["params"]["quantType"] == "Spectral Counts"


def test_a_missing_input_file_is_still_recorded_as_missing(tmp_path, run_id):
    """Distinguishing a scalar from a path must not silence the genuine
    "the input was not there" signal."""
    import parse

    @parse._parse_stage
    def fake_entry(proteinGroups, outputPath):
        return 1, 0

    fake_entry(str(tmp_path / "absent.txt"), str(tmp_path))

    inputs = _only_run(tmp_path)["stages"][0]["inputs"]
    assert inputs[0]["role"] == "proteinGroups"
    assert inputs[0]["missing"] is True


def test_build_info_is_read_from_the_file_setup_datasets_writes(tmp_path):
    setup_datasets.write_build_info(str(tmp_path), {"biogrid": True})

    assert "Build date" in provenance.read_build_info(str(tmp_path))


def test_missing_build_info_reads_as_none(tmp_path):
    assert provenance.read_build_info(str(tmp_path)) is None


# ---------------------------------------------------------------------------
# Locating a dataset's annotation databases
# ---------------------------------------------------------------------------

def test_the_biogrid_summary_is_looked_for_where_setup_datasets_writes_it(tmp_path):
    """The summary is built per organism, so there is no copy at the top of the
    datasets directory to fall back on."""
    path = provenance.biogrid_summary_path("mouse", str(tmp_path))

    assert path == os.path.join(str(tmp_path), "mouse",
                                setup_datasets.BIOGRID_SUMMARY_FILENAME)


def test_the_organism_is_read_from_the_run_that_annotated_the_dataset(tmp_path, run_id):
    with provenance.stage(tmp_path, "annotate") as record:
        record.extra(organism="yeast")

    assert provenance.dataset_organism(tmp_path) == "yeast"


def test_a_dataset_with_no_manifest_falls_back_to_human(tmp_path):
    """Datasets scored before the organism was recorded still have to resolve to
    something, and every one of them was human."""
    assert provenance.dataset_organism(tmp_path) == "human"


def test_a_manifest_recording_no_organism_falls_back_to_human(tmp_path, run_id):
    with provenance.stage(tmp_path, "parse"):
        pass

    assert provenance.dataset_organism(tmp_path) == "human"


def test_the_most_recently_recorded_organism_wins(tmp_path, monkeypatch):
    """A dataset re-scored against another organism is annotated against that one."""
    monkeypatch.setenv("PROXIMATE_RUN_ID", "20260831T120000Z-0badcafe")
    with provenance.stage(tmp_path, "annotate") as record:
        record.extra(organism="human")
    monkeypatch.setenv("PROXIMATE_RUN_ID", "20260901T120000Z-1badcafe")
    with provenance.stage(tmp_path, "annotate") as record:
        record.extra(organism="mouse")

    assert provenance.dataset_organism(tmp_path) == "mouse"


def test_an_unreadable_manifest_falls_back_rather_than_raising(tmp_path):
    """QC plots resolve this on every redraw; a corrupt manifest must not take the
    panel down with it."""
    (tmp_path / provenance.RUN_JSON_FILENAME).write_text("{not json", encoding="utf-8")

    assert provenance.dataset_organism(tmp_path) == "human"
