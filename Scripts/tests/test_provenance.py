"""Tests for the run manifest written alongside each dataset's results, and for the
logging configuration that puts a ``proximate.log`` beside it."""

import hashlib
import json
import os

import pytest

import log_config
import provenance
import setup_datasets


def _manifest(output_dir):
    with open(os.path.join(output_dir, provenance.RUN_JSON_FILENAME), encoding="utf-8") as handle:
        return json.load(handle)


def _only_run(output_dir):
    runs = _manifest(output_dir)["runs"]
    assert len(runs) == 1
    return next(iter(runs.values()))


def _lines(path):
    with open(path, encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle if line.strip()]


@pytest.fixture
def run_id(monkeypatch):
    """Pin the run ID so stages in a test share one run without minting."""
    value = "20260831T120000Z-0badcafe"
    monkeypatch.setenv("PROXIMATE_RUN_ID", value)
    return value


# ---------------------------------------------------------------------------
# Stage lifecycle
# ---------------------------------------------------------------------------

def test_failing_stage_records_the_error_and_reraises(tmp_path, run_id):
    with pytest.raises(ValueError, match="boom"):
        with provenance.stage(tmp_path, "score"):
            raise ValueError("boom")

    stage = _only_run(tmp_path)["stages"][0]
    assert stage["status"] == "error"
    assert stage["error"]["type"] == "ValueError"
    assert "boom" in stage["error"]["message"]
    assert stage["error"]["traceback"]


@pytest.mark.parametrize("code, status", [(0, "ok"), (1, "error")])
def test_sys_exit_is_recorded_by_its_exit_code(tmp_path, run_id, code, status):
    """score.py and annotator.py signal failure with sys.exit(1); catching only
    Exception would file those runs as successful."""
    with pytest.raises(SystemExit):
        with provenance.stage(tmp_path, "score"):
            raise SystemExit(code)

    stage = _only_run(tmp_path)["stages"][0]
    assert stage["status"] == status
    assert stage["exit_code"] == code


def test_three_stages_accumulate_under_one_run(tmp_path, run_id):
    for name in ("parse", "score", "annotate"):
        with provenance.stage(tmp_path, name):
            pass

    run = _only_run(tmp_path)
    assert [s["stage"] for s in run["stages"]] == ["parse", "score", "annotate"]
    assert [s["status"] for s in run["stages"]] == ["ok"] * 3
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


# ---------------------------------------------------------------------------
# Recorded content
# ---------------------------------------------------------------------------

def test_stage_records_inputs_outputs_metrics_and_arguments(tmp_path, run_id):
    source = tmp_path / "ED.csv"
    source.write_bytes(b"Experiment Name,Type\ne1,T\n")

    with provenance.stage(tmp_path, "parse", cli_args={"path": tmp_path}) as record:
        record.add_input(source, role="experimentalDesign")
        record.add_input(tmp_path / "absent.txt", role="proteinGroups")
        record.add_output(source, rows=1)
        record.metric("n_experiments", 12)
        record.extra(quant_type="Intensity")

    stage = _only_run(tmp_path)["stages"][0]
    present, absent = stage["inputs"]
    assert present["role"] == "experimentalDesign"
    assert present["bytes"] == source.stat().st_size
    assert present["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert absent["role"] == "proteinGroups"
    assert absent["missing"] is True
    assert stage["outputs"][0]["rows"] == 1
    assert stage["metrics"]["n_experiments"] == 12
    assert stage["extra"]["quant_type"] == "Intensity"
    # A Path in the cli args must not abort the JSON write.
    assert isinstance(stage["cli_args"]["path"], str)


@pytest.mark.parametrize("version, label", [
    ({"version": "9f2c1ab", "source": "env"}, "version 9f2c1ab"),
    ({"version": "9f2c1ab", "source": "git"}, "version 9f2c1ab (source checkout)"),
    ({"version": None, "source": "unknown"}, "version unknown"),
])
def test_version_label(monkeypatch, version, label):
    """A commit read from a working tree may include uncommitted edits, so it must not
    be mistaken for the build that commit produced."""
    monkeypatch.setattr(provenance, "proximate_version", lambda: version)

    assert provenance.version_label() == label


# ---------------------------------------------------------------------------
# Robustness — provenance must never fail a scientific run
# ---------------------------------------------------------------------------

def test_a_corrupt_manifest_is_replaced_rather_than_raising(tmp_path, run_id):
    (tmp_path / provenance.RUN_JSON_FILENAME).write_text("{ this is not json")

    with provenance.stage(tmp_path, "parse"):
        pass

    assert _only_run(tmp_path)["stages"][0]["stage"] == "parse"


def test_a_failed_manifest_write_does_not_mask_the_real_work(tmp_path, run_id, monkeypatch):
    monkeypatch.setattr(provenance.json, "dump",
                        lambda *a, **k: (_ for _ in ()).throw(OSError("disk full")))

    with provenance.stage(tmp_path, "parse"):
        pass  # must not raise


# ---------------------------------------------------------------------------
# parse.py's stage decorator
# ---------------------------------------------------------------------------

def test_parse_stage_classifies_arguments(tmp_path, run_id):
    """A file argument is fingerprinted as an input, a scalar such as ``quantType``
    is a parameter, and a file that is not there is recorded as missing."""
    import parse

    source = tmp_path / "proteinGroups.txt"
    source.write_text("data\n")

    @parse._parse_stage
    def fake_entry(proteinGroups, experimentalDesign, quantType, outputPath):
        return 3, 1

    fake_entry(str(source), str(tmp_path / "absent.csv"), "LFQ", str(tmp_path))

    stage = _only_run(tmp_path)["stages"][0]
    by_role = {entry["role"]: entry for entry in stage["inputs"]}
    assert set(by_role) == {"proteinGroups", "experimentalDesign"}
    assert by_role["proteinGroups"]["sha256"]
    assert by_role["experimentalDesign"]["missing"] is True
    assert stage["params"]["quantType"] == "LFQ"
    assert stage["metrics"] == {"n_experiments": 3, "n_controls": 1}


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

    assert "inputs rejected" in (tmp_path / "proximate.log").read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# Locating a dataset's annotation databases
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("exclude_hcm, filename", [
    (False, setup_datasets.BIOGRID_SUMMARY_FILENAME),
    (True, setup_datasets.BIOGRID_NO_HCM_SUMMARY_FILENAME),
])
def test_the_biogrid_summary_is_looked_for_where_setup_datasets_writes_it(
        tmp_path, exclude_hcm, filename):
    """The summary is built per organism, so there is no copy at the top of the
    datasets directory to fall back on."""
    path = provenance.biogrid_summary_path("mouse", str(tmp_path), exclude_hcm=exclude_hcm)

    assert path == os.path.join(str(tmp_path), "mouse", filename)


def test_annotation_settings_come_from_the_latest_run_that_recorded_them(
        tmp_path, monkeypatch):
    """A dataset with no manifest, or a corrupt one, resolves to human against the
    full summary: every dataset scored before either became a parameter was.  Once
    recorded, the most recent annotating run wins."""
    assert provenance.annotation_settings(tmp_path) == {
        "organism": "human", "exclude_hcm": False}

    monkeypatch.setenv("PROXIMATE_RUN_ID", "20260831T120000Z-0badcafe")
    with provenance.stage(tmp_path, "annotate") as record:
        record.extra(organism="human", exclude_hcm=True)
    monkeypatch.setenv("PROXIMATE_RUN_ID", "20260901T120000Z-1badcafe")
    with provenance.stage(tmp_path, "annotate") as record:
        record.extra(organism="mouse")

    assert provenance.annotation_settings(tmp_path) == {"organism": "mouse", "exclude_hcm": True}

    (tmp_path / provenance.RUN_JSON_FILENAME).write_text("{not json", encoding="utf-8")
    assert provenance.annotation_settings(tmp_path) == {
        "organism": "human", "exclude_hcm": False}


# ---------------------------------------------------------------------------
# Logging: the per-dataset log and the run ID stamped on it
# ---------------------------------------------------------------------------

def test_add_file_handler_is_idempotent_for_the_same_file(tmp_path):
    """The long-lived GUI process calls this per action; a second call, or one
    through an unnormalized spelling of the path, must not duplicate every
    subsequent line."""
    target = tmp_path / "proximate.log"

    log_config.add_file_handler(target)
    log_config.add_file_handler(target)
    log_config.add_file_handler(tmp_path / "." / "proximate.log")
    log_config.get_logger("score").info("only once")

    assert sum("only once" in line for line in _lines(target)) == 1


def test_nested_dataset_log_blocks_share_one_handler_until_the_outer_exits(tmp_path):
    """app.py wraps an action in ``dataset_log`` and the parse stage opens another on
    the same directory inside it.  The inner exit must not detach the handler the
    outer block is still using, and once the outer block has exited (here by
    raising) nothing further may reach the file."""
    logger = log_config.get_logger("app")
    target = tmp_path / "proximate.log"

    with pytest.raises(ValueError):
        with log_config.dataset_log(tmp_path):
            with log_config.dataset_log(tmp_path):
                logger.info("inside both")
            logger.info("after inner")
            raise ValueError("boom")
    logger.info("after outer")

    lines = _lines(target)
    assert sum("inside both" in line for line in lines) == 1
    assert any("after inner" in line for line in lines)
    assert not any("after outer" in line for line in lines)
    assert not log_config._file_handlers


def test_run_id_is_minted_exported_and_overridden_by_run_context(tmp_path):
    """Children inherit os.environ, so exporting the minted ID is what shares it;
    ``run_context`` overrides it for a block, including one that raises."""
    target = tmp_path / "proximate.log"
    log_config.add_file_handler(target)
    logger = log_config.get_logger("app")

    ambient = log_config.get_run_id()
    assert os.environ["PROXIMATE_RUN_ID"] == ambient

    with pytest.raises(ValueError):
        with log_config.run_context("20260831T000000Z-scoped01"):
            logger.info("during")
            raise ValueError("boom")
    logger.info("after")

    assert log_config.get_run_id() == ambient
    lines = _lines(target)
    assert "20260831T000000Z-scoped01" in next(l for l in lines if "during" in l)
    assert ambient in next(l for l in lines if "after" in l)
