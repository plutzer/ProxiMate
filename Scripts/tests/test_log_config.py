"""Tests for the centralized logging configuration.

These tests all take the ``clean_logging`` fixture, which isolates the
process-wide ProxiMate logger and resets ``log_config``'s one-shot setup guard.
"""

import logging
import logging.handlers
import os
import re
import sys

import pytest


def _messages(path):
    """Read a log file, returning its lines with trailing newlines stripped."""
    with open(path, encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle if line.strip()]


# ---------------------------------------------------------------------------
# Logger scoping
# ---------------------------------------------------------------------------

def test_get_logger_returns_a_child_of_the_package_logger(clean_logging):
    logger = clean_logging.get_logger("score")
    assert logger.name == f"{clean_logging.PACKAGE}.score"


def test_a_script_run_directly_is_named_after_its_file(clean_logging, monkeypatch):
    """`python Scripts/parse.py` gives __name__ == "__main__"; without this the
    three CLI stages would all log under the same uninformative name."""
    monkeypatch.setattr(sys, "argv", ["/Scripts/parse.py", "--outputPath", "x"])

    assert clean_logging.get_logger("__main__").name == f"{clean_logging.PACKAGE}.parse"


def test_package_logger_does_not_propagate_to_root(clean_logging):
    clean_logging.setup_logging()
    assert logging.getLogger(clean_logging.PACKAGE).propagate is False


def test_setup_does_not_touch_the_root_logger(clean_logging):
    """Handlers belong to the ProxiMate logger so third-party output is untouched."""
    root = logging.getLogger()
    before = list(root.handlers)

    clean_logging.setup_logging()

    assert root.handlers == before


def test_debug_level_does_not_enable_third_party_loggers(clean_logging, monkeypatch):
    """The whole point of package scoping: LOG_LEVEL=DEBUG must not flood with
    matplotlib/urllib3 chatter."""
    monkeypatch.setenv("LOG_LEVEL", "DEBUG")

    clean_logging.setup_logging()

    assert logging.getLogger(clean_logging.PACKAGE).level == logging.DEBUG
    assert not logging.getLogger("matplotlib").isEnabledFor(logging.DEBUG)


# ---------------------------------------------------------------------------
# Level handling
# ---------------------------------------------------------------------------

def test_log_level_env_var_is_honored(clean_logging, monkeypatch):
    monkeypatch.setenv("LOG_LEVEL", "WARNING")

    clean_logging.setup_logging()

    assert logging.getLogger(clean_logging.PACKAGE).level == logging.WARNING


def test_invalid_log_level_warns_and_falls_back_to_info(clean_logging, monkeypatch, tmp_path):
    """A typo in LOG_LEVEL must be reported, not silently swallowed."""
    monkeypatch.setenv("LOG_LEVEL", "VERBOSE")
    target = tmp_path / "proximate.log"

    clean_logging.setup_logging()
    clean_logging.add_file_handler(target)

    assert logging.getLogger(clean_logging.PACKAGE).level == logging.INFO
    assert any("VERBOSE" in line for line in _messages(target))


# ---------------------------------------------------------------------------
# File handlers
# ---------------------------------------------------------------------------

def test_add_file_handler_writes_formatted_records(clean_logging, tmp_path):
    target = tmp_path / "proximate.log"
    clean_logging.add_file_handler(target)

    clean_logging.get_logger("score").info("SAINTexpress completed")

    lines = _messages(target)
    assert any("SAINTexpress completed" in line for line in lines)
    assert any("proximate.score" in line for line in lines)


def test_log_records_carry_the_run_id_and_pid(clean_logging, tmp_path, monkeypatch):
    """Three processes append to one dataset log; lines must be attributable."""
    monkeypatch.setenv("PROXIMATE_RUN_ID", "20260831T000000Z-deadbeef")
    target = tmp_path / "proximate.log"
    clean_logging.add_file_handler(target)

    clean_logging.get_logger("parse").info("parsing")

    line = next(l for l in _messages(target) if "parsing" in l)
    assert "20260831T000000Z-deadbeef" in line
    assert f"pid={os.getpid()}" in line


def test_add_file_handler_is_idempotent_for_the_same_path(clean_logging, tmp_path):
    """The long-lived GUI process calls this per action; a second call must not
    duplicate every subsequent line."""
    target = tmp_path / "proximate.log"

    clean_logging.add_file_handler(target)
    clean_logging.add_file_handler(target)
    clean_logging.get_logger("score").info("only once")

    assert sum("only once" in line for line in _messages(target)) == 1


def test_add_file_handler_normalizes_paths_before_deduplicating(clean_logging, tmp_path):
    """``dir/proximate.log`` and ``dir/./proximate.log`` are the same file."""
    target = tmp_path / "proximate.log"

    clean_logging.add_file_handler(target)
    clean_logging.add_file_handler(tmp_path / "." / "proximate.log")
    clean_logging.get_logger("score").info("only once")

    assert sum("only once" in line for line in _messages(target)) == 1


def test_remove_file_handler_stops_writing(clean_logging, tmp_path):
    target = tmp_path / "proximate.log"
    clean_logging.add_file_handler(target)

    clean_logging.remove_file_handler(target)
    clean_logging.get_logger("score").info("after removal")

    assert not any("after removal" in line for line in _messages(target))


# ---------------------------------------------------------------------------
# dataset_log context manager
# ---------------------------------------------------------------------------

def test_dataset_log_captures_inside_the_block_only(clean_logging, tmp_path):
    logger = clean_logging.get_logger("app")

    with clean_logging.dataset_log(tmp_path):
        logger.info("inside")
    logger.info("outside")

    lines = _messages(tmp_path / "proximate.log")
    assert any("inside" in line for line in lines)
    assert not any("outside" in line for line in lines)


def test_dataset_log_detaches_when_the_block_raises(clean_logging, tmp_path):
    """A failed scoring run must not leave the GUI writing to that dataset forever."""
    logger = clean_logging.get_logger("app")

    with pytest.raises(ValueError):
        with clean_logging.dataset_log(tmp_path):
            raise ValueError("boom")
    logger.info("after the failure")

    assert not any("after the failure" in line
                   for line in _messages(tmp_path / "proximate.log"))


def test_dataset_log_creates_the_output_directory(clean_logging, tmp_path):
    """annotator.py attaches its log before creating its output directory."""
    target = tmp_path / "not_yet_created"

    with clean_logging.dataset_log(target):
        clean_logging.get_logger("annotator").info("annotating")

    assert (target / "proximate.log").exists()


# ---------------------------------------------------------------------------
# Run IDs
# ---------------------------------------------------------------------------

def test_run_id_is_taken_from_the_environment_when_present(clean_logging, monkeypatch):
    monkeypatch.setenv("PROXIMATE_RUN_ID", "20260831T000000Z-abcd1234")

    assert clean_logging.get_run_id() == "20260831T000000Z-abcd1234"


def test_run_id_is_minted_and_exported_when_absent(clean_logging):
    """Children inherit os.environ, so exporting it is what shares the ID."""
    run_id = clean_logging.get_run_id()

    assert os.environ["PROXIMATE_RUN_ID"] == run_id
    assert clean_logging.get_run_id() == run_id


def test_run_context_overrides_the_ambient_run_id(clean_logging, monkeypatch):
    monkeypatch.setenv("PROXIMATE_RUN_ID", "20260831T000000Z-ambient0")

    with clean_logging.run_context("20260831T000000Z-scoped01"):
        assert clean_logging.get_run_id() == "20260831T000000Z-scoped01"


def test_run_context_restores_the_previous_run_id(clean_logging, monkeypatch):
    monkeypatch.setenv("PROXIMATE_RUN_ID", "20260831T000000Z-ambient0")

    with clean_logging.run_context("20260831T000000Z-scoped01"):
        pass

    assert clean_logging.get_run_id() == "20260831T000000Z-ambient0"


def test_run_context_restores_the_run_id_when_the_block_raises(clean_logging, monkeypatch):
    """A failed scoring run must not leave later lines stamped with its ID."""
    monkeypatch.setenv("PROXIMATE_RUN_ID", "20260831T000000Z-ambient0")

    with pytest.raises(ValueError):
        with clean_logging.run_context("20260831T000000Z-scoped01"):
            raise ValueError("boom")

    assert clean_logging.get_run_id() == "20260831T000000Z-ambient0"


def test_run_context_stamps_records_written_inside_it(clean_logging, tmp_path):
    target = tmp_path / "proximate.log"
    clean_logging.add_file_handler(target)
    logger = clean_logging.get_logger("app")

    with clean_logging.run_context("20260831T000000Z-scoped01"):
        logger.info("during")
    logger.info("after")

    lines = _messages(target)
    during = next(l for l in lines if "during" in l)
    after = next(l for l in lines if "after" in l)
    assert "20260831T000000Z-scoped01" in during
    assert "20260831T000000Z-scoped01" not in after


def test_new_run_id_has_the_documented_shape(clean_logging):
    assert re.fullmatch(r"\d{8}T\d{6}Z-[0-9a-f]{8}", clean_logging.new_run_id())


def test_new_run_ids_are_unique(clean_logging):
    ids = {clean_logging.new_run_id() for _ in range(50)}

    assert len(ids) == 50


# ---------------------------------------------------------------------------
# Unified operational log
# ---------------------------------------------------------------------------

def test_unified_log_is_written_to_the_configured_directory(clean_logging, log_dir):
    clean_logging.setup_logging()

    clean_logging.get_logger("app").info("server started")

    assert any("server started" in line
               for line in _messages(log_dir / "proximate-server.log"))


def test_unified_log_handler_rotates(clean_logging, log_dir):
    """The server log is unbounded in time, so it must not grow without limit."""
    clean_logging.setup_logging()

    handlers = [h for h in logging.getLogger(clean_logging.PACKAGE).handlers
                if isinstance(h, logging.handlers.RotatingFileHandler)]

    assert handlers, "expected a RotatingFileHandler for the unified log"
    assert handlers[0].maxBytes > 0
    assert handlers[0].backupCount > 0


def test_unusable_explicit_log_dir_warns_without_raising(clean_logging, tmp_path, monkeypatch):
    """An explicitly requested destination that cannot be used is an error worth
    reporting, but must not stop the app from starting."""
    blocker = tmp_path / "not_a_directory"
    blocker.write_text("")
    monkeypatch.setenv("PROXIMATE_LOG_DIR", str(blocker))
    dataset_log = tmp_path / "proximate.log"

    clean_logging.setup_logging()
    clean_logging.add_file_handler(dataset_log)
    clean_logging.get_logger("app").info("still running")

    lines = _messages(dataset_log)
    assert any("still running" in line for line in lines)
    assert any(str(blocker) in line and "WARNING" in line for line in lines)


def test_missing_default_log_dir_is_not_a_warning(clean_logging, tmp_path, monkeypatch):
    """Running from a source checkout, /Outputs does not exist.  That is normal
    and must not warn on every start."""
    monkeypatch.delenv("PROXIMATE_LOG_DIR", raising=False)
    monkeypatch.setattr(clean_logging, "DEFAULT_LOG_DIR", str(tmp_path / "absent"))
    dataset_log = tmp_path / "proximate.log"

    clean_logging.setup_logging()
    clean_logging.add_file_handler(dataset_log)
    clean_logging.get_logger("app").info("started")

    assert not any("WARNING" in line for line in _messages(dataset_log))


def test_log_file_env_var_still_works(clean_logging, tmp_path, monkeypatch):
    """Backward compatibility with the documented LOG_FILE knob."""
    target = tmp_path / "explicit.log"
    monkeypatch.setenv("LOG_FILE", str(target))

    clean_logging.setup_logging()
    clean_logging.get_logger("score").info("via LOG_FILE")

    assert any("via LOG_FILE" in line for line in _messages(target))
