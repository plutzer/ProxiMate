"""
Run manifests: a structured record of how a dataset's results were produced.

Each output directory gets a ``run.json`` recording, for every pipeline stage
that touched it, the parameters it ran with, the identity of its inputs, the
shapes of its outputs, and whether it succeeded.  Together with the
``proximate.log`` beside it, that is enough to answer "how was this made?"
months later, and it travels with the results.

Usage from a pipeline stage::

    with provenance.stage(args.outputPath, "score", cli_args=vars(args)) as record:
        record.add_input(interaction_path, role="interaction")
        ...
        record.metric("merged_rows", len(merged))
        record.add_output(merged_path, rows=len(merged))

Recording is best-effort by design: a manifest that cannot be written produces a
warning, never an exception.  Provenance must not be able to fail a scientific
run that otherwise succeeded.
"""

import copy
import getpass
import hashlib
import json
import os
import platform
import socket
import subprocess
import sys
import tempfile
import time
import traceback
from datetime import datetime, timezone

from log_config import get_logger, get_run_id

logger = get_logger(__name__)

SCHEMA_VERSION = 1
RUN_JSON_FILENAME = "run.json"

DEFAULT_DATASETS_DIR = "/Datasets"

# Recorded for every run so a result can be reproduced against the same stack.
TRACKED_PACKAGES = ("pandas", "numpy", "scipy", "statsmodels", "scikit-learn",
                    "matplotlib", "shiny", "plotly")

# Hashing a very large input costs more than the provenance is worth.
MAX_HASH_BYTES = 2 * 1024 * 1024 * 1024

# Only the tail of a traceback is worth keeping; the rest is noise in a manifest.
MAX_TRACEBACK_CHARS = 4000


def _utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _jsonable(value):
    """Coerce a value into something ``json.dump`` accepts.

    argparse namespaces and metrics carry Paths and numpy scalars, which would
    otherwise abort the write and lose the whole manifest.
    """
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(v) for v in value]
    if hasattr(value, "item"):  # numpy scalar
        try:
            return _jsonable(value.item())
        except (ValueError, AttributeError):
            pass
    return str(value)


def proximate_version():
    """Identify the running code.

    ``PROXIMATE_VERSION`` is baked into the image at build time; the container
    has no ``.git``, so the git fallback only ever applies to a source checkout.
    """
    baked = os.environ.get("PROXIMATE_VERSION")
    if baked:
        return {"version": baked, "source": "env"}

    try:
        sha = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=os.path.dirname(os.path.abspath(__file__)),
            capture_output=True, text=True, timeout=5)
        if sha.returncode == 0 and sha.stdout.strip():
            return {"version": sha.stdout.strip(), "source": "git"}
    except (OSError, subprocess.SubprocessError):
        pass
    return {"version": None, "source": "unknown"}


def version_label():
    """Describe the running build in one line, for display to a user.

    A commit read from a working tree may carry uncommitted edits, so it is marked as a
    checkout rather than presented as the build that commit produced.  An image built
    without ``--build-arg PROXIMATE_VERSION`` carries the literal "unknown" the
    Dockerfile defaults to, which identifies the build no better than an absent version
    and so reads the same way.
    """
    version = proximate_version()
    if not version["version"]:
        return "version unknown"
    if version["source"] == "git":
        return f"version {version['version']} (source checkout)"
    return f"version {version['version']}"


def _package_versions():
    from importlib import metadata
    versions = {}
    for name in TRACKED_PACKAGES:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    return versions


def environment_snapshot():
    """Describe the interpreter and stack this run executed on."""
    try:
        user = getpass.getuser()
    except (KeyError, OSError):
        # No password-file entry for the container's UID.
        user = None
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "hostname": socket.gethostname(),
        "user": user,
        "executable": sys.executable,
        "argv": [str(a) for a in sys.argv],
        "packages": _package_versions(),
    }


def file_record(path, role=None, rows=None):
    """Describe a file by size, digest and modification time.

    A missing file is recorded as missing rather than omitted: "the input was not
    there" is itself worth knowing when reading a manifest later.
    """
    entry = {"role": role, "path": str(path), "name": os.path.basename(str(path))}
    if rows is not None:
        entry["rows"] = rows

    try:
        size = os.path.getsize(path)
    except OSError:
        entry["missing"] = True
        return entry

    entry["bytes"] = size
    entry["mtime_utc"] = datetime.fromtimestamp(
        os.path.getmtime(path), timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    if size > MAX_HASH_BYTES:
        entry["sha256"] = None
        entry["hash_skipped"] = "too_large"
        return entry

    try:
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        entry["sha256"] = digest.hexdigest()
    except OSError as error:
        entry["sha256"] = None
        entry["hash_skipped"] = str(error)
    return entry


def read_build_info(datasets_dir=DEFAULT_DATASETS_DIR):
    """Return the annotation datasets' build stamp, or None if absent.

    The filename comes from ``setup_datasets`` so the two stay in step.
    """
    from setup_datasets import BUILD_INFO_FILENAME

    try:
        with open(os.path.join(datasets_dir, BUILD_INFO_FILENAME), encoding="utf-8") as handle:
            return handle.read()
    except OSError:
        return None


def biogrid_summary_path(organism, datasets_dir=None, exclude_hcm=False):
    """Where ``setup_datasets`` writes the BioGRID summary for one organism.

    The summary is filtered by taxonomy id, so each organism gets its own; there is no
    copy at the top of the datasets directory.  ``exclude_hcm`` selects the variant with
    Human Cell Map evidence removed.  The filenames come from ``setup_datasets`` so the
    two stay in step.
    """
    from setup_datasets import BIOGRID_NO_HCM_SUMMARY_FILENAME, BIOGRID_SUMMARY_FILENAME

    if datasets_dir is None:
        datasets_dir = DEFAULT_DATASETS_DIR
    filename = BIOGRID_NO_HCM_SUMMARY_FILENAME if exclude_hcm else BIOGRID_SUMMARY_FILENAME
    return os.path.join(datasets_dir, organism, filename)


def dataset_organism(output_dir, default="human"):
    """The organism a dataset's results were annotated against.

    Falls back to ``default`` for a dataset whose manifest records none, which covers
    every run scored before the organism became a parameter.
    """
    return _annotation_setting(output_dir, "organism", default)


def dataset_excludes_hcm(output_dir):
    """Whether a dataset was annotated against the BioGRID summary without Human Cell
    Map evidence.  A manifest that records nothing means the full summary was used."""
    return _annotation_setting(output_dir, "exclude_hcm", False)


def _annotation_setting(output_dir, key, default):
    """The last value of ``key`` recorded in any stage's ``extra`` of the manifest.

    A manifest that cannot be read yields ``default``: callers resolve these to draw a
    plot, and an unreadable manifest is not a reason to fail one.
    """
    document = _load(os.path.join(str(output_dir), RUN_JSON_FILENAME))

    value = default
    for run in document["runs"].values():
        for entry in run.get("stages", []):
            extra = entry.get("extra", {})
            if key in extra:
                value = extra[key]
    return value


def _load(manifest_path):
    """Read an existing manifest, starting fresh if it is unusable."""
    try:
        with open(manifest_path, encoding="utf-8") as handle:
            document = json.load(handle)
    except FileNotFoundError:
        return {"schema_version": SCHEMA_VERSION, "runs": {}}
    except (OSError, ValueError) as error:
        logger.warning("Could not read %s (%s); starting a new manifest.",
                       manifest_path, error)
        return {"schema_version": SCHEMA_VERSION, "runs": {}}

    document.setdefault("schema_version", SCHEMA_VERSION)
    document.setdefault("runs", {})
    return document


def _write(manifest_path, document):
    """Replace the manifest atomically, so a reader never sees a partial file."""
    directory = os.path.dirname(manifest_path) or "."
    handle = tempfile.NamedTemporaryFile(
        "w", dir=directory, prefix=".run.json.", delete=False, encoding="utf-8")
    try:
        with handle:
            json.dump(document, handle, indent=2)
            handle.write("\n")
        os.replace(handle.name, manifest_path)
    except BaseException:
        try:
            os.unlink(handle.name)
        except OSError:
            pass
        raise


class StageRecord:
    """Accumulates what one pipeline stage did. Created by :func:`stage`."""

    def __init__(self, stage_name, entrypoint=None, params=None, cli_args=None):
        self.entry = {
            "stage": stage_name,
            "entrypoint": entrypoint,
            "pid": os.getpid(),
            "started_utc": _utc_now(),
            "ended_utc": None,
            "wall_seconds": None,
            "status": "running",
            "exit_code": None,
            "error": None,
            # Copied, not referenced: `vars(args)` hands over argparse's live dict.
            "params": _jsonable(copy.deepcopy(params)) if params else {},
            "cli_args": _jsonable(copy.deepcopy(cli_args)) if cli_args else None,
            "inputs": [],
            "outputs": [],
            "metrics": {},
            "extra": {},
        }

    def add_input(self, path, role=None):
        self.entry["inputs"].append(file_record(path, role=role))

    def add_output(self, path, rows=None, role=None):
        self.entry["outputs"].append(file_record(path, role=role, rows=rows))

    def metric(self, key, value):
        self.entry["metrics"][key] = _jsonable(value)

    def extra(self, **fields):
        self.entry["extra"].update({k: _jsonable(v) for k, v in fields.items()})


class stage:
    """Context manager recording one pipeline stage into ``run.json``.

    Implemented as a class rather than via ``@contextmanager`` so that
    ``SystemExit`` — which ``score.py`` and ``annotator.py`` raise through
    ``sys.exit(1)`` in around ten places — is caught and recorded.  A generator
    context manager catching only ``Exception`` would file those runs as
    successful.
    """

    def __init__(self, output_dir, stage_name, entrypoint=None,
                 params=None, cli_args=None, run_id=None):
        self.output_dir = str(output_dir)
        self.stage_name = stage_name
        self.run_id = run_id or get_run_id()
        self.record = StageRecord(stage_name, entrypoint=entrypoint,
                                  params=params, cli_args=cli_args)
        self._started = time.time()

    def __enter__(self):
        return self.record

    def __exit__(self, exc_type, exc_value, tb):
        entry = self.record.entry
        entry["ended_utc"] = _utc_now()
        entry["wall_seconds"] = round(time.time() - self._started, 3)

        if exc_type is None:
            entry["status"] = "ok"
        elif isinstance(exc_value, SystemExit):
            code = exc_value.code
            code = 0 if code is None else code
            entry["exit_code"] = code if isinstance(code, int) else 1
            entry["status"] = "ok" if entry["exit_code"] == 0 else "error"
        else:
            entry["status"] = "error"
            entry["error"] = {
                "type": exc_type.__name__,
                "message": str(exc_value),
                "traceback": "".join(traceback.format_exception(
                    exc_type, exc_value, tb))[-MAX_TRACEBACK_CHARS:],
            }

        self._merge()
        return False  # never suppress; the caller's failure is still a failure

    def _merge(self):
        """Fold this stage into the manifest, warning rather than raising."""
        try:
            os.makedirs(self.output_dir, exist_ok=True)
            manifest_path = os.path.join(self.output_dir, RUN_JSON_FILENAME)
            document = _load(manifest_path)

            run = document["runs"].setdefault(self.run_id, {})
            if not run:
                run.update({
                    "run_id": self.run_id,
                    "created_utc": _utc_now(),
                    "proximate": proximate_version(),
                    "environment": environment_snapshot(),
                    "datasets_build_info": read_build_info(),
                    "stages": [],
                })
            run["stages"].append(self.record.entry)

            _write(manifest_path, document)
            logger.info("Recorded %s stage in run.json (status=%s, %.1fs)",
                        self.stage_name, self.record.entry["status"],
                        self.record.entry["wall_seconds"])
        except Exception:
            logger.warning("Could not record the %s stage in %s.",
                           self.stage_name, self.output_dir, exc_info=True)
