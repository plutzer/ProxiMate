"""The operations layer: every dataset action and analysis the GUI and the MCP
server share, as plain functions over the output directory and the datasets store.

Nothing here touches Shiny.  Callers own progress bars and notifications; this
module owns validation, the job lock per dataset, provenance and the store.

Operations come in three kinds:

- **dataset**: ``run_parse``, ``run_score``, ``write_feature_enrichment``,
  ``delete_dataset``, ``rename_dataset``, ``load_session``.  They create or change a
  dataset directory, take that dataset's job lock (a second caller is refused with
  ``BusyError`` rather than racing) and update the datasets store, so every open
  browser session sees the result.
- **sandbox**: ``threshold_metrics``, ``feature_analysis``, ``compare_networks``.  They
  read a dataset directory and return their result.  They write nothing under it, so
  what a run directory holds is exactly what its dataset operations produced.
- **read**: ``dataset_info`` and the store's listings.

Every threshold argument is an explicit ``{'SaintScore', 'BFDR', 'WD', 'WDFDR'}``
dict; no operation reads a setting from the GUI.
"""

import contextlib
import datetime
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading

import pandas as pd

import cytoscape_ctl
import dataset_store as store
import log_config
import parse
import provenance
import session_archive
from Ann_Enrichment import process_refactored
from network_comparison import calculate_volcano_data, load_and_filter_bait_data
from QC_plots import apply_score_thresholds, calculate_threshold_metrics
from setup_datasets import ORGANISMS
from log_config import get_logger

logger = get_logger(__name__)

SCRIPTS_DIR = os.path.dirname(os.path.abspath(parse.__file__))
THRESHOLD_KEYS = ('SaintScore', 'BFDR', 'WD', 'WDFDR')
THRESHOLD_RANGES = {'SaintScore': (0.0, 1.0), 'BFDR': (0.0, 1.0),
                    'WD': (0.0, None), 'WDFDR': (0.0, 1.0)}
INPUT_FORMATS = {
    # format: (required file keys, quant type fixed by the format or None)
    'MaxQuant': (('pg', 'ed'), None),
    'DIA-NN': (('matrix', 'ed'), 'Intensity'),
    'Pioneer': (('matrix', 'ed'), 'Intensity'),
    'FragPipe': (('fp', 'ed'), None),
    'MSstats': (('msstats', 'ed'), 'Intensity'),
    'SAINT': (('bait', 'prey', 'interaction'), None),
}
QUANT_TYPES = ('Intensity', 'LFQ', 'Spectral Counts')
IMPUTATION_LABELS = {0: 'Default', 1: 'Prey-specific', 2: 'Two-component AFT', 3: 'One-component AFT'}
PI_METHODS = ('weighted_average', 'single_bait')
FEATURE_TYPES = ('GO_CC', 'GO_BP', 'GO_MF', 'Motifs', 'Regions', 'Repeats', 'Compositions', 'Domains')
RESULT_FILES = ('ED.csv', 'interaction.txt', 'merged.csv', 'annotated_scores.csv',
                'Feature_enrichment.csv', 'run.json')
# The messages each long operation passes to its ``progress`` callback, in order, so a
# caller can list them ahead of the run and mark the one under way.
PARSE_STAGES = ('Parsing the input files', 'Counting what was written', 'Recording the dataset')
SCORE_STAGES = ('Running SAINT and CompPASS', 'Adding protein annotation', 'Updating the datasets table')
FEATURE_STAGES = ('Loading the scored table', 'Testing features for enrichment',
                  'Writing Feature_enrichment.csv')

OUT_DIR = None


class StageError(RuntimeError):
    """A pipeline stage exited nonzero; the message carries what it logged."""


class BusyError(RuntimeError):
    """The dataset is being worked on by another caller."""


def configure(out_dir):
    """Point the backend and the datasets store at one output directory."""
    global OUT_DIR
    OUT_DIR = out_dir
    os.makedirs(out_dir, exist_ok=True)
    store.configure(out_dir)


# --- validation and paths ---------------------------------------------------------

def imputation_choices(quant_type):
    """Imputation modes the GUI offers for a dataset, keyed as ``score.py`` takes them.
    AFT imputation models intensities, so spectral counts get the default only."""
    choices = {0: IMPUTATION_LABELS[0]}
    if quant_type != 'Spectral Counts':
        choices.update({k: IMPUTATION_LABELS[k] for k in (2, 3)})
    return choices


def validate_thresholds(thresholds):
    """The four score cutoffs as floats; anything missing, extra or out of range raises."""
    if not isinstance(thresholds, dict):
        raise ValueError(f"thresholds must be a dict with keys {THRESHOLD_KEYS}")
    missing = [k for k in THRESHOLD_KEYS if k not in thresholds]
    extra = [k for k in thresholds if k not in THRESHOLD_KEYS]
    if missing or extra:
        raise ValueError(f"thresholds need exactly {THRESHOLD_KEYS}; "
                         f"missing {missing}, unexpected {extra}")
    out = {}
    for key in THRESHOLD_KEYS:
        try:
            value = float(thresholds[key])
        except (TypeError, ValueError):
            raise ValueError(f"threshold {key} must be a number, got {thresholds[key]!r}")
        low, high = THRESHOLD_RANGES[key]
        if value < low or (high is not None and value > high):
            raise ValueError(f"threshold {key}={value} is outside [{low}, {high if high is not None else 'inf'}]")
        out[key] = value
    return out


def dataset_dir(name):
    if OUT_DIR is None:
        raise RuntimeError("backend.configure(out_dir) has not been called")
    name = str(name)
    if not re.fullmatch(r"[A-Za-z0-9_]+", name):
        raise ValueError(f"not a dataset name: {name!r}")
    return os.path.join(OUT_DIR, name)


def results_path(name):
    return os.path.join(dataset_dir(name), 'annotated_scores.csv')


def _require_scored(name):
    path = results_path(name)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"dataset {name!r} has no annotated_scores.csv; score it first")
    return path


def _require_bait(scores, bait):
    known = sorted(scores['Experiment.ID'].astype(str).unique())
    if str(bait) not in known:
        raise ValueError(f"no bait {bait!r} in this dataset; baits are {known}")


# --- jobs ---------------------------------------------------------------------------

_JOBS = {'lock': threading.Lock(), 'running': {}}


@contextlib.contextmanager
def _job(name, actor, what):
    with _JOBS['lock']:
        holder = _JOBS['running'].get(name)
        if holder:
            raise BusyError(f"dataset {name!r} is busy: {holder['what']} by {holder['actor']} "
                            f"since {holder['since']}")
        _JOBS['running'][name] = {'actor': actor, 'what': what,
                                  'since': datetime.datetime.now().isoformat(timespec='seconds')}
    try:
        yield
    finally:
        with _JOBS['lock']:
            _JOBS['running'].pop(name, None)


def running_jobs():
    with _JOBS['lock']:
        return {name: dict(job) for name, job in _JOBS['running'].items()}


# --- parse -----------------------------------------------------------------------

def _write_frame(frame):
    handle = tempfile.NamedTemporaryFile("w", suffix=".csv", delete=False, newline="", encoding="utf-8")
    with handle:
        frame.to_csv(handle, index=False)
    return handle.name


def run_parse(name, input_format, files, quant_type, actor='gui', ed_frame=None, bait_frame=None,
              run_id=None, progress=None):
    """Parse one dataset into ``<out_dir>/<name>/`` and add its row to the store.

    ``files`` maps the format's file keys (see ``INPUT_FORMATS``) to paths.  ``ed_frame``
    replaces the design file with an edited table; ``bait_frame`` does the same for
    the SAINT bait table.  An empty frame means the file.  ``progress(message, fraction)`` is called with each of
    ``PARSE_STAGES``.  Returns the new row plus ``n_preys`` and ``n_interactions``, the
    rows written to prey.txt and interaction.txt.
    """
    if input_format not in INPUT_FORMATS:
        raise ValueError(f"unknown input format {input_format!r}; one of {list(INPUT_FORMATS)}")
    required, fixed_quant = INPUT_FORMATS[input_format]
    missing = [k for k in required if not files.get(k)]
    if missing:
        raise ValueError(f"{input_format} needs the {', '.join(missing)} file(s)")
    for key in required:
        if not os.path.isfile(files[key]):
            raise ValueError(f"{key} file not found: {files[key]}")
    quant_type = fixed_quant or quant_type
    if quant_type not in QUANT_TYPES:
        raise ValueError(f"quant_type must be one of {QUANT_TYPES}, got {quant_type!r}")

    taken = set(store.names()) | set(session_archive.dataset_directories(OUT_DIR))
    problem = parse.validate_name(name, taken)
    if problem != 0:
        raise ValueError(problem)
    output_path = dataset_dir(name)
    run_id = run_id or log_config.new_run_id()
    report = progress or (lambda message, value: None)

    def design_path():
        if ed_frame is not None and not ed_frame.empty:
            logger.info("Parsing the edited experimental design table (%d rows)", len(ed_frame))
            return _write_frame(ed_frame), True
        return files['ed'], False

    with _job(name, actor, 'parse'), log_config.run_context(run_id), log_config.dataset_log(output_path):
        logger.info("Parsing dataset '%s' (format=%s, actor=%s)", name, input_format, actor)
        report(PARSE_STAGES[0], 0.1)
        if input_format == 'SAINT':
            if bait_frame is None or bait_frame.empty:
                bait_frame = pd.read_csv(files['bait'], sep="\t", header=None,
                                         names=["Experiment Name", "Bait", "Type"])
                bait_frame['Bait ID'] = 'None'
            os.makedirs(output_path, exist_ok=True)
            n_exp, n_ctrl = parse.parse_from_saint(bait_frame, files['prey'], files['interaction'],
                                                   output_path)
            for key in ('bait', 'prey', 'interaction'):
                shutil.copy(files[key], os.path.join(output_path, f'{key}.txt'))
        else:
            ed_path, temporary = design_path()
            try:
                if input_format == 'MaxQuant':
                    n_exp, n_ctrl = parse.parse_ed_pg(files['pg'], ed_path, quant_type, output_path)
                elif input_format == 'DIA-NN':
                    n_exp, n_ctrl = parse.parse_diann(files['matrix'], ed_path, quant_type, output_path)
                elif input_format == 'Pioneer':
                    n_exp, n_ctrl = parse.parse_pioneer(files['matrix'], ed_path, quant_type, output_path)
                elif input_format == 'FragPipe':
                    n_exp, n_ctrl = parse.parse_fragpipe(files['fp'], ed_path, quant_type, output_path)
                else:
                    n_exp, n_ctrl = parse.parse_msstats(files['msstats'], ed_path, output_path)
            finally:
                if temporary:
                    os.unlink(ed_path)
        report(PARSE_STAGES[1], 0.8)
        counts = {'n_preys': _line_count(os.path.join(output_path, 'prey.txt')),
                  'n_interactions': _line_count(os.path.join(output_path, 'interaction.txt'))}
        report(PARSE_STAGES[2], 0.9)
        row = {'Dataset Name': name, 'Input Type': input_format, 'Quant Type': quant_type,
               'Experiments': n_exp, 'Controls': n_ctrl, 'Scored': '', 'Imputation': '',
               'WDFDR iterations': ''}
        store.append(row)
        logger.info("Parsed '%s': %d experiments, %d controls, %d preys, %d interactions",
                    name, n_exp, n_ctrl, counts['n_preys'], counts['n_interactions'])
    return {**row, **counts}


def _line_count(path):
    with open(path, encoding='utf-8') as handle:
        return sum(1 for _ in handle)


# --- score -----------------------------------------------------------------------

def _log_size(dataset_path):
    try:
        return os.path.getsize(os.path.join(dataset_path, "proximate.log"))
    except OSError:
        return 0


def _log_tail_since(dataset_path, offset, max_chars=1500):
    """The log text a child process appended after ``offset``.

    Reading from an offset rather than the end of the file matters when the child
    fails before it attaches its own handler: the tail would then be the previous
    run's output, a plausible but wrong explanation of this failure.
    """
    try:
        with open(os.path.join(dataset_path, "proximate.log"), encoding="utf-8") as handle:
            handle.seek(offset)
            new_text = handle.read().strip()
    except OSError:
        return ""
    if len(new_text) > max_chars:
        new_text = "..." + new_text[-max_chars:]
    return new_text


def _run_stage_subprocess(command, dataset_path, run_id):
    """Run a pipeline stage, streaming its output to this process's terminal.

    The child configures logging the same way this process does and writes to the
    dataset's own log, so nothing needs re-logging here.  Returns (returncode, text
    the child appended to the dataset log).
    """
    child_env = dict(os.environ, PROXIMATE_RUN_ID=run_id)
    logger.info("Running: %s", " ".join(command))
    offset = _log_size(dataset_path)
    result = subprocess.run(command, env=child_env)
    return result.returncode, _log_tail_since(dataset_path, offset)


def run_score(name, imputation, wdfdr_iterations, organism, exclude_hcm, pi_method=None,
              pi_bait=None, aft_min_obs=None, seed=None, actor='gui', run_id=None, progress=None):
    """Score and annotate a parsed dataset; on success its row reads Scored = Yes.

    ``imputation`` is 0-3 as ``score.py`` takes it; ``pi_method``/``pi_bait`` apply to
    imputation 2 only and ``aft_min_obs`` (the sigma floor, None for score.py's default
    of 4) to imputation 2 and 3.  ``progress`` is called with each of ``SCORE_STAGES``.
    Raises ``StageError`` with the stage's log when either script fails, and
    ``BusyError`` when the dataset is already being worked on.
    """
    row = store.row(name)
    imputation = int(imputation)
    if imputation not in IMPUTATION_LABELS:
        raise ValueError(f"imputation must be one of {sorted(IMPUTATION_LABELS)}, got {imputation}")
    if row['Quant Type'] == 'Spectral Counts' and imputation != 0:
        raise ValueError(
            f"dataset {name!r} holds spectral counts, and AFT imputation applies to "
            "intensity data only; score it with imputation 0 (Default)")
    if organism not in ORGANISMS:
        raise ValueError(f"organism must be one of {sorted(ORGANISMS)}, got {organism!r}")
    if pi_method is not None and pi_method not in PI_METHODS:
        raise ValueError(f"pi_method must be one of {PI_METHODS}, got {pi_method!r}")
    if pi_method == 'single_bait' and not pi_bait:
        raise ValueError("pi_method 'single_bait' needs pi_bait")
    if aft_min_obs is not None and int(aft_min_obs) < 0:
        raise ValueError("aft_min_obs must be >= 0")
    wdfdr_iterations = int(wdfdr_iterations)
    if wdfdr_iterations < 0:
        raise ValueError("wdfdr_iterations must be >= 0")
    dataset_path = dataset_dir(name)
    ed_path = os.path.join(dataset_path, 'ED.csv')
    if not os.path.isfile(ed_path):
        raise FileNotFoundError(f"dataset {name!r} has no ED.csv; parse it first")
    run_id = run_id or log_config.new_run_id()
    report = progress or (lambda message, value: None)

    score_cmd = [sys.executable, os.path.join(SCRIPTS_DIR, 'score.py'),
                 '--experimentalDesign', ed_path, '--scoreInputs', dataset_path,
                 '--outputPath', dataset_path, '--n-iterations', str(wdfdr_iterations),
                 '--imputation', str(imputation), '--quantType', str(row['Quant Type'])]
    if imputation == 2 and pi_method:
        score_cmd += ['--pi-method', pi_method]
        if pi_method == 'single_bait':
            score_cmd += ['--pi-bait', str(pi_bait)]
    if imputation in (2, 3) and aft_min_obs is not None:
        score_cmd += ['--aft-min-obs', str(int(aft_min_obs))]
    if seed is not None:
        score_cmd += ['--seed', str(int(seed))]
    ann_cmd = [sys.executable, os.path.join(SCRIPTS_DIR, 'annotator.py'),
               '--organism', organism, '--scoreFile', os.path.join(dataset_path, 'merged.csv'),
               '--outputDir', dataset_path]
    if organism == 'human' and exclude_hcm:
        ann_cmd.append('--excludeHCM')

    with _job(name, actor, 'score'), log_config.run_context(run_id), log_config.dataset_log(dataset_path):
        logger.info("Starting scoring for dataset '%s' (quant=%s, imputation=%d, actor=%s)",
                    name, row['Quant Type'], imputation, actor)
        report(SCORE_STAGES[0], 0.1)
        code, tail = _run_stage_subprocess(score_cmd, dataset_path, run_id)
        if code != 0:
            logger.error("score.py failed (exit code %d)", code)
            raise StageError(f"Scoring failed for '{name}' (exit code {code}):\n"
                             f"{tail or 'No output was logged; check the server log.'}")
        if not os.path.isfile(os.path.join(dataset_path, 'merged.csv')):
            raise StageError(f"Scoring failed for '{name}': merged.csv was not produced")
        report(SCORE_STAGES[1], 0.65)
        code, tail = _run_stage_subprocess(ann_cmd, dataset_path, run_id)
        if code != 0:
            logger.error("annotator.py failed (exit code %d)", code)
            raise StageError(f"Annotation failed for '{name}' (exit code {code}):\n"
                             f"{tail or 'No output was logged; check the server log.'}")
        report(SCORE_STAGES[2], 0.9)
        store.update(name, Scored='Yes', Imputation=IMPUTATION_LABELS[imputation],
                     **{'WDFDR iterations': wdfdr_iterations},
                     **session_archive.scored_fields(dataset_path))
        logger.info("Scoring and annotation completed for '%s'", name)
    return {'dataset': name, 'run_id': run_id, 'organism': organism, 'exclude_hcm': bool(exclude_hcm),
            'imputation': IMPUTATION_LABELS[imputation], 'wdfdr_iterations': wdfdr_iterations}


# --- dataset actions ----------------------------------------------------------------

def write_feature_enrichment(name, thresholds, actor='gui', progress=None):
    """Feature enrichment of every bait's passing preys at ``thresholds``, written to
    the dataset's Feature_enrichment.csv.  ``progress`` is called with each of
    ``FEATURE_STAGES``.  Returns the table."""
    thresholds = validate_thresholds(thresholds)
    path = _require_scored(name)
    report = progress or (lambda message, value: None)
    with _job(name, actor, 'feature enrichment'), log_config.dataset_log(dataset_dir(name)):
        logger.info("Starting feature enrichment for '%s' (thresholds=%s, actor=%s)",
                    name, thresholds, actor)
        report(FEATURE_STAGES[0], 0.1)
        scores = pd.read_csv(path)
        report(FEATURE_STAGES[1], 0.3)
        result = process_refactored(scores, list(FEATURE_TYPES), thresholds)
        report(FEATURE_STAGES[2], 0.9)
        result.to_csv(os.path.join(dataset_dir(name), 'Feature_enrichment.csv'), index=False)
        logger.info("Feature enrichment for '%s': %d rows written", name, len(result))
    return result


def delete_dataset(name, actor='gui'):
    """Remove the dataset's directory and its row.  A network drawn from it stays in
    Cytoscape, but the controller no longer names it as the drawn dataset."""
    store.row(name)
    path = dataset_dir(name)
    with _job(name, actor, 'delete'):
        logger.info("Deleting dataset '%s' (actor=%s)", name, actor)
        shutil.rmtree(path)
        store.remove(name)
    cytoscape_ctl.forget_dataset(name, actor=actor)


def _rename_paths(value, old_path, new_path):
    """``value`` with every string under ``old_path`` moved under ``new_path``."""
    if isinstance(value, dict):
        return {k: _rename_paths(v, old_path, new_path) for k, v in value.items()}
    if isinstance(value, list):
        return [_rename_paths(v, old_path, new_path) for v in value]
    if isinstance(value, str) and (value == old_path or value.startswith(old_path + os.sep)):
        return new_path + value[len(old_path):]
    return value


def rename_dataset(old, new, actor='gui'):
    """Rename a dataset: its directory, the paths its run.json records, its row and,
    when it is the drawn one, the Cytoscape network title.  Returns the row."""
    store.row(old)
    taken = set(store.names()) | set(session_archive.dataset_directories(OUT_DIR))
    problem = parse.validate_name(new, taken)
    if problem != 0:
        raise ValueError(problem)
    old_path, new_path = dataset_dir(old), dataset_dir(new)
    with _job(old, actor, 'rename'):
        logger.info("Renaming dataset '%s' to '%s' (actor=%s)", old, new, actor)
        os.rename(old_path, new_path)
        manifest = os.path.join(new_path, provenance.RUN_JSON_FILENAME)
        if os.path.isfile(manifest):
            with open(manifest, encoding='utf-8') as handle:
                document = json.load(handle)
            with open(manifest, 'w', encoding='utf-8') as handle:
                json.dump(_rename_paths(document, old_path, new_path), handle, indent=2)
        store.update(old, **{'Dataset Name': new})
    cytoscape_ctl.rename_dataset(old, new, actor=actor)
    return store.row(new)


# --- session -----------------------------------------------------------------------

def clear_datasets():
    """Remove every dataset directory and file under the output directory.

    The operational log stays: it must outlive the datasets it describes, and its
    handler is still writing to it.
    """
    if running_jobs():
        raise BusyError(f"jobs are running: {running_jobs()}")
    store.clear()
    logger.info("Clearing all datasets under %s", OUT_DIR)
    keep = (log_config.SERVER_LOG_FILENAME,)
    for root, dirs, files in os.walk(OUT_DIR):
        for file in files:
            if root == OUT_DIR and file.startswith(keep):
                continue
            path = os.path.join(root, file)
            try:
                os.remove(path)
            except OSError:
                logger.exception("Could not remove %s", path)
        for directory in dirs:
            path = os.path.join(root, directory)
            try:
                shutil.rmtree(path)
            except OSError:
                logger.exception("Could not remove %s", path)


def load_session(zip_path, actor='gui'):
    """Replace the session with the datasets in a session archive.

    The archive is checked before anything is removed, so a wrong file leaves the
    session as it was.  Returns the restored table.
    """
    session_archive.inspect_session_archive(zip_path)
    clear_datasets()
    logger.info("Restoring session from archive %s (actor=%s)", zip_path, actor)
    table = session_archive.extract_session_archive(zip_path, OUT_DIR)
    store.replace(table)
    logger.info("Session restored: %d datasets", len(table))
    return table


# --- read ----------------------------------------------------------------------------

def dataset_info(name):
    """The store row, which result files exist, the baits, and the run manifest's stages."""
    row = store.row(name)
    root = dataset_dir(name)
    files = {f: os.path.isfile(os.path.join(root, f)) for f in RESULT_FILES}
    baits = []
    if files['annotated_scores.csv']:
        scores = pd.read_csv(results_path(name), usecols=['Experiment.ID'])
        baits = sorted(scores['Experiment.ID'].astype(str).unique())
    settings = provenance.annotation_settings(root) if files['run.json'] else None
    return {'row': row, 'files': files, 'baits': baits, 'annotation': settings,
            'busy': running_jobs().get(name)}


STAGE_FIELDS = ('stage', 'entrypoint', 'status', 'started_utc', 'wall_seconds', 'params',
                'metrics', 'extra', 'outputs')


def run_info(name, last_n=5):
    """The last ``last_n`` runs in the dataset's manifest, oldest first, each with its
    stages' status, timing, parameters, metrics and error message.  Environment
    snapshots and tracebacks stay in run.json."""
    path = os.path.join(dataset_dir(name), provenance.RUN_JSON_FILENAME)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"dataset {name!r} has no run.json")
    with open(path, encoding='utf-8') as handle:
        document = json.load(handle)
    runs = []
    for run in list(document['runs'].values())[-int(last_n):]:
        stages = []
        for entry in run.get('stages', []):
            stage = {key: entry.get(key) for key in STAGE_FIELDS}
            stage['outputs'] = [o.get('path') for o in entry.get('outputs', [])]
            stage['error'] = (entry.get('error') or {}).get('message')
            stages.append(stage)
        runs.append({'run_id': run.get('run_id'), 'created_utc': run.get('created_utc'),
                     'proximate': run.get('proximate'), 'stages': stages})
    return {'dataset': name, 'path': path, 'n_runs': len(document['runs']), 'runs': runs}


def log_tail(name, n_lines=50):
    """The last ``n_lines`` of the dataset's proximate.log."""
    path = os.path.join(dataset_dir(name), 'proximate.log')
    if not os.path.isfile(path):
        raise FileNotFoundError(f"dataset {name!r} has no proximate.log")
    with open(path, encoding='utf-8', errors='replace') as handle:
        lines = [line.rstrip('\n') for line in handle]
    return {'dataset': name, 'path': path, 'n_lines_total': len(lines), 'lines': lines[-int(n_lines):]}


# --- sandbox -----------------------------------------------------------------------------

SCORE_COLUMNS = ('Experiment.ID', 'Prey.ID', 'First_ID', 'First_Prey_Gene', 'SaintScore', 'BFDR',
                 'WD', 'WDFDR', 'FoldChange', 'AvgIntensity', 'AvgSpec', 'In.BioGRID', 'Multivalidated')
# The QC tab's results table: display name -> annotated_scores.csv column.
RESULTS_TABLE_COLUMNS = {'Bait': 'Experiment.ID', 'Prey gene': 'First_Prey_Gene', 'Prey ID': 'First_ID',
                         'SAINT': 'SaintScore', 'BFDR': 'BFDR', 'WD': 'WD', 'WDFDR': 'WDFDR',
                         'Fold change': 'FoldChange', 'Abundance': None, 'In BioGRID': 'In.BioGRID',
                         'Multivalidated': 'Multivalidated'}


def passing_scores(name, thresholds, baits=None):
    """The scored rows passing ``thresholds``, restricted to ``baits`` when given,
    ordered by bait, then SAINT score descending, then BFDR."""
    thresholds = validate_thresholds(thresholds)
    scores = pd.read_csv(_require_scored(name))
    if baits:
        for bait in baits:
            _require_bait(scores, bait)
        scores = scores[scores['Experiment.ID'].astype(str).isin([str(b) for b in baits])]
    passing = apply_score_thresholds(scores, thresholds)
    columns = [c for c in SCORE_COLUMNS if c in passing.columns]
    return (passing[columns]
            .sort_values(['Experiment.ID', 'SaintScore', 'BFDR'], ascending=[True, False, True], kind='stable')
            .reset_index(drop=True))


def results_table(rows):
    """``passing_scores`` rows under the display names of ``RESULTS_TABLE_COLUMNS``,
    numbers rounded to 3 decimals.  Abundance is AvgSpec or AvgIntensity, whichever
    the scores carry."""
    abundance = [c for c in ('AvgSpec', 'AvgIntensity') if c in rows.columns]
    if not abundance:
        raise ValueError("the scores carry neither AvgSpec nor AvgIntensity")
    columns = {**RESULTS_TABLE_COLUMNS, 'Abundance': abundance[0]}
    return pd.DataFrame({label: rows[column] for label, column in columns.items()}).round(3)


def threshold_metrics(name, thresholds, bait=None, biogrid_path=None):
    """The QC tab's three metrics at one threshold set, for every bait or one."""
    thresholds = validate_thresholds(thresholds)
    path = _require_scored(name)
    experiments = None
    if bait is not None:
        _require_bait(pd.read_csv(path, usecols=['Experiment.ID']), bait)
        experiments = [str(bait)]
    metrics = calculate_threshold_metrics(path, thresholds, ctrl_experiments=experiments,
                                          biogrid_path=biogrid_path)
    out = {k: (None if v is None else float(v)) for k, v in metrics.items()}
    out.update(dataset=name, bait=bait, thresholds=thresholds)
    return out


def feature_analysis(name, thresholds, feature_types=None):
    """Feature enrichment for every bait's foreground at these thresholds, as a table.

    Nothing is written: the GUI's own run is what lands in ``Feature_enrichment.csv``.
    """
    thresholds = validate_thresholds(thresholds)
    feature_types = list(feature_types or FEATURE_TYPES)
    unknown = [f for f in feature_types if f not in FEATURE_TYPES]
    scores = pd.read_csv(_require_scored(name))
    absent = [f for f in feature_types if f not in scores.columns]
    if unknown and absent:
        raise ValueError(f"unknown feature types {unknown}; available: {list(FEATURE_TYPES)}")
    if absent:
        raise ValueError(f"feature types not in this dataset's annotation: {absent}")
    return process_refactored(scores, feature_types, thresholds)


PREY_ANNOTATION_COLUMNS = ('first_SCL', 'Main location', 'GO_CC', 'Human_Complex')


def prey_annotations(name, thresholds, ids=None):
    """One row per prey: identifiers, the baits it passes ``thresholds`` under, the
    baits BioGRID already links it to, its best scores and its annotation columns.

    ``ids`` (accessions or gene symbols, case-insensitive) restrict the rows.  Returns
    ``(frame, unmatched)``; an id absent from the dataset is a normal answer, so it is
    listed rather than raised.  Annotation columns an organism lacks are left out.
    """
    thresholds = validate_thresholds(thresholds)
    scores = pd.read_csv(_require_scored(name))
    scores['First_ID'] = scores['First_ID'].astype(str)
    unmatched = []
    if ids:
        keys = [scores[c].astype(str).str.lower()
                for c in ('Prey.ID', 'First_ID', 'First_Prey_Gene')]
        wanted = {str(i).lower() for i in ids}
        hit = keys[0].isin(wanted) | keys[1].isin(wanted) | keys[2].isin(wanted)
        found = set().union(*(set(k[hit]) for k in keys))
        unmatched = [i for i in ids if str(i).lower() not in found]
        scores = scores[hit]

    def baits_of(frame):
        return frame.groupby('First_ID')['Experiment.ID'].agg(lambda s: sorted(set(s)))

    groups = scores.groupby('First_ID', sort=False)
    rows = groups.agg(accession=('Prey.ID', 'first'), gene=('First_Prey_Gene', 'first'),
                      n_baits_seen=('Experiment.ID', 'nunique'),
                      max_saint=('SaintScore', 'max'), max_fold_change=('FoldChange', 'max'))
    lists = {'passing_baits': baits_of(apply_score_thresholds(scores, thresholds))}
    if 'In.BioGRID' in scores.columns:
        lists['known_baits'] = baits_of(scores[scores['In.BioGRID'].eq(True)])
    for column, baits in lists.items():
        rows[column] = [baits.get(prey, []) for prey in rows.index]
    for column in PREY_ANNOTATION_COLUMNS:
        if column in scores.columns:
            rows[column] = groups[column].first()
    rows['n_passing'] = rows['passing_baits'].str.len()
    rows = rows.sort_values(['n_passing', 'max_saint'], ascending=False, kind='stable')
    return rows.reset_index(), unmatched


def compare_sets(data_a, data_b):
    """Gene lists for the Venn regions of two filtered bait tables."""
    set_a = set(data_a['Prey.ID'].unique()) if len(data_a) else set()
    set_b = set(data_b['Prey.ID'].unique()) if len(data_b) else set()

    def genes(frame, ids):
        if not ids or not len(frame):
            return []
        return sorted(frame[frame['Prey.ID'].isin(ids)]['First_Prey_Gene'].astype(str).unique())
    return {'a_only': genes(data_a, set_a - set_b), 'b_only': genes(data_b, set_b - set_a),
            'both': genes(data_a, set_a & set_b)}


def compare_networks(name, bait_a, bait_b, thresholds_a, thresholds_b):
    """The Network Comparison tab's volcano table and Venn gene lists."""
    thresholds_a = validate_thresholds(thresholds_a)
    thresholds_b = validate_thresholds(thresholds_b)
    scores = pd.read_csv(_require_scored(name), usecols=['Experiment.ID'])
    _require_bait(scores, bait_a)
    _require_bait(scores, bait_b)
    data_a = load_and_filter_bait_data(name, bait_a, thresholds_a, OUT_DIR)
    data_b = load_and_filter_bait_data(name, bait_b, thresholds_b, OUT_DIR)
    volcano = calculate_volcano_data(name, bait_a, bait_b, thresholds_a, thresholds_b, OUT_DIR)
    return {'dataset': name, 'bait_a': bait_a, 'bait_b': bait_b,
            'thresholds_a': thresholds_a, 'thresholds_b': thresholds_b,
            'genes': compare_sets(data_a, data_b), 'volcano': volcano}
