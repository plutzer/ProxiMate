"""
Centralized logging configuration for ProxiMate.

Usage in any module:
    from log_config import get_logger
    logger = get_logger(__name__)
    logger.info("message")

Logs are emitted on three destinations:

- **stderr**, captured by ``docker logs`` and visible in the terminal the
  container was started from.
- **A unified operational log**, ``proximate-server.log`` in ``PROXIMATE_LOG_DIR``
  (default ``/Outputs``), rotated so it cannot grow without bound.  This is the
  cross-dataset record of what the server did.
- **A per-dataset provenance log**, ``<outputDir>/proximate.log``, attached with
  :func:`add_file_handler` or :func:`dataset_log` once an output directory is
  known.  Parse, score and annotate all append to the same file, so it is a
  complete account of how one dataset's results were produced, and it travels
  with the results.

Environment variables:

- ``LOG_LEVEL`` — DEBUG, INFO, WARNING, ...  Defaults to INFO.  An unrecognized
  value is reported and INFO is used.
- ``PROXIMATE_LOG_DIR`` — directory for the unified operational log.
- ``LOG_FILE`` — an additional rotating log file at an explicit path.
- ``PROXIMATE_RUN_ID`` — identifier shared by every process of one run.  It is
  minted on first use and exported, so subprocesses inherit it through
  ``os.environ`` and their records can be correlated with the parent's.
"""

import contextlib
import contextvars
import logging
import logging.handlers
import os
import sys
import uuid
from datetime import datetime, timezone

PACKAGE = "proximate"

DEFAULT_LOG_DIR = "/Outputs"
SERVER_LOG_FILENAME = "proximate-server.log"
DATASET_LOG_FILENAME = "proximate.log"

RUN_ID_ENV_VAR = "PROXIMATE_RUN_ID"

# Keep the rotating logs bounded: the server log accumulates for the life of a
# container, and LOG_FILE may point at a long-lived path.
MAX_LOG_BYTES = 10 * 1024 * 1024
LOG_BACKUP_COUNT = 5

LOG_FORMAT = ("%(asctime)s [%(levelname)-7s] %(run_id)s pid=%(process)d "
              "%(name)s: %(message)s")
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

_initialized = False

# Set by `run_context` for the duration of one GUI action; see `get_run_id`.
_scoped_run_id = contextvars.ContextVar("proximate_run_id", default=None)

# Attached file handlers, keyed by normalized path, with a reference count so
# that overlapping `dataset_log` blocks do not detach each other's handler.
_file_handlers = {}
_file_handler_refs = {}

# Problems found while configuring logging, before any file handler exists.
# Replayed into each new file handler so a dataset's log records them too.
_startup_diagnostics = []


class _RunIDFilter(logging.Filter):
    """Stamp every record with the current run ID.

    This is a handler filter rather than a logger filter: records reach the
    package logger's handlers by propagation from module loggers, and logger
    filters are not applied to propagated records.
    """

    def filter(self, record):
        record.run_id = get_run_id()
        return True


def _formatter():
    return logging.Formatter(LOG_FORMAT, datefmt=LOG_DATE_FORMAT)


def _normalize(path):
    return os.path.realpath(os.fspath(path))


def _diagnose(level, message, *args):
    """Record and emit a problem with the logging configuration itself."""
    _startup_diagnostics.append((level, message, args))
    logging.getLogger(f"{PACKAGE}.log_config").log(level, message, *args)


def _resolve_level():
    """Return the configured level, reporting an unrecognized LOG_LEVEL."""
    name = os.environ.get("LOG_LEVEL", "INFO").upper()
    level = getattr(logging, name, None)
    if not isinstance(level, int):
        _diagnose(logging.WARNING,
                  "Unrecognized LOG_LEVEL %r; using INFO instead.", name)
        return logging.INFO
    return level


def _build_handler(handler, level):
    handler.setLevel(level)
    handler.setFormatter(_formatter())
    handler.addFilter(_RunIDFilter())
    return handler


def _attach_server_log(logger, level):
    """Attach the unified operational log.

    An explicitly configured directory that cannot be used is a misconfiguration
    worth reporting.  The default directory simply not existing is normal when
    running from a source checkout rather than the container, so it is not.
    """
    configured = os.environ.get("PROXIMATE_LOG_DIR")
    log_dir = configured or DEFAULT_LOG_DIR

    try:
        # An explicit destination is created on demand; the default is used only
        # where it already exists, so running from a source checkout does not
        # scatter an /Outputs directory across the host.
        if configured:
            os.makedirs(log_dir, exist_ok=True)
        elif not os.path.isdir(log_dir):
            raise OSError(f"{log_dir} does not exist")
        handler = logging.handlers.RotatingFileHandler(
            os.path.join(log_dir, SERVER_LOG_FILENAME),
            maxBytes=MAX_LOG_BYTES, backupCount=LOG_BACKUP_COUNT,
            encoding="utf-8")
    except OSError as error:
        if configured:
            _diagnose(logging.WARNING,
                      "Cannot write the unified log to %s (%s); "
                      "continuing with stderr only.", log_dir, error)
        else:
            _diagnose(logging.DEBUG,
                      "No unified log: %s is unavailable (%s).", log_dir, error)
        return
    logger.addHandler(_build_handler(handler, level))


def setup_logging():
    """Configure the ProxiMate logger. Safe to call multiple times."""
    global _initialized
    if _initialized:
        return
    _initialized = True

    logger = logging.getLogger(PACKAGE)
    # Handlers hang off this logger rather than the root logger, so LOG_LEVEL
    # controls ProxiMate's own output without also turning on DEBUG for
    # matplotlib, urllib3 and py4cytoscape.
    logger.propagate = False
    del logger.handlers[:]
    del _startup_diagnostics[:]

    level = _resolve_level()
    logger.setLevel(level)

    logger.addHandler(_build_handler(logging.StreamHandler(sys.stderr), level))
    _attach_server_log(logger, level)

    log_file = os.environ.get("LOG_FILE")
    if log_file:
        logger.addHandler(_build_handler(
            logging.handlers.RotatingFileHandler(
                log_file, maxBytes=MAX_LOG_BYTES,
                backupCount=LOG_BACKUP_COUNT, encoding="utf-8"),
            level))


def add_file_handler(log_path):
    """Append ProxiMate's log to `log_path` for the rest of the process.

    Repeated calls for the same file are ignored, so the long-lived GUI process
    can attach a dataset's log on every action without multiplying every
    subsequent line.  Use :func:`dataset_log` where the attachment should end
    with the action.
    """
    setup_logging()
    key = _normalize(log_path)
    if key in _file_handlers:
        return _file_handlers[key]

    level = _resolve_level()
    handler = _build_handler(
        logging.FileHandler(key, mode="a", encoding="utf-8"), level)
    logging.getLogger(PACKAGE).addHandler(handler)
    _file_handlers[key] = handler
    _file_handler_refs[key] = 1

    # A dataset's log should show any problem with the logging setup itself,
    # even though it was found before this file existed.
    for diag_level, message, args in _startup_diagnostics:
        handler.handle(logging.getLogger(f"{PACKAGE}.log_config").makeRecord(
            f"{PACKAGE}.log_config", diag_level, __file__, 0, message, args, None))
    return handler


def remove_file_handler(log_path):
    """Detach the handler for `log_path`, whatever its reference count."""
    key = _normalize(log_path)
    handler = _file_handlers.pop(key, None)
    _file_handler_refs.pop(key, None)
    if handler is not None:
        logging.getLogger(PACKAGE).removeHandler(handler)
        handler.close()


@contextlib.contextmanager
def dataset_log(output_dir):
    """Append to `output_dir`'s provenance log for the duration of the block.

    The handler is detached on the way out, including when the block raises, so
    a long-lived process does not keep writing to a dataset it has finished with.
    """
    setup_logging()
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, DATASET_LOG_FILENAME)
    key = _normalize(path)

    if key in _file_handlers:
        _file_handler_refs[key] += 1
    else:
        add_file_handler(path)
    try:
        yield path
    finally:
        _file_handler_refs[key] -= 1
        if _file_handler_refs[key] <= 0:
            remove_file_handler(key)


def new_run_id():
    """Mint a run identifier: a UTC timestamp plus a short random suffix."""
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"{stamp}-{uuid.uuid4().hex[:8]}"


def get_run_id():
    """Return the identifier of the run this code is part of.

    A run set by :func:`run_context` wins; otherwise the value inherited through
    the environment; otherwise one is minted and exported.  Exporting matters:
    subprocesses inherit ``os.environ``, which is how the GUI and the parse,
    score and annotate processes it launches share an ID.
    """
    scoped = _scoped_run_id.get()
    if scoped:
        return scoped

    run_id = os.environ.get(RUN_ID_ENV_VAR)
    if not run_id:
        run_id = new_run_id()
        os.environ[RUN_ID_ENV_VAR] = run_id
    return run_id


@contextlib.contextmanager
def run_context(run_id):
    """Stamp records emitted in this block with `run_id`.

    Held in a context variable rather than the environment so that concurrent
    Shiny sessions, which share one process and one event loop, do not overwrite
    each other's run.
    """
    token = _scoped_run_id.set(run_id)
    try:
        yield run_id
    finally:
        _scoped_run_id.reset(token)


def get_logger(name):
    """Get a logger under the ProxiMate package, initializing on first call."""
    setup_logging()
    if name == "__main__":
        # parse.py, score.py and annotator.py all run as scripts, so they would
        # otherwise share the name "__main__" in a log they all append to.
        name = os.path.splitext(os.path.basename(sys.argv[0]))[0] or "main"
    return logging.getLogger(f"{PACKAGE}.{name}")
