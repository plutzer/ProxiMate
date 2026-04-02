"""
Centralized logging configuration for ProxiMate.

Usage in any module:
    from log_config import get_logger
    logger = get_logger(__name__)
    logger.info("message")

The first call to get_logger() initializes the root logger with:
- A StreamHandler writing to stderr (captured by Docker)
- An optional FileHandler if LOG_FILE env var is set

Log level defaults to INFO; set LOG_LEVEL env var to override (DEBUG, WARNING, etc.).
"""

import logging
import os
import sys

_initialized = False

LOG_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def setup_logging():
    """Configure the root logger. Safe to call multiple times (only runs once)."""
    global _initialized
    if _initialized:
        return
    _initialized = True

    level_name = os.environ.get("LOG_LEVEL", "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)

    root = logging.getLogger()
    root.setLevel(level)

    # Stderr handler — Docker captures this via `docker logs`
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(level)
    stderr_handler.setFormatter(logging.Formatter(LOG_FORMAT, datefmt=LOG_DATE_FORMAT))
    root.addHandler(stderr_handler)

    # Optional file handler — set LOG_FILE env var to enable
    log_file = os.environ.get("LOG_FILE")
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(logging.Formatter(LOG_FORMAT, datefmt=LOG_DATE_FORMAT))
        root.addHandler(file_handler)


def add_file_handler(log_path):
    """
    Add a FileHandler to the root logger that writes to log_path (append mode).

    Call this once an output directory is known, e.g.:
        add_file_handler(f"{output_dir}/proximate.log")

    Multiple scripts (parse, score, annotator) can append to the same file
    for a complete per-dataset log.
    """
    setup_logging()
    level_name = os.environ.get("LOG_LEVEL", "INFO").upper()
    level = getattr(logging, level_name, logging.INFO)

    file_handler = logging.FileHandler(log_path, mode='a')
    file_handler.setLevel(level)
    file_handler.setFormatter(logging.Formatter(LOG_FORMAT, datefmt=LOG_DATE_FORMAT))
    logging.getLogger().addHandler(file_handler)


def get_logger(name):
    """Get a named logger, initializing the root logger on first call."""
    setup_logging()
    return logging.getLogger(name)
