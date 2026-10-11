"""Logging to stdout and to a log file, configured once by each script."""

import logging
import sys
from pathlib import Path


def configure_logging(log_path, log_level=logging.INFO):
    """Send every logger's records to stdout and to `log_path` (appended to)."""
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    root_logger.handlers.clear()

    formatter = logging.Formatter(fmt="%(asctime)s %(name)s %(levelname)s %(message)s",
                                  datefmt="%Y-%m-%d %H:%M:%S")

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)

    root_logger.addHandler(file_handler)
    root_logger.addHandler(stream_handler)
