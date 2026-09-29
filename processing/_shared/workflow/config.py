"""Portable locations shared by the Python processing modules."""

import os
from pathlib import Path


def work_root() -> Path:
    """Return the isolated work tree selected by the processing launcher."""
    value = os.environ.get("BIOLOGY_WORK_ROOT")
    if not value:
        raise RuntimeError("BIOLOGY_WORK_ROOT is required")
    return Path(value).resolve()


def rscript_command() -> str:
    """Use the caller's Rscript executable or the executable on PATH."""
    return os.environ.get("BIOLOGY_RSCRIPT", "Rscript")
