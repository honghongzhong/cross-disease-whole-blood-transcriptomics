"""Run this comparison through the shared processing workflow."""

from pathlib import Path
import runpy
import sys


if __name__ == "__main__":
    comparison = Path(__file__).resolve().parent
    sys.argv[1:1] = ["--contrast", comparison.name]
    runpy.run_path(
        str(comparison.parent / "_shared" / "reproduce.py"),
        run_name="__main__",
    )
