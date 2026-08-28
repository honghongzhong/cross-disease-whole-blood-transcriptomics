"""Portable entry point for formal or reduced EU57 workflows."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CORE = ROOT / "code/modeling/run_experiment.py"
PANEL = ROOT / "code/modeling/run_panel_sensitivity.py"
MANIFEST = ROOT / "data/input_manifest.tsv"
AUDIT = ROOT / "metadata/input_audit.tsv"
CLASS_MAP = ROOT / "metadata/disease_class_map.tsv"
RESOURCES = ROOT / "resources/annotation_sources"


def prepare_annotations(output: Path) -> None:
    destination = output / "annotation_sources"
    destination.mkdir(parents=True, exist_ok=True)
    for source in RESOURCES.glob("*.gmt"):
        shutil.copy2(source, destination / source.name)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the formal EU57 workflow or a reduced executability check."
    )
    parser.add_argument(
        "--mode",
        choices=("formal", "smoke"),
        required=True,
        help="Use 'formal' for manuscript settings and 'smoke' only for a quick check.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional output directory; defaults remain inside this package.",
    )
    args = parser.parse_args()

    # Packaged smoke artifacts under validation/ are immutable evidence from
    # release preparation. New runs always write under outputs/.
    default = ROOT / ("outputs/formal_results" if args.mode == "formal" else "outputs/smoke_results")
    output = (args.output or default).resolve()
    formal_default = (ROOT / "outputs/formal_results").resolve()
    if args.mode == "smoke" and output == formal_default:
        raise RuntimeError("Reduced workflow cannot write to the formal output directory")
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Output directory is not empty: {output}")

    output.mkdir(parents=True, exist_ok=True)
    prepare_annotations(output)
    command = [
        sys.executable,
        str(CORE),
        "--manifest", str(MANIFEST),
        "--audit", str(AUDIT),
        "--class-map", str(CLASS_MAP),
        "--output", str(output),
    ]
    if args.mode == "smoke":
        command.append("--quick")
    subprocess.run(command, cwd=ROOT, check=True)

    if args.mode == "formal":
        subprocess.run(
            [sys.executable, str(PANEL), "--output", str(output)],
            cwd=ROOT,
            check=True,
        )
    print(f"Completed {args.mode} workflow: {output}")


if __name__ == "__main__":
    main()
