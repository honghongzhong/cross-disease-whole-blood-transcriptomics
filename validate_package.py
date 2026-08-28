"""Static integrity and portability checks for the EU57 package."""

from __future__ import annotations

import hashlib
import json
import py_compile
import re
import tempfile
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
PUBLIC_TEXT_SUFFIXES = {".json", ".md", ".py", ".tsv", ".txt", ".yaml", ".yml"}
# Construct the detector from fragments so this source file does not contain a
# literal example that would trigger its own public-text scan.
ABSOLUTE_PATH = re.compile(
    r"(?:(?<![A-Za-z0-9])[A-Za-z]"
    + r":[\\/]"
    + r"|/"
    + r"(?:Users|home)/"
    + r"|\\\\"
    + r"\?\\)"
)


def label(path: Path) -> str:
    """Return a package-relative path for validation reports."""
    return path.relative_to(ROOT).as_posix()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def main() -> None:
    errors: list[str] = []
    manifest = pd.read_csv(ROOT / "data/input_manifest.tsv", sep="\t", dtype=str)
    if len(manifest) != 57 or manifest["EU_ID"].nunique() != 57:
        errors.append("Input manifest does not contain 57 unique contrasts")
    for row in manifest.itertuples(index=False):
        path = ROOT / "data" / row.DatasetPath
        if not path.is_file():
            errors.append(f"Missing standardized input: {label(path)}")
        elif sha256(path) != row.SHA256.upper():
            errors.append(f"Input hash mismatch: {label(path)}")

    process_manifest = pd.read_csv(
        ROOT / "metadata/dataset_process_code_manifest.tsv", sep="\t", dtype=str
    )
    if process_manifest["dataset_directory"].nunique() != 57:
        errors.append("Dataset processing records do not cover all 57 directories")
    for row in process_manifest.itertuples(index=False):
        path = ROOT / row.relative_path
        if not path.is_file() or sha256(path) != row.sha256.upper():
            errors.append(f"Dataset-process file missing or changed: {label(path)}")

    python_files = sorted((ROOT / "code").rglob("*.py"))
    # Compile into a temporary directory so validation does not add
    # __pycache__ files to the release tree.
    with tempfile.TemporaryDirectory(prefix="eu57_compile_") as compile_dir:
        compile_root = Path(compile_dir)
        for index, path in enumerate(python_files):
            try:
                py_compile.compile(
                    str(path),
                    cfile=str(compile_root / f"{index:03d}.pyc"),
                    doraise=True,
                )
            except py_compile.PyCompileError as exc:
                errors.append(f"Python compilation failed: {label(path)}: {exc.msg}")

    reference = ROOT / "reference_results"
    if (reference / "SMOKE_TEST_ONLY.txt").exists():
        errors.append("Reference formal results contain a reduced-workflow marker")
    audit = json.loads((reference / "reproducibility_audit.json").read_text(encoding="utf-8"))
    if audit.get("formal_output_has_smoke_marker") is not False:
        errors.append("Reference audit does not confirm formal/reduced separation")
    artifact_manifest = pd.read_csv(reference / "artifact_manifest.tsv", sep="\t", dtype=str)
    for row in artifact_manifest.itertuples(index=False):
        path = reference / Path(row.path)
        if not path.is_file() or sha256(path) != row.sha256.upper():
            errors.append(f"Reference result missing or changed: {label(path)}")

    packaged_smoke = ROOT / "validation/smoke_results"
    if not (packaged_smoke / "SMOKE_TEST_ONLY.txt").is_file():
        errors.append("Packaged smoke results are missing their non-reportable marker")

    # Public text artifacts must not expose a developer's local directory.
    absolute_path_files = []
    for path in sorted(p for p in ROOT.rglob("*") if p.is_file()):
        if path.suffix.lower() not in PUBLIC_TEXT_SUFFIXES or path.name == "PACKAGE_MANIFEST.tsv":
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        if ABSOLUTE_PATH.search(text):
            absolute_path_files.append(label(path))
    if absolute_path_files:
        errors.extend(f"Machine-specific absolute path found: {path}" for path in absolute_path_files)

    report = {
        "schema": "eu57_package_validation_v1",
        "input_contrasts": len(manifest),
        "dataset_process_directories": int(process_manifest["dataset_directory"].nunique()),
        "dataset_process_files": len(process_manifest),
        "python_files_compiled": len(python_files),
        "reference_artifacts_checked": len(artifact_manifest),
        "absolute_path_files": absolute_path_files,
        "formal_and_smoke_outputs_separated": not (reference / "SMOKE_TEST_ONLY.txt").exists(),
        "errors": errors,
        "passed": not errors,
    }
    output = ROOT / "validation/package_validation.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
