"""Write a relative SHA-256 manifest for the complete package."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "PACKAGE_MANIFEST.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def main() -> None:
    rows = []
    for path in sorted(p for p in ROOT.rglob("*") if p.is_file() and p != OUTPUT):
        if "__pycache__" in path.parts or path.suffix == ".pyc":
            continue
        rows.append(
            {
                "relative_path": path.relative_to(ROOT).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    pd.DataFrame(rows).to_csv(OUTPUT, sep="\t", index=False)
    print(f"Wrote {len(rows)} entries to {OUTPUT}")


if __name__ == "__main__":
    main()
