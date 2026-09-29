"""Record source provenance for the processing outputs."""

from pathlib import Path
from config import work_root
import json, hashlib, shutil, os, sys
import pandas as pd

ROOT = work_root()
OUT = ROOT / "outputs/stage1_three_effects_20260903_v1"
HERE = Path(__file__).parent


def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda: f.read(2**20), b""):
            h.update(b)
    return h.hexdigest()


prov = OUT / "provenance"
prov.mkdir(exist_ok=True)
additional_sources = [
    ROOT / "GPL570.annot.gz",
    ROOT
    / "phase1b_EU_coordination_v1_2026-08-23/01_candidates/GSE83632/GPL5175_matrix_probe_mapping.tsv",
]
additional_sources += list(
    (ROOT / "tmp/eu58_full_model_20260827/annotations").glob("*.annot.gz")
)
additional_sources += [
    ROOT / "tmp/eu58_full_model_20260827/annotations" / name
    for name in ["GPL8136.txt", "GPL16384.txt"]
]
if "--artifacts-only" in sys.argv:
    assert (prov / "source_file_manifest.tsv").exists()
    source_table = pd.read_csv(prov / "source_file_manifest.tsv", sep="\t")
    known = set(source_table.path)
    extra = [
        {"path": str(p), "bytes": p.stat().st_size, "sha256": sha(p)}
        for p in additional_sources
        if p.is_file() and str(p) not in known
    ]
    if extra:
        pd.concat([source_table, pd.DataFrame(extra)], ignore_index=True).sort_values(
            "path"
        ).to_csv(prov / "source_file_manifest.tsv", sep="\t", index=False)
    codes = [
        {"path": str(p), "sha256": sha(p), "bytes": p.stat().st_size}
        for p in HERE.rglob("*")
        if p.is_file() and "__pycache__" not in p.parts
    ]
    pd.DataFrame(codes).to_csv(prov / "code_manifest.tsv", sep="\t", index=False)
    art = []
    for p in sorted(OUT.rglob("*")):
        if (
            not p.is_file()
            or "preparation_revisions" in p.parts
            or p.name == "artifact_manifest.tsv"
            or p.suffix == ".log"
            or "before_" in p.name
        ):
            continue
        art.append(
            {
                "path": str(p.relative_to(OUT)),
                "bytes": p.stat().st_size,
                "sha256": sha(p),
            }
        )
    pd.DataFrame(art).to_csv(OUT / "artifact_manifest.tsv", sep="\t", index=False)
    print("ARTIFACT_PROVENANCE_COMPLETE", len(art), flush=True)
    sys.exit(0)
for src, name in [
    ("dataset/DATASET_MANIFEST.tsv", "original_dataset_manifest.tsv"),
    (
        "outputs/dataset_contrast_audit_20260903_v1/notion_project_snapshot.json",
        "notion_project_snapshot.json",
    ),
    (
        "outputs/dataset_contrast_audit_20260903_v1/contrast_inventory.tsv",
        "screening_inventory.tsv",
    ),
]:
    shutil.copy2(ROOT / src, prov / name)
infos = [
    json.loads(p.read_text(encoding="utf8"))
    for p in (OUT / "contrasts").glob("*/input.json")
]
paths = {
    ROOT / "tools/r461_bioc323_gse153315/library/org.Hs.eg.db/extdata/org.Hs.eg.sqlite"
}
paths.add(ROOT / "external/hgu133plus2.sqlite")
paths.update(p for p in additional_sources if p.is_file())
for info in infos:
    base = ROOT / "dataset" / info["directory"]
    paths.update(
        p
        for p in base.rglob("*")
        if p.is_file()
        and ("raw" in p.parts or "process" in p.parts)
        and "__pycache__" not in p.parts
        and "input_E" not in p.name
        and "input_U" not in p.name
        and "input_effects" not in p.name
    )
    for part in info["source"].split(";"):
        p = Path(part)
        if p.is_file():
            paths.add(p)
    for part in info.get("evidence_paths", "").split(";"):
        p = Path(part)
        if p.is_file():
            paths.add(p)
rows = []
cache = {}
for i, p in enumerate(sorted(paths)):
    st = p.stat()
    key = (st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns)
    if key not in cache:
        cache[key] = sha(p)
    rows.append({"path": str(p), "bytes": st.st_size, "sha256": cache[key]})
    if i % 50 == 0:
        print("SOURCE_HASH", i + 1, "/", len(paths), flush=True)
pd.DataFrame(rows).to_csv(prov / "source_file_manifest.tsv", sep="\t", index=False)
codes = [
    {"path": str(p), "sha256": sha(p), "bytes": p.stat().st_size}
    for p in HERE.rglob("*")
    if p.is_file() and "__pycache__" not in p.parts
]
pd.DataFrame(codes).to_csv(prov / "code_manifest.tsv", sep="\t", index=False)
if "--sources-only" in sys.argv:
    print("SOURCE_PROVENANCE_COMPLETE", len(rows), flush=True)
    sys.exit(0)
art = []
for p in sorted(OUT.rglob("*")):
    if (
        not p.is_file()
        or any(part in ["preparation_revisions"] for part in p.parts)
        or p.name == "artifact_manifest.tsv"
        or p.suffix == ".log"
        or "before_" in p.name
    ):
        continue
    art.append(
        {"path": str(p.relative_to(OUT)), "bytes": p.stat().st_size, "sha256": sha(p)}
    )
pd.DataFrame(art).to_csv(OUT / "artifact_manifest.tsv", sep="\t", index=False)
print(
    "PROVENANCE_COMPLETE", len(rows), "source files", len(art), "artifacts", flush=True
)
