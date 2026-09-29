"""Source-specific readers retained for the deposited study formats."""

from __future__ import annotations
import csv
import gzip
import io
import json
import os
import re
import sqlite3
from pathlib import Path
import pandas as pd

ROOT = Path(os.environ["BIOLOGY_WORK_ROOT"])
TMP = ROOT / "tmp/eu58_full_model_20260827"
ANN_DIR = TMP / "annotations"
ORG_DB = (
    ROOT / "tools/r461_bioc323_gse153315/library/org.Hs.eg.db/extdata/org.Hs.eg.sqlite"
)


def read_geo_series_matrix(path: Path) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    opener = gzip.open if path.suffix == ".gz" else open
    meta: dict[str, list[str]] = {}
    table_lines: list[str] = []
    in_table = False
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                table_lines.append(line)
            elif line.startswith("!Sample_"):
                parts = next(csv.reader([line.rstrip("\r\n")], delimiter="\t"))
                meta.setdefault(parts[0], []).append([x.strip('"') for x in parts[1:]])
    if not table_lines:
        raise ValueError(f"No series matrix table: {path}")
    frame = pd.read_csv(io.StringIO("".join(table_lines)), sep="\t", quotechar='"')
    frame = frame.rename(columns={frame.columns[0]: "feature_id"})
    return (frame, meta)


def read_annotation(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
    lines = []
    in_table = False
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("!platform_table_begin", "!annotation_table_begin")):
                in_table = True
                continue
            if line.startswith(("!platform_table_end", "!annotation_table_end")):
                break
            if in_table:
                lines.append(line)
    if not lines:
        raise ValueError(f"No annotation table: {path}")
    return pd.read_csv(
        io.StringIO("".join(lines)), sep="\t", dtype=str, keep_default_na=False
    )


def load_org_maps() -> dict[str, dict[str, str]]:
    con = sqlite3.connect(ORG_DB)
    maps: dict[str, dict[str, str]] = {}
    for name, table, col in (
        ("symbol", "gene_info", "symbol"),
        ("ensembl", "ensembl", "ensembl_id"),
        ("ensembl_trans", "ensembl_trans", "trans_id"),
        ("refseq", "refseq", "accession"),
    ):
        rows = con.execute(
            f"select t.{col}, g.gene_id from {table} t join genes g on t._id=g._id where t.{col} is not null"
        ).fetchall()
        bucket: dict[str, set[str]] = {}
        for key, gene in rows:
            bucket.setdefault(str(key).upper(), set()).add(str(gene))
        maps[name] = {k: next(iter(v)) for k, v in bucket.items() if len(v) == 1}
    symbols = con.execute(
        "select g.gene_id, i.symbol from genes g join gene_info i on g._id=i._id"
    ).fetchall()
    maps["entrez_symbol"] = {str(g): str(s) for g, s in symbols}
    maps["valid_entrez"] = {str(g): str(g) for g, _ in symbols}
    con.close()
    return maps


_cache_path = os.environ.get("BIOLOGY_ANNOTATION_CACHE")
if _cache_path:
    with gzip.open(_cache_path, "rt", encoding="utf8") as _cache_file:
        _annotation_cache = json.load(_cache_file)
    ORG = _annotation_cache["ORG"]
else:
    _annotation_cache = None
    ORG = load_org_maps()


def annotation_probe_map(platform: str) -> dict[str, str]:
    candidates = {
        "GPL570": ROOT / "GPL570.annot.gz",
        "GPL2986": ANN_DIR / "GPL2986.annot.gz",
        "GPL6947": ANN_DIR / "GPL6947.annot.gz",
        "GPL13158": ANN_DIR / "GPL13158.annot.gz",
        "GPL6480": ANN_DIR / "GPL6480.annot.gz",
        "GPL10558": ANN_DIR / "GPL10558.annot.gz",
        "GPL8136": ANN_DIR / "GPL8136.txt",
        "GPL16384": ANN_DIR / "GPL16384.txt",
    }
    if platform == "GPL5175":
        p = (
            ROOT
            / "phase1b_EU_coordination_v1_2026-08-23/01_candidates/GSE83632/GPL5175_matrix_probe_mapping.tsv"
        )
        d = pd.read_csv(
            p,
            sep="\t",
            dtype=str,
            keep_default_na=False,
            usecols=["ID", "gene_assignment"],
        )
        out = {}
        for probe, assignment in zip(d["ID"], d["gene_assignment"]):
            ids = set()
            for part in str(assignment).split(" /// "):
                tokens = [token.strip() for token in part.split(" // ")]
                if tokens and re.fullmatch("\\d+", tokens[-1]):
                    ids.add(tokens[-1])
            if len(ids) == 1:
                out[str(probe)] = next(iter(ids))
        return out
    path = candidates.get(platform)
    if path is None or not path.exists():
        return {}
    d = read_annotation(path)
    id_col = "ID" if "ID" in d.columns else d.columns[0]
    gene_col = next(
        (
            c
            for c in d.columns
            if c.lower()
            in {"gene id", "entrez gene", "entrez_gene_id", "entrezid", "orf"}
        ),
        None,
    )
    symbol_col = next(
        (
            c
            for c in d.columns
            if c.lower()
            in {"gene symbol", "symbol", "mirna_id", "transcript id(array design)"}
        ),
        None,
    )
    out: dict[str, str] = {}
    for _, row in d.iterrows():
        probe = str(row[id_col])
        ids: set[str] = set()
        if gene_col:
            ids = set(re.findall("(?<!\\d)(\\d+)(?!\\d)", str(row[gene_col])))
        if not ids and symbol_col:
            symbols = re.split("///|//|,|;|\\s+", str(row[symbol_col]))
            ids = {ORG["symbol"].get(x.strip().upper(), "") for x in symbols}
            ids.discard("")
        if len(ids) == 1:
            out[probe] = next(iter(ids))
    return out


PROBE_MAPS = (
    _annotation_cache["PROBE_MAPS"]
    if _annotation_cache is not None
    else {
        p: annotation_probe_map(p)
        for p in (
            "GPL570",
            "GPL2986",
            "GPL6947",
            "GPL5175",
            "GPL13158",
            "GPL6480",
            "GPL10558",
            "GPL8136",
        )
    }
)


def map_feature(
    value: object, feature_type: str, platform: str | None = None
) -> str | None:
    raw = str(value).strip().strip('"')
    if platform:
        hit = PROBE_MAPS.get(platform, {}).get(raw)
        if hit:
            return hit
    if feature_type == "entrez":
        key = re.sub("\\.0$", "", raw)
        return key if key in ORG["valid_entrez"] else None
    if feature_type == "ensembl":
        key = raw.split(".")[0].upper()
        return ORG["ensembl"].get(key)
    if feature_type == "ensembl_trans":
        key = raw.split(".")[0].upper()
        return ORG["ensembl_trans"].get(key)
    if feature_type == "refseq":
        key = raw.split(".")[0].split(">")[0].upper()
        return ORG["refseq"].get(key)
    if feature_type == "symbol":
        key = raw.split(">")[-1].strip().upper()
        return ORG["symbol"].get(key)
    return None


def capture_selection(
    contrast_id,
    frame,
    cases,
    controls,
    feature_col,
    feature_type,
    kind,
    source,
    platform=None,
    symbol_col=None,
    feature_space=None,
):
    """Return a selected source matrix for the shared effect estimator."""
    cases = list(cases)
    controls = list(controls)
    return {
        "x": frame.set_index(feature_col)[cases + controls],
        "cases": cases,
        "controls": controls,
        "kind": kind,
        "source": source,
        "platform": platform,
        "feature_type": feature_type,
    }


def _declared_source_path(row: pd.Series, filename: str) -> Path:
    """Resolve an external input by its declared filename, not an audit folder."""
    matches = [
        ROOT / item["relative"]
        for item in row.get("source_inputs", [])
        if Path(item["relative"]).name == filename
    ]
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one declared source named {filename}; found {len(matches)}"
        )
    path = matches[0]
    if not path.is_file():
        raise FileNotFoundError(f"Required source input is missing: {filename}")
    return path


def load_declared_source(row: pd.Series, sample_selection: Path) -> dict | None:
    """Load exceptional source formats using the comparison's sample contract."""
    routes = {
        "GSE83563_EVD_vs_HC": (
            "GSE83563_sample_expression.txt.gz",
            "geneID",
            "ensembl",
            "normalized",
            None,
        ),
        "GSE236442_EH_vs_HC": (
            "GSE236442_hypertension_EH_vs_healthy_expression_matrix.txt",
            "gene_id",
            "ensembl",
            "normalized",
            None,
        ),
        "GSE269497_AIT_vs_HC": (
            "GSE269497_genes_FPKM_expression.txt.gz",
            "Gene ID",
            "entrez",
            "normalized",
            "Gene Symbol",
        ),
        "GSE262659_pediatric_SCD_vs_HC": (
            "GSE262659_SALMON_Counts.txt.gz",
            "gene_id",
            "ensembl",
            "counts",
            None,
        ),
        "GSE244401_SCA_T1_vs_HC": (
            "GSE244401_processed_data.txt.gz",
            "ID",
            "ensembl",
            "counts",
            None,
        ),
        "GSE51405_CVID_vs_HC": (
            "series_matrix.txt.gz",
            "feature_id",
            "probe",
            "normalized",
            None,
        ),
        "GSE44969_CGD_vs_HC": (
            "GSE44969_non-normalized.txt.gz",
            None,
            "probe",
            "normalized",
            None,
        ),
    }
    contrast_id = str(row["contrast_id"])
    route = routes.get(contrast_id)
    if route is None:
        return None

    filename, feature_col, feature_type, kind, symbol_col = route
    path = _declared_source_path(row, filename)
    selection = pd.read_csv(sample_selection, sep="\t", dtype=str)
    cases = selection.loc[selection.group.eq("case"), "sample_id"].tolist()
    controls = selection.loc[selection.group.eq("control"), "sample_id"].tolist()
    if contrast_id == "GSE83563_EVD_vs_HC":
        frame = pd.read_csv(path, sep="\t")
    elif contrast_id == "GSE51405_CVID_vs_HC":
        frame, _ = read_geo_series_matrix(path)
    elif contrast_id == "GSE262659_pediatric_SCD_vs_HC":
        frame = pd.read_csv(path, index_col=0).reset_index(names="gene_id")
    elif contrast_id == "GSE44969_CGD_vs_HC":
        frame = pd.read_csv(path, sep="\t")
        keep = [
            frame.columns[0],
            *[
                column
                for column in frame.columns[1:]
                if "detection" not in str(column).lower()
                and "pval" not in str(column).lower()
                and not str(column).lower().startswith("unnamed")
            ],
        ]
        frame = frame[keep]
        frame = frame.rename(
            columns={
                column: str(column).replace(".AVG_Signal", "")
                for column in frame.columns[1:]
            }
        )
    else:
        frame = pd.read_csv(path, sep="\t")

    available = set(map(str, frame.columns))
    missing = (set(cases) | set(controls)) - available
    if missing:
        raise ValueError(
            f"Declared samples are absent from {filename}: {sorted(missing)[:5]}"
        )
    return capture_selection(
        contrast_id,
        frame,
        cases,
        controls,
        feature_col or frame.columns[0],
        feature_type,
        kind,
        str(path),
        "GPL6947" if feature_type == "probe" else None,
        symbol_col,
    )
