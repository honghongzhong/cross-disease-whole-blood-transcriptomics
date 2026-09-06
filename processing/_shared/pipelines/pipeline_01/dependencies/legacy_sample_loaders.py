from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import re
import sqlite3
import sys
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

sys.path.append(r"C:\Users\peter\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\Lib\site-packages")

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(r"E:\Biology")
BOOK = ROOT / "outputs/eu58_model_base_selection_reduced_gse100150_20260827/eu58_model_base_selection_reduced_gse100150.xlsx"
TMP = ROOT / "tmp/eu58_full_model_20260827"
OUT = ROOT / "outputs/eu58_full_disease_model_20260827"
EU_OUT = OUT / "eu"
ANN_DIR = TMP / "annotations"
ORG_DB = ROOT / "tools/r461_bioc323_gse153315/library/org.Hs.eg.db/extdata/org.Hs.eg.sqlite"
SEED = 20260827


@dataclass
class EUResult:
    eu_id: str
    genes: np.ndarray
    symbols: np.ndarray
    effect: np.ndarray
    uncertainty: np.ndarray
    p_value: np.ndarray
    adj_p_value: np.ndarray
    n_case: int
    n_control: int
    method: str
    source: str
    feature_space: str = "mRNA_Entrez"


def read_table(path: Path, **kwargs) -> pd.DataFrame:
    return pd.read_csv(path, **kwargs)


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
    return frame, meta


def read_family_soft_values(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
    values: dict[str, dict[str, float]] = {}
    current = None
    in_table = False
    id_idx = value_idx = None
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                current = line.split("=", 1)[1].strip()
            elif line.startswith("!sample_table_begin"):
                in_table = True
                id_idx = value_idx = None
            elif line.startswith("!sample_table_end"):
                in_table = False
            elif in_table and current:
                parts = line.split("\t")
                if id_idx is None:
                    header = [x.strip('"') for x in parts]
                    id_idx = header.index("ID_REF") if "ID_REF" in header else 0
                    value_idx = header.index("VALUE") if "VALUE" in header else 1
                elif len(parts) > max(id_idx, value_idx):
                    try:
                        values.setdefault(parts[id_idx].strip('"'), {})[current] = float(parts[value_idx].strip('"'))
                    except ValueError:
                        pass
    if not values:
        raise ValueError(f"No numeric sample tables: {path}")
    frame = pd.DataFrame.from_dict(values, orient="index").reset_index(names="feature_id")
    return frame


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
    return pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str, keep_default_na=False)


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
        p = ROOT / "phase1b_EU_coordination_v1_2026-08-23/01_candidates/GSE83632/GPL5175_matrix_probe_mapping.tsv"
        d = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False, usecols=["ID", "gene_assignment"])
        out = {}
        for probe, assignment in zip(d["ID"], d["gene_assignment"]):
            ids = set()
            for part in str(assignment).split(" /// "):
                tokens = [token.strip() for token in part.split(" // ")]
                if tokens and re.fullmatch(r"\d+", tokens[-1]):
                    ids.add(tokens[-1])
            if len(ids) == 1:
                out[str(probe)] = next(iter(ids))
        return out
    path = candidates.get(platform)
    if path is None or not path.exists():
        return {}
    d = read_annotation(path)
    id_col = "ID" if "ID" in d.columns else d.columns[0]
    gene_col = next((c for c in d.columns if c.lower() in {"gene id", "entrez gene", "entrez_gene_id", "entrezid", "orf"}), None)
    symbol_col = next((c for c in d.columns if c.lower() in {"gene symbol", "symbol", "mirna_id", "transcript id(array design)"}), None)
    out: dict[str, str] = {}
    for _, row in d.iterrows():
        probe = str(row[id_col])
        ids: set[str] = set()
        if gene_col:
            ids = set(re.findall(r"(?<!\d)(\d+)(?!\d)", str(row[gene_col])))
        if not ids and symbol_col:
            symbols = re.split(r"///|//|,|;|\s+", str(row[symbol_col]))
            ids = {ORG["symbol"].get(x.strip().upper(), "") for x in symbols}
            ids.discard("")
        if len(ids) == 1:
            out[probe] = next(iter(ids))
    return out


PROBE_MAPS = {p: annotation_probe_map(p) for p in ("GPL570", "GPL2986", "GPL6947", "GPL5175", "GPL13158", "GPL6480", "GPL10558", "GPL8136")}


def map_feature(value: object, feature_type: str, platform: str | None = None) -> str | None:
    raw = str(value).strip().strip('"')
    if platform:
        hit = PROBE_MAPS.get(platform, {}).get(raw)
        if hit:
            return hit
    if feature_type == "entrez":
        key = re.sub(r"\.0$", "", raw)
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


def bh_adjust(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adj = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(adj)
    out[order] = np.minimum(adj, 1.0)
    return out


def transform_expression(x: np.ndarray, data_kind: str) -> tuple[np.ndarray, str]:
    x = np.asarray(x, dtype=float)
    if data_kind == "counts":
        x = np.maximum(x, 0)
        lib = np.sum(x, axis=0)
        lib[lib <= 0] = 1
        return np.log2((x / lib[None, :]) * 1e6 + 0.5), "log2_CPM_plus_0.5"
    finite = x[np.isfinite(x)]
    if finite.size == 0:
        return x, "none"
    if np.nanmin(finite) >= 0 and np.nanpercentile(finite, 95) > 100:
        return np.log2(x + 1), "log2_x_plus_1"
    return x, "as_deposited_scale"


def build_eu(
    eu_id: str,
    frame: pd.DataFrame,
    case_cols: list[str],
    control_cols: list[str],
    feature_col: str,
    feature_type: str,
    data_kind: str,
    source: str,
    platform: str | None = None,
    symbol_col: str | None = None,
    feature_space: str = "mRNA_Entrez",
) -> EUResult:
    if len(case_cols) < 2 or len(control_cols) < 2:
        raise ValueError(f"insufficient groups after manifest mapping: case={len(case_cols)} control={len(control_cols)}")
    cols = case_cols + control_cols
    missing = [c for c in cols if c not in frame.columns]
    if missing:
        raise ValueError(f"{eu_id}: missing matrix columns {missing[:8]}")
    raw = frame[cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    features = frame[feature_col].astype(str).tolist()
    if feature_space == "miRNA_probe":
        mapped = np.array(features, dtype=object)
        symbols = np.array(features, dtype=object)
    else:
        mapped = np.array([map_feature(v, feature_type, platform) for v in features], dtype=object)
        if symbol_col and symbol_col in frame.columns:
            symbols = frame[symbol_col].astype(str).to_numpy(object)
        else:
            symbols = np.array([ORG["entrez_symbol"].get(str(g), "") if g else "" for g in mapped], dtype=object)
    keep = np.array([m is not None and str(m) != "" for m in mapped]) & np.isfinite(raw).any(axis=1)
    raw, mapped, symbols = raw[keep], mapped[keep], symbols[keep]
    if raw.shape[0] == 0:
        raise ValueError(f"{eu_id}: no mapped numeric features")
    # Outcome-blind collapse: retain the feature with highest all-sample mean per gene.
    means = np.nanmean(raw, axis=1)
    chosen: dict[str, int] = {}
    for i, gene in enumerate(mapped.astype(str)):
        if gene not in chosen or means[i] > means[chosen[gene]]:
            chosen[gene] = i
    idx = np.array(sorted(chosen.values()), dtype=int)
    raw, mapped, symbols = raw[idx], mapped[idx].astype(str), symbols[idx].astype(str)
    x, transform = transform_expression(raw, data_kind)
    n_case, n_control = len(case_cols), len(control_cols)
    case = x[:, :n_case]
    control = x[:, n_case:]
    effect = np.nanmean(case, axis=1) - np.nanmean(control, axis=1)
    vc = np.nanvar(case, axis=1, ddof=1)
    vh = np.nanvar(control, axis=1, ddof=1)
    se = np.sqrt(vc / n_case + vh / n_control)
    positive = se[np.isfinite(se) & (se > 0)]
    floor = np.nanpercentile(positive, 1) if positive.size else 1e-6
    se[~np.isfinite(se) | (se <= 0)] = floor
    t = effect / se
    num = (vc / n_case + vh / n_control) ** 2
    den = (vc / n_case) ** 2 / max(n_case - 1, 1) + (vh / n_control) ** 2 / max(n_control - 1, 1)
    df = np.divide(num, den, out=np.full_like(num, n_case + n_control - 2, dtype=float), where=den > 0)
    p = 2 * stats.t.sf(np.abs(t), df=np.maximum(df, 1))
    p[~np.isfinite(p)] = 1
    adj = bh_adjust(p)
    finite = np.isfinite(effect) & np.isfinite(se) & (se > 0)
    return EUResult(
        eu_id=eu_id,
        genes=mapped[finite],
        symbols=symbols[finite],
        effect=effect[finite],
        uncertainty=se[finite],
        p_value=p[finite],
        adj_p_value=adj[finite],
        n_case=n_case,
        n_control=n_control,
        method=f"unified exploratory Welch E/U; {transform}; outcome-blind highest-mean feature collapse",
        source=source,
        feature_space=feature_space,
    )


def manifest_groups(path: Path, key: str, label: Callable[[pd.Series], str | None]) -> dict[str, str]:
    d = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    out = {}
    for _, row in d.iterrows():
        group = label(row)
        if group:
            out[str(row[key]).strip()] = group
    return out


def cols_from_map(columns: list[object], mapping: dict[str, str]) -> tuple[list[str], list[str]]:
    cases, controls = [], []
    for col in map(str, columns):
        key = col.strip().strip('"')
        group = mapping.get(key)
        if group == "case":
            cases.append(col)
        elif group == "control":
            controls.append(col)
    return cases, controls


def extract_signal_columns(frame: pd.DataFrame) -> pd.DataFrame:
    keep = [frame.columns[0]]
    for i, col in enumerate(frame.columns[1:], start=1):
        low = str(col).lower()
        if "detection" not in low and "pval" not in low and not low.startswith("unnamed"):
            keep.append(col)
    return frame[keep]


def load_existing_eu(row: pd.Series) -> EUResult:
    path = Path(str(row["证据路径"]).replace("E:/", "E:\\"))
    d = pd.read_csv(path, sep="\t", compression="infer", dtype=str)
    lower = {c.lower().strip(): c for c in d.columns}
    gene_col = lower.get("gene_id") or lower.get("entrez") or d.columns[0]
    if "effect" in lower and ("se" in lower or "u" in lower):
        effect_col = lower["effect"]
        u_col = lower.get("se") or lower["u"]
        p_col = lower.get("p_value")
        q_col = lower.get("adjusted_p_value")
        effect = pd.to_numeric(d[effect_col], errors="coerce").to_numpy(float)
        u = pd.to_numeric(d[u_col], errors="coerce").to_numpy(float)
        p = pd.to_numeric(d[p_col], errors="coerce").fillna(1).to_numpy(float) if p_col else np.ones(len(d))
        q = pd.to_numeric(d[q_col], errors="coerce").fillna(1).to_numpy(float) if q_col else bh_adjust(p)
    else:
        effect_col = str(row["EU_ID"]) if str(row["EU_ID"]) in d.columns else d.columns[-1]
        effect = pd.to_numeric(d[effect_col], errors="coerce").to_numpy(float)
        name = path.name
        if "_E_locked_" in name:
            u_path = path.with_name(name.replace("_E_locked_", "_U_locked_"))
        elif "_E_" in name:
            u_path = path.with_name(name.replace("_E_", "_U_", 1))
        else:
            raise ValueError(f"Cannot derive U path for {path}")
        ud = pd.read_csv(u_path, sep="\t", compression="infer", dtype=str)
        u_col = str(row["EU_ID"]) if str(row["EU_ID"]) in ud.columns else ud.columns[-1]
        u = pd.to_numeric(ud[u_col], errors="coerce").to_numpy(float)
        p = np.ones(len(d))
        q = np.ones(len(d))
    raw_genes = d[gene_col].astype(str).str.replace(r"\.0$", "", regex=True).to_numpy(object)
    genes = np.array([
        str(g) if str(g) in ORG["valid_entrez"]
        else ORG["symbol"].get(str(g).upper(), "")
        for g in raw_genes
    ], dtype=object)
    symbols = np.array([ORG["entrez_symbol"].get(str(g), "") for g in genes], dtype=object)
    keep = np.array([str(g) in ORG["valid_entrez"] for g in genes]) & np.isfinite(effect) & np.isfinite(u) & (u > 0)
    if not np.any(keep):
        raise ValueError(f"existing E/U contains no mappable numeric features: {path}")
    n_case = int(pd.to_numeric(d[lower["n_case"]], errors="coerce").dropna().iloc[0]) if "n_case" in lower else int(row["病例N"])
    n_control = int(pd.to_numeric(d[lower["n_control"]], errors="coerce").dropna().iloc[0]) if "n_control" in lower else int(row["对照N"])
    return EUResult(str(row["EU_ID"]), genes[keep].astype(str), symbols[keep], effect[keep], u[keep], p[keep], q[keep], n_case, n_control, "existing locked E/U", str(path))


def load_missing_eu(row: pd.Series) -> EUResult:
    i = int(row["序号"])
    eu = str(row["EU_ID"])
    base = Path(str(row["证据路径"]).replace("E:/", "E:\\")).parent
    manifest = base / "sample_manifest.tsv"
    if not manifest.exists():
        manifest = base / "selection_manifest.tsv"

    def finish(frame, cases, controls, feature_col, feature_type, kind, source, platform=None, symbol_col=None, feature_space="mRNA_Entrez"):
        return build_eu(eu, frame, cases, controls, feature_col, feature_type, kind, str(source), platform, symbol_col, feature_space)

    if i == 1:
        frame = read_family_soft_values(base / "raw_evidence/GSE41055_family.soft.gz")
        mp = manifest_groups(manifest, "gsm", lambda r: "case" if r["contrast_role"] == "case" and r["inclusion_status"] == "primary_included" else ("control" if r["contrast_role"] == "control" and r["inclusion_status"] == "primary_included" else None))
        cases, controls = cols_from_map(list(frame.columns[1:]), mp)
        return finish(frame, cases, controls, "feature_id", "probe", "normalized", base / "raw_evidence/GSE41055_family.soft.gz", "GPL5175")
    if i == 3:
        p = base / "raw_evidence/GSE158592_norm_SIG_Mar18_070920.txt.gz"
        frame = pd.read_csv(p, sep=r"\s+", quotechar='"')
        mp = manifest_groups(manifest, "matrix_id", lambda r: "case" if r["group"] == "influenza" else ("control" if r["group"] == "healthy_control" else None))
        cases, controls = cols_from_map(list(frame.columns), mp)
        return finish(frame, cases, controls, "entrezgene_id", "entrez", "normalized", p, symbol_col="gene_symbol")
    if i == 6:
        p = ROOT / "phase1b_EU_coordination_v1_2026-08-23/02_worker_packages/worker_gse83563_20260824/GSE83563_sample_expression.txt.gz"
        frame = pd.read_csv(p, sep="\t")
        cases = [c for c in frame.columns if re.fullmatch(r"L\d+", str(c))]
        controls = [c for c in frame.columns if re.fullmatch(r"C\d+", str(c))]
        return finish(frame, cases, controls, "geneID", "ensembl", "normalized", p)
    if i in {8, 11, 18, 20, 27, 36, 51, 53, 56}:
        series_cfg = {
            8: (base / "GSE93272_series_matrix.txt.gz", "GPL570", "gsm", lambda r: "case" if r["primary"] == "yes" and r["disease_state"] == "RA" else ("control" if r["primary"] == "yes" and "healthy" in r["disease_state"].lower() else None)),
            11: (base / "raw_evidence/GSE25101_series_matrix.txt.gz", "GPL6947", "gsm", lambda r: "case" if r["group"] == "case" else ("control" if r["group"] == "control" else None)),
            18: (base / "GSE11545_series_matrix.txt.gz", "GPL2986", "gsm", lambda r: "case" if r["group"] == "breast_cancer" else ("control" if r["group"] == "control" else None)),
            20: (base / "GSE83632_series_matrix.txt.gz", "GPL5175", "accession", lambda r: "case" if r["group"] == "case" else "control"),
            27: (base / "GSE109597_series_matrix.txt.gz", "GPL570", "geo_accession", lambda r: "case" if r["case_control"] == "case" else "control"),
            36: (base / "GSE69683_series_matrix.txt.gz", "GPL13158", "gsm", lambda r: "case" if r["group"] == "moderate_asthma" else ("control" if r["group"] == "control" else None)),
            51: (base / "GSE14655_series_matrix.txt.gz", "GPL8136", "gsm", lambda r: "case" if r["group"] == "CR" else ("control" if r["group"] == "HC" else None)),
            53: (base / "GSE22229_series_matrix.txt.gz", "GPL570", "gsm", lambda r: "case" if r["group"] == "case" and r["include_primary"] == "yes" else ("control" if r["group"] == "control" and r["include_primary"] == "yes" else None)),
            56: (base / "GSE19743_series_matrix.txt.gz", "GPL570", "sample_id", lambda r: "case" if r["group"] == "burn" else ("control" if r["group"] == "healthy_control" else None)),
        }
        p, platform, key, labeller = series_cfg[i]
        frame, _ = read_geo_series_matrix(p)
        mp = manifest_groups(manifest, key, labeller)
        cases, controls = cols_from_map(list(frame.columns[1:]), mp)
        return finish(frame, cases, controls, "feature_id", "probe", "normalized", p, platform)
    if i == 10:
        p = base / "raw_evidence/GSE224849_processed_data.csv.gz"
        frame = pd.read_csv(p)
        mp = manifest_groups(manifest, "title", lambda r: "case" if r["formal_inclusion"] == "yes" and r["formal_group"] == "PreSSc" else ("control" if r["formal_inclusion"] == "yes" and r["formal_group"] == "healthy_control" else None))
        cases, controls = cols_from_map(list(frame.columns), mp)
        return finish(frame, cases, controls, "ensembl_id", "ensembl", "normalized", p, symbol_col="gene_symbol")
    if i in {19, 43, 44}:
        p = base / "raw_evidence/GSE26049_series_matrix.txt.gz"
        frame, meta = read_geo_series_matrix(p)
        gsm = meta["!Sample_geo_accession"][0]
        titles = meta["!Sample_title"][0]
        target = {19: "RMA PMF_", 43: "RMA ET_", 44: "RMA PV_"}[i]
        mp = {g: ("case" if t.startswith(target) else ("control" if t.startswith("RMA Control_") else None)) for g, t in zip(gsm, titles)}
        cases, controls = cols_from_map(list(frame.columns[1:]), mp)
        return finish(frame, cases, controls, "feature_id", "probe", "normalized", p, "GPL570")
    if i == 17:
        p = ROOT / "phase1b_EU_coordination_v1_2026-08-23/02_worker_packages/worker_gse198048_20260824/GSE198048_Normalized_mRNA_matrix.txt.gz"
        frame = pd.read_csv(p, sep="\t", index_col=0).reset_index(names="feature_id")
        mp = manifest_groups(manifest, "sample_id", lambda r: "case" if r["primary_role"] == "case" else "control")
        cases, controls = cols_from_map(list(frame.columns[1:]), mp)
        return finish(frame, cases, controls, "feature_id", "ensembl", "normalized", p)
    if i == 21:
        p = ROOT / "phase1b_EU_coordination_v1_2026-08-23/02_worker_packages/worker_gse244827_20260824/GSE244827_CHAVArawcounts.txt.gz"
        frame = pd.read_csv(p, sep="\t")
        mp = manifest_groups(manifest, "sample_id", lambda r: "case" if r["primary_role"] == "case" else ("control" if r["primary_role"] == "control" else None))
        cases, controls = cols_from_map(list(frame.columns[6:]), mp)
        return finish(frame, cases, controls, "Geneid", "ensembl", "counts", p)
    if i == 24:
        p = base / "raw_evidence/GSE68004_family.soft.gz"
        frame = read_family_soft_values(p)
        m = pd.read_csv(manifest, sep="\t", dtype=str)
        r = m[m["contrast"] == "complete_KD_vs_healthy"].iloc[0]
        cases = [x for x in str(r["case_gsm_list"]).split(";") if x in frame.columns]
        controls = [x for x in str(r["control_gsm_list"]).split(";") if x in frame.columns]
        return finish(frame, cases, controls, "feature_id", "probe", "normalized", p, "GPL10558")
    if i == 25:
        p = base / "GSE236442_hypertension_EH_vs_healthy_expression_matrix.txt"
        frame = pd.read_csv(p, sep="\t")
        mp = manifest_groups(manifest, "sample_id", lambda r: "case" if r["group"] == "essential_hypertension" else "control")
        cases, controls = cols_from_map(list(frame.columns[1:]), mp)
        return finish(frame, cases, controls, "gene_id", "ensembl", "normalized", p)
    if i == 26:
        p = base / "raw_evidence/GSE280402_counts.tsv.gz"
        frame = pd.read_csv(p, sep="\t")
        mp = manifest_groups(manifest, "matrix_id", lambda r: "case" if r["group"] == "T2D" and r["mapping_status"] == "PASS" else ("control" if r["group"] == "healthy_control" and r["mapping_status"] == "PASS" else None))
        cases, controls = cols_from_map(list(frame.columns), mp)
        return finish(frame, cases, controls, "ENSEMBL", "ensembl_trans", "counts", p, symbol_col="SYMBOL")
    if i == 28:
        p = base / "GSE145412_rnaseq.txt.gz"
        frame = pd.read_csv(p, sep="\t")
        mp = manifest_groups(manifest, "sample_id", lambda r: "case" if r["group"] == "DL" else ("control" if r["group"] == "HL" else None))
        cases, controls = cols_from_map(list(frame.columns), mp)
        return finish(frame, cases, controls, "EntrezID", "entrez", "normalized", p, symbol_col="External_gene_name")
    if i == 29:
        p = base / "raw_evidence/GSE269497_genes_FPKM_expression.txt.gz"
        frame = pd.read_csv(p, sep="\t")
        cases = [c for c in frame.columns if str(c).startswith("AIT-")]
        controls = [c for c in frame.columns if str(c).startswith("Con-")]
        return finish(frame, cases, controls, "Gene ID", "entrez", "normalized", p, symbol_col="Gene Symbol")
    if i == 30:
        p = base / "GSE176153_hypothyroidism_vs_healthy_expression_matrix.txt"
        frame = pd.read_csv(p, sep="\t")
        mp = manifest_groups(manifest, "sample_id", lambda r: "case" if r["group"] == "hypothyroidism" else "control")
        cases, controls = cols_from_map(list(frame.columns[1:]), mp)
        return finish(frame, cases, controls, "transcript", "refseq", "counts", p)
    if i == 31:
        p = base / "raw_evidence/GSE249477_raw_count_normalize_04-10-2025.csv.gz"
        frame = pd.read_csv(p, low_memory=False)
        mp = manifest_groups(manifest, "subject_id", lambda r: "case" if r["include_primary_ad_vs_healthy"] == "YES" and r["group"] == "AD" else ("control" if r["include_primary_ad_vs_healthy"] == "YES" and r["group"] == "healthy_control" else None))
        total_cols = [c for c in frame.columns if str(c).endswith(" - Total counts")]
        renamed = {c: str(c).split(" (GE)")[0] for c in total_cols}
        frame = frame[["Identifier", "Database symbol"] + total_cols].rename(columns=renamed)
        cases, controls = cols_from_map(list(frame.columns), mp)
        return finish(frame, cases, controls, "Identifier", "ensembl", "counts", p, symbol_col="Database symbol")
    if i == 33:
        p = ROOT / "phase1b_EU_coordination_v1_2026-08-23/02_worker_packages/worker_gse234297_20260824/GSE234297_gene_raw_counts.txt.gz"
        frame = pd.read_csv(p, sep="\t")
        cases = [c for c in frame.columns if str(c).startswith("Case")]
        controls = [c for c in frame.columns if str(c).startswith("Control")]
        return finish(frame, cases, controls, "EntrezGeneID", "entrez", "counts", p)
    if i == 34:
        p = base / "raw_evidence/GSE186505_fpkm.txt"
        frame = pd.read_csv(p, sep="\t")
        cases = [c for c in frame.columns if str(c).startswith("TN-")]
        controls = [c for c in frame.columns if str(c).startswith("Control-")]
        return finish(frame, cases, controls, "GeneID", "symbol", "normalized", p)
    if i == 37:
        p = base / "raw_evidence/GSE33566_series_matrix.txt.gz"
        frame, meta = read_geo_series_matrix(p)
        gsm = meta["!Sample_geo_accession"][0]
        titles = meta["!Sample_title"][0]
        mp = {g: ("case" if t.startswith("IPF_PBB") else ("control" if t.startswith("N_PBB") else None)) for g, t in zip(gsm, titles)}
        cases, controls = cols_from_map(list(frame.columns[1:]), mp)
        return finish(frame, cases, controls, "feature_id", "probe", "normalized", p, "GPL6480")
    if i == 38:
        p = ROOT / "phase1b_EU_coordination_v1_2026-08-23/02_worker_packages/worker_gse231691_20260824/GSE231691_count.tsv.gz"
        with gzip.open(p, "rt", encoding="utf-8") as f:
            headers = f.readline().rstrip().split("\t")
        frame = pd.read_csv(p, sep="\t", skiprows=1, header=None, names=["gene_id"] + headers)
        cases = [c for c in headers if c.startswith("SSc_")]
        controls = [c for c in headers if c.startswith("Normal_")]
        return finish(frame, cases, controls, "gene_id", "entrez", "counts", p)
    if i == 40:
        p = base / "raw_evidence/GSE262659_SALMON_Counts.txt.gz"
        frame = pd.read_csv(p, index_col=0).reset_index(names="gene_id")
        cases = [c for c in frame.columns if str(c).startswith("PSCD")]
        controls = [c for c in frame.columns if str(c).startswith("PCON")]
        return finish(frame, cases, controls, "gene_id", "ensembl", "counts", p)
    if i == 41:
        p = base / "raw_evidence/GSE191274_Raw_counts.txt.gz"
        with gzip.open(p, "rt", encoding="utf-8") as f:
            headers = next(csv.reader([f.readline()], delimiter="\t"))
        headers = [h.strip('"') for h in headers]
        frame = pd.read_csv(p, sep="\t", skiprows=1, header=None, names=["gene_id"] + headers)
        mp = manifest_groups(manifest, "title", lambda r: "case" if r["group"] == "Case" else "control")
        cases, controls = cols_from_map(list(frame.columns[1:]), mp)
        return finish(frame, cases, controls, "gene_id", "ensembl", "counts", p)
    if i == 42:
        p = base / "GSE174056_Exon_Normalized_Counts.csv.gz"
        frame = pd.read_csv(p, index_col=0).reset_index(names="gene_id")
        cases = [c for c in frame.columns if str(c).startswith("SMA")]
        controls = [c for c in frame.columns if str(c).startswith("Ctrl")]
        return finish(frame, cases, controls, "gene_id", "ensembl", "normalized", p)
    if i == 45:
        p = base / "GSE65581_series_matrix.txt.gz"
        frame, _ = read_geo_series_matrix(p)
        mp = manifest_groups(manifest, "gsm", lambda r: r.get("case_control") if r.get("case_control") in {"case", "control"} else None)
        cases, controls = cols_from_map(list(frame.columns[1:]), mp)
        return finish(frame, cases, controls, "feature_id", "probe", "normalized", p, "GPL16384", feature_space="miRNA_probe")
    if i == 46:
        p = base / "raw_evidence/GSE244401_processed_data.txt.gz"
        frame = pd.read_csv(p, sep="\t")
        mp = manifest_groups(manifest, "sample_title", lambda r: "case" if r["group"] == "SCA" and r["timepoint"] == "T1" else ("control" if r["group"] == "healthy_control" and r["timepoint"] == "T1" else None))
        cases, controls = cols_from_map(list(frame.columns[1:]), mp)
        return finish(frame, cases, controls, "ID", "ensembl", "counts", p)
    if i == 47:
        p = base / "GSE51405_series_matrix.txt.gz"
        frame, _ = read_geo_series_matrix(p)
        mp = manifest_groups(manifest, "gsm", lambda r: "case" if r["illness"].startswith("CVID") else ("control" if r["illness"] == "Healthy" else None))
        cases, controls = cols_from_map(list(frame.columns[1:]), mp)
        return finish(frame, cases, controls, "feature_id", "probe", "normalized", p, "GPL6947")
    if i in {48, 50, 54}:
        cfg = {
            48: (base / "GSE51404_non-normalized_CVID_TST.txt.gz", "GPL6947", "array_id", lambda r: "case" if r["group"] == "XLA" else "control"),
            50: (base / "raw_evidence/GSE44969_non-normalized.txt.gz", "GPL6947", "title", lambda r: "case" if r["selected_primary"] == "TRUE" and "patient" in r["group"].lower() else ("control" if r["selected_primary"] == "TRUE" and "control" in r["group"].lower() else None)),
            54: (base / "GSE56495_non-normalized.txt.gz", "GPL10558", "array_id", lambda r: "case" if r["case_control"] == "case" else "control"),
        }
        p, platform, key, labeller = cfg[i]
        frame = extract_signal_columns(pd.read_csv(p, sep="\t"))
        mp = manifest_groups(manifest, key, labeller)
        normalized = {}
        for c in frame.columns[1:]:
            key_name = str(c).replace(".AVG_Signal", "")
            normalized[str(c)] = mp.get(key_name) or mp.get(str(c))
        cases = [c for c, g in normalized.items() if g == "case"]
        controls = [c for c, g in normalized.items() if g == "control"]
        return finish(frame, cases, controls, frame.columns[0], "probe", "normalized", p, platform)
    if i == 55:
        p = base / "GSE226238_Morrison_et_al_processed_data.xlsx"
        frame0 = pd.read_excel(p, sheet_name="processed data final")
        trans = frame0.set_index(frame0.columns[0]).T.reset_index(names="feature_id")
        mp = manifest_groups(manifest, "library_name", lambda r: "case" if r["selection"] == "primary" and r["group"] == "SCI" else ("control" if r["selection"] == "primary" and r["group"] == "CTL" else None))
        cases = [c for c in trans.columns[1:] if mp.get(str(c)) == "case"]
        controls = [c for c in trans.columns[1:] if mp.get(str(c)) == "control"]
        return finish(trans, cases, controls, "feature_id", "refseq", "normalized", p)
    if i == 58:
        p = ROOT / "phase1b_EU_coordination_v1_2026-08-23/02_worker_packages/worker_gse202518_20260824/GSE202518_RAW.tar"
        columns = {}
        with tarfile.open(p) as tar:
            for member in tar.getmembers():
                raw = tar.extractfile(member).read()
                with gzip.GzipFile(fileobj=io.BytesIO(raw)) as z:
                    d = pd.read_csv(z, sep=None, engine="python")
                columns[member.name.split("_")[0]] = pd.Series(pd.to_numeric(d.iloc[:, 1], errors="coerce").values, index=d.iloc[:, 0].astype(str))
        frame = pd.DataFrame(columns).reset_index(names="feature_id")
        mp = manifest_groups(manifest, "geo_accession", lambda r: "case" if r["primary_role"] == "case" else "control")
        cases, controls = cols_from_map(list(frame.columns[1:]), mp)
        return finish(frame, cases, controls, "feature_id", "symbol", "normalized", p)
    raise NotImplementedError(f"No loader for row {i} {eu}")


def save_eu(result: EUResult) -> Path:
    EU_OUT.mkdir(parents=True, exist_ok=True)
    path = EU_OUT / f"{result.eu_id}_E_U.tsv.gz"
    d = pd.DataFrame({
        "gene_id": result.genes,
        "canonical_symbol": result.symbols,
        "effect": result.effect,
        "uncertainty": result.uncertainty,
        "p_value": result.p_value,
        "adjusted_p_value": result.adj_p_value,
        "n_case": result.n_case,
        "n_control": result.n_control,
        "contrast_id": result.eu_id,
        "analysis_method": result.method,
        "source_path": result.source,
        "feature_space": result.feature_space,
    })
    d.to_csv(path, sep="\t", index=False, compression="gzip")
    return path


def load_saved_eu(eu_id: str) -> EUResult:
    path = EU_OUT / f"{eu_id}_E_U.tsv.gz"
    d = pd.read_csv(path, sep="\t", compression="gzip")
    first = d.iloc[0]
    return EUResult(
        eu_id,
        d["gene_id"].astype(str).to_numpy(object),
        d["canonical_symbol"].fillna("").astype(str).to_numpy(object),
        d["effect"].to_numpy(float),
        d["uncertainty"].to_numpy(float),
        d["p_value"].to_numpy(float),
        d["adjusted_p_value"].to_numpy(float),
        int(first["n_case"]),
        int(first["n_control"]),
        str(first["analysis_method"]),
        str(first["source_path"]),
        str(first["feature_space"]),
    )


def make_folds(manifest: pd.DataFrame, n_folds: int = 5) -> dict[int, list[int]]:
    groups = manifest.groupby("独立性组").agg(category=("大类", "first"), n=("EU_ID", "size")).reset_index()
    groups["hash"] = groups["独立性组"].map(lambda x: hashlib.sha256((str(SEED) + str(x)).encode()).hexdigest())
    folds = {i: [] for i in range(n_folds)}
    load = np.zeros(n_folds)
    for _, g in groups.sort_values(["category", "hash"]).iterrows():
        f = int(np.argmin(load))
        idx = manifest.index[manifest["独立性组"] == g["独立性组"]].tolist()
        folds[f].extend(idx)
        load[f] += len(idx)
    return folds


def weighted_factor_fit(E: np.ndarray, U: np.ndarray, row_w: np.ndarray, k: int, iterations: int = 30) -> tuple[np.ndarray, np.ndarray]:
    mask = np.isfinite(E) & np.isfinite(U) & (U > 0)
    y = np.nan_to_num(E, nan=0.0)
    weights = np.where(mask, row_w[:, None] / np.maximum(U, 1e-12) ** 2, 0.0)
    pos = weights[weights > 0]
    if pos.size:
        weights = np.minimum(weights, np.quantile(pos, 0.99))
    filled = y.copy()
    gene_means = np.divide(np.sum(weights * y, axis=0), np.sum(weights, axis=0), out=np.zeros(y.shape[1]), where=np.sum(weights, axis=0) > 0)
    filled[~mask] = np.broadcast_to(gene_means, filled.shape)[~mask]
    _, _, vt = np.linalg.svd(filled - gene_means, full_matrices=False)
    H = vt[:k].copy()
    W = np.zeros((y.shape[0], k))
    ridge = 1e-6
    for _ in range(iterations):
        for i in range(y.shape[0]):
            z = weights[i]
            W[i] = np.linalg.solve((H * z) @ H.T + ridge * np.eye(k), (H * z) @ y[i])
        for g in range(y.shape[1]):
            z = weights[:, g]
            H[:, g] = np.linalg.solve((W.T * z) @ W + ridge * np.eye(k), (W.T * z) @ y[:, g])
        norms = np.linalg.norm(H, axis=1)
        norms[norms == 0] = 1
        H /= norms[:, None]
        W *= norms[None, :]
    return W, H


def weighted_project(H: np.ndarray, E: np.ndarray, U: np.ndarray, row_w: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    k = H.shape[0]
    mask = np.isfinite(E) & np.isfinite(U) & (U > 0)
    y = np.nan_to_num(E, nan=0.0)
    weights = np.where(mask, row_w[:, None] / np.maximum(U, 1e-12) ** 2, 0.0)
    pos = weights[weights > 0]
    if pos.size:
        weights = np.minimum(weights, np.quantile(pos, 0.99))
    W = np.zeros((len(E), k))
    pred = np.zeros_like(y)
    for i in range(len(E)):
        z = weights[i]
        W[i] = np.linalg.solve((H * z) @ H.T + 1e-6 * np.eye(k), (H * z) @ y[i])
        pred[i] = W[i] @ H
    num = np.sum(weights * (y - pred) ** 2)
    den = np.sum(weights * y**2)
    return W, pred, float(num / den) if den > 0 else math.nan


def run_model(results: dict[str, EUResult], manifest: pd.DataFrame) -> dict[str, object]:
    mrna = manifest[manifest["EU_ID"].map(lambda x: results[str(x)].feature_space == "mRNA_Entrez")].copy().reset_index(drop=True)
    gene_counts: dict[str, int] = {}
    for eu in mrna["EU_ID"]:
        for g in set(results[str(eu)].genes.astype(str)):
            gene_counts[g] = gene_counts.get(g, 0) + 1
    threshold = math.ceil(0.80 * len(mrna))
    genes = sorted(g for g, n in gene_counts.items() if n >= threshold)
    if len(genes) < 300:
        threshold = math.ceil(0.65 * len(mrna))
        genes = sorted(g for g, n in gene_counts.items() if n >= threshold)
    if len(genes) < 100:
        raise ValueError(f"Insufficient harmonized genes: {len(genes)}")
    gpos = {g: j for j, g in enumerate(genes)}
    E = np.full((len(mrna), len(genes)), np.nan)
    U = np.full_like(E, np.nan)
    for i, eu in enumerate(mrna["EU_ID"]):
        r = results[str(eu)]
        for g, e, u in zip(r.genes.astype(str), r.effect, r.uncertainty):
            if g in gpos:
                E[i, gpos[g]] = e
                U[i, gpos[g]] = u
    # Put heterogeneous deposited platforms on a comparable contrast-wise robust scale.
    # The transform is outcome-blind within each already-defined case-control contrast.
    for i in range(len(mrna)):
        valid = np.isfinite(E[i]) & np.isfinite(U[i]) & (U[i] > 0)
        center = float(np.nanmedian(E[i, valid]))
        mad = float(np.nanmedian(np.abs(E[i, valid] - center)))
        scale = max(1.4826 * mad, 1e-6)
        E[i, valid] = np.clip((E[i, valid] - center) / scale, -8.0, 8.0)
        U[i, valid] = U[i, valid] / scale
        u_floor = float(np.nanquantile(U[i, valid], 0.10))
        U[i, valid] = np.maximum(U[i, valid], max(u_floor, 1e-6))
    model_gene_cap = 1000
    if len(genes) > model_gene_cap:
        coverage_score = np.sum(np.isfinite(E) & np.isfinite(U) & (U > 0), axis=0)
        precision_score = np.nanmedian(np.where(np.isfinite(U) & (U > 0), 1.0 / U, np.nan), axis=0)
        precision_score = np.nan_to_num(precision_score, nan=0.0, posinf=0.0, neginf=0.0)
        order = np.lexsort((np.array(genes, dtype=str), -precision_score, -coverage_score))[:model_gene_cap]
        genes = [genes[j] for j in order]
        E = E[:, order]
        U = U[:, order]
    group_sizes = mrna.groupby("独立性组")["EU_ID"].transform("size").to_numpy(float)
    row_w = 1.0 / group_sizes
    folds = make_folds(mrna)
    cv_rows = []
    max_k = min(8, len(mrna) - 2, len(genes))
    for k in range(1, max_k + 1):
        fold_nre = []
        for fold, test_idx in folds.items():
            test = np.array(test_idx, dtype=int)
            train = np.array([i for i in range(len(mrna)) if i not in set(test_idx)], dtype=int)
            _, H = weighted_factor_fit(E[train], U[train], row_w[train], k, iterations=10)
            _, _, nre = weighted_project(H, E[test], U[test], row_w[test])
            fold_nre.append(nre)
            cv_rows.append({"K": k, "fold": fold + 1, "NRE": nre, "test_contrasts": len(test), "test_groups": mrna.iloc[test]["独立性组"].nunique()})
    cv = pd.DataFrame(cv_rows)
    summary = cv.groupby("K")["NRE"].agg(["mean", "std", "count"]).reset_index()
    best = summary.loc[summary["mean"].idxmin()]
    one_se = float(best["std"] / math.sqrt(max(int(best["count"]), 1)))
    eligible = summary[summary["mean"] <= float(best["mean"]) + one_se]
    chosen_k = int(eligible["K"].min())
    W, H = weighted_factor_fit(E, U, row_w, chosen_k, iterations=30)
    _, pred, full_nre = weighted_project(H, E, U, row_w)
    scores = mrna[["序号", "大类", "所选EU病种", "EU_ID", "Accession", "独立性组"]].copy()
    for k in range(chosen_k):
        scores[f"C{k+1}_score"] = W[:, k]
    loadings = pd.DataFrame({"gene_id": genes, "symbol": [ORG["entrez_symbol"].get(g, "") for g in genes]})
    for k in range(chosen_k):
        loadings[f"C{k+1}_loading"] = H[k]
        loadings[f"C{k+1}_abs_rank"] = stats.rankdata(-np.abs(H[k]), method="ordinal").astype(int)
    coverage = pd.DataFrame({"gene_id": genes, "symbol": [ORG["entrez_symbol"].get(g, "") for g in genes], "contrast_count": [gene_counts[g] for g in genes]})
    category_rows = []
    for cat, sub in scores.groupby("大类"):
        row = {"大类": cat, "contrast_count": len(sub), "independence_groups": sub["独立性组"].nunique()}
        for k in range(chosen_k):
            row[f"C{k+1}_mean"] = sub[f"C{k+1}_score"].mean()
            row[f"C{k+1}_sd"] = sub[f"C{k+1}_score"].std(ddof=1)
        category_rows.append(row)
    categories = pd.DataFrame(category_rows)
    # Leave-one-independence-group reconstruction diagnostics.
    logo = []
    for group in mrna["独立性组"].drop_duplicates():
        test = np.where(mrna["独立性组"].to_numpy() == group)[0]
        train = np.where(mrna["独立性组"].to_numpy() != group)[0]
        _, h = weighted_factor_fit(E[train], U[train], row_w[train], chosen_k, iterations=10)
        _, _, nre = weighted_project(h, E[test], U[test], row_w[test])
        logo.append({"independence_group": group, "heldout_contrasts": len(test), "NRE": nre})
    logo = pd.DataFrame(logo)
    OUT.mkdir(parents=True, exist_ok=True)
    cv.to_csv(OUT / "model_cv_by_fold.tsv", sep="\t", index=False)
    summary.to_csv(OUT / "model_rank_selection.tsv", sep="\t", index=False)
    scores.to_csv(OUT / "model_contrast_scores.tsv", sep="\t", index=False)
    loadings.to_csv(OUT / "model_gene_loadings.tsv.gz", sep="\t", index=False, compression="gzip")
    coverage.to_csv(OUT / "model_gene_coverage.tsv", sep="\t", index=False)
    categories.to_csv(OUT / "model_category_scores.tsv", sep="\t", index=False)
    logo.to_csv(OUT / "model_leave_group_out.tsv", sep="\t", index=False)
    np.savetxt(OUT / "model_effect_matrix.tsv.gz", E, delimiter="\t")
    np.savetxt(OUT / "model_uncertainty_matrix.tsv.gz", U, delimiter="\t")
    return {
        "main_model_contrasts": len(mrna),
        "main_model_independence_groups": int(mrna["独立性组"].nunique()),
        "isolated_miRNA_contrasts": int(len(manifest) - len(mrna)),
        "gene_coverage_threshold_contrasts": threshold,
        "model_gene_count": len(genes),
        "chosen_K": chosen_k,
        "best_mean_cv_NRE": float(best["mean"]),
        "chosen_K_mean_cv_NRE": float(summary.loc[summary["K"] == chosen_k, "mean"].iloc[0]),
        "full_fit_NRE": full_nre,
        "median_leave_group_out_NRE": float(logo["NRE"].median()),
        "claim_boundary": "exploratory shared whole-blood bulk-transcriptional structure; not a universal, cell-intrinsic, diagnostic, or confirmatory disease mechanism",
    }


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest().upper()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    table = pd.read_excel(BOOK, sheet_name="58条EU选择", header=2)
    results: dict[str, EUResult] = {}
    audit = []
    for _, row in table.iterrows():
        eu_id = str(row["EU_ID"])
        try:
            cached = EU_OUT / f"{eu_id}_E_U.tsv.gz"
            result = load_saved_eu(eu_id) if cached.exists() else (load_existing_eu(row) if "E/U数值锁" in str(row["锁定层级"]) else load_missing_eu(row))
            out_path = save_eu(result)
            results[eu_id] = result
            audit.append({
                "序号": int(row["序号"]), "EU_ID": eu_id, "status": "PASS", "feature_space": result.feature_space,
                "gene_count": len(result.genes), "n_case": result.n_case, "n_control": result.n_control,
                "source_state": "existing_locked" if "E/U数值锁" in str(row["锁定层级"]) else "new_exploratory_EU",
                "method": result.method, "source_path": result.source, "output_path": str(out_path), "limitation": "",
            })
            print(f"PASS {int(row['序号']):02d} {eu_id} genes={len(result.genes)} n={result.n_case}/{result.n_control}", flush=True)
        except Exception as exc:
            audit.append({"序号": int(row["序号"]), "EU_ID": eu_id, "status": "FAIL", "feature_space": "", "gene_count": 0, "n_case": row["病例N"], "n_control": row["对照N"], "source_state": "", "method": "", "source_path": row["证据路径"], "output_path": "", "limitation": repr(exc)})
            print(f"FAIL {int(row['序号']):02d} {eu_id}: {exc!r}", flush=True)
    audit_df = pd.DataFrame(audit)
    audit_df.to_csv(OUT / "eu58_computation_audit.tsv", sep="\t", index=False)
    if len(results) != len(table):
        raise RuntimeError(f"E/U completed {len(results)}/{len(table)}; inspect audit")
    model_summary = run_model(results, table)
    model_summary.update({
        "schema": "eu58_full_disease_model_v1",
        "date": "2026-08-27",
        "requested_slots": len(table),
        "eu_computed_or_loaded": len(results),
        "new_exploratory_eu": int((audit_df["source_state"] == "new_exploratory_EU").sum()),
        "existing_locked_eu": int((audit_df["source_state"] == "existing_locked").sum()),
        "seed": SEED,
    })
    (OUT / "model_summary.json").write_text(json.dumps(model_summary, ensure_ascii=False, indent=2), encoding="utf-8")
    report = [
        "# 58 槽位统一 E/U 与疾病共有结构模型（探索性）",
        "",
        f"- E/U 完成或载入：{len(results)}/58；其中新补算 {model_summary['new_exploratory_eu']}，复用既有数值锁 {model_summary['existing_locked_eu']}。",
        f"- 共有 mRNA 模型：{model_summary['main_model_contrasts']} 条、{model_summary['main_model_independence_groups']} 个独立性组、{model_summary['model_gene_count']} 个高覆盖 Entrez 基因。",
        f"- 交叉验证选择 K={model_summary['chosen_K']}；该 K 平均 NRE={model_summary['chosen_K_mean_cv_NRE']:.4f}；整组留出中位 NRE={model_summary['median_leave_group_out_NRE']:.4f}。",
        "- GSE65581 为 miRNA 芯片：已计算其 E/U，但无法与 57 条 mRNA E/U 在同一基因特征空间联合分解，因此单独保留，不虚假拼接。",
        "- 新补算 E/U 使用统一探索性流程（按输入尺度转换、病例-对照效应、Welch 不确定性、结局盲最高均值探针折叠）；它们不是对 Phase 1b 正式 EU/Frozen-2 锁的替代。",
        f"- 结论边界：{model_summary['claim_boundary']}。",
    ]
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    files = sorted(p for p in OUT.rglob("*") if p.is_file())
    manifest_rows = [{"path": p.relative_to(OUT).as_posix(), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in files]
    pd.DataFrame(manifest_rows).to_csv(OUT / "artifact_manifest.tsv", sep="\t", index=False)
    print(json.dumps(model_summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
