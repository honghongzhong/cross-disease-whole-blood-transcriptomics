"""Audit and rebuild the 57 mRNA E/U contrasts from sample-level data.

This is a versioned repair workflow.  It never overwrites the existing dataset
E/U files or their mapping manifest.  Repaired effects are Hedges' g and their
sampling standard errors, produced after study-local expression preprocessing.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import runpy
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any
from unittest.mock import patch

sys.path.append(r"C:\Users\peter\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\Lib\site-packages")

import numpy as np
import pandas as pd
from scipy import special, stats


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATASET = ROOT / "dataset"
STATUS = ROOT / "outputs/eu58_full_disease_model_20260827/workbook_58_status.tsv"
ORIGINAL_MAP = DATASET / "EU57_MRNA_DATASET_MAPPING_20260827.tsv"
LEGACY_LOADER = ROOT / "tmp/run_eu58_full_model.py"
F0_LEDGER = ROOT / "census_step0/f0_canonical_cohort_2026-08-22_v1/canonical_contrast_ledger_v1.tsv"
PHASE1B_LEDGER = ROOT / "phase1b_EU_coordination_v1_2026-08-23/00_control/candidate_ledger.tsv"
DEFAULT_EU_OUTPUT = DATASET / "eu57_repaired_20260827_v1"
DEFAULT_REPORT_OUTPUT = HERE / "results"


@dataclass
class RepairedEU:
    eu_id: str
    genes: np.ndarray
    symbols: np.ndarray
    effect: np.ndarray
    uncertainty: np.ndarray
    p_value: np.ndarray
    adjusted_p_value: np.ndarray
    n_case: int
    n_control: int
    method: str
    expression_transform: str
    source_path: str
    selection_source: str
    rebuild_route: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eu-output", type=Path, default=DEFAULT_EU_OUTPUT)
    parser.add_argument("--report-output", type=Path, default=DEFAULT_REPORT_OUTPUT)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def import_legacy() -> ModuleType:
    spec = importlib.util.spec_from_file_location("eu57_legacy_loader", LEGACY_LOADER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import legacy loader: {LEGACY_LOADER}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def bh_adjust(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def quantile_normalize(matrix: np.ndarray) -> np.ndarray:
    """Quantile-normalize columns while preserving missing cells."""
    frame = pd.DataFrame(matrix)
    ranks = frame.rank(method="min", na_option="keep").astype("Int64")
    sorted_values = np.sort(matrix, axis=0)
    rank_means = np.nanmean(sorted_values, axis=1)
    out = np.full(matrix.shape, np.nan, dtype=float)
    for column in range(matrix.shape[1]):
        valid = ranks.iloc[:, column].notna().to_numpy()
        rank_index = ranks.iloc[valid, column].astype(int).to_numpy() - 1
        out[valid, column] = rank_means[rank_index]
    return out


def preprocess_expression(raw: np.ndarray, data_kind: str, source: str) -> tuple[np.ndarray, str]:
    raw = np.asarray(raw, dtype=float)
    if data_kind == "counts":
        raw = np.maximum(raw, 0)
        library = np.nansum(raw, axis=0)
        library[library <= 0] = 1.0
        return np.log2((raw / library[None, :]) * 1e6 + 0.5), "library_size_log2_CPM_plus_0.5"
    finite = raw[np.isfinite(raw)]
    if finite.size == 0:
        return raw, "no_finite_values"
    lower_source = source.lower()
    if "non-normalized" in lower_source or "avg_signal" in lower_source:
        logged = np.log2(np.maximum(raw, 1.0))
        return quantile_normalize(logged), "log2_then_quantile_normalized"
    if np.nanmin(finite) >= 0 and np.nanpercentile(finite, 95) > 100:
        return np.log2(raw + 1.0), "log2_x_plus_1"
    return raw, "as_deposited_normalized_expression"


def hedges_from_matrix(
    eu_id: str,
    genes: np.ndarray,
    symbols: np.ndarray,
    case: np.ndarray,
    control: np.ndarray,
    method: str,
    transform: str,
    source: str,
    selection_source: str,
    route: str,
) -> RepairedEU:
    case = np.asarray(case, dtype=float)
    control = np.asarray(control, dtype=float)
    if case.ndim != 2 or control.ndim != 2 or case.shape[0] != control.shape[0]:
        raise ValueError(f"{eu_id}: case/control matrices are incompatible")
    n_case = np.sum(np.isfinite(case), axis=1)
    n_control = np.sum(np.isfinite(control), axis=1)
    mean_case = np.nanmean(case, axis=1)
    mean_control = np.nanmean(control, axis=1)
    var_case = np.nanvar(case, axis=1, ddof=1)
    var_control = np.nanvar(control, axis=1, ddof=1)
    degrees = n_case + n_control - 2
    pooled_variance = np.divide(
        (n_case - 1) * var_case + (n_control - 1) * var_control,
        degrees,
        out=np.full(case.shape[0], np.nan),
        where=degrees > 0,
    )
    pooled_sd = np.sqrt(pooled_variance)
    mean_difference = mean_case - mean_control
    cohen_d = np.divide(mean_difference, pooled_sd, out=np.full(case.shape[0], np.nan), where=pooled_sd > 0)
    log_j = special.gammaln(degrees / 2) - 0.5 * np.log(degrees / 2) - special.gammaln((degrees - 1) / 2)
    correction = np.exp(log_j)
    effect = correction * cohen_d
    variance = np.square(correction) * (
        np.divide(n_case + n_control, n_case * n_control, out=np.full(case.shape[0], np.nan), where=(n_case * n_control) > 0)
        + np.divide(np.square(cohen_d), 2 * degrees, out=np.full(case.shape[0], np.nan), where=degrees > 0)
    )
    uncertainty = np.sqrt(variance)
    welch = stats.ttest_ind(case, control, axis=1, equal_var=False, nan_policy="omit")
    p_value = np.asarray(welch.pvalue, dtype=float)
    keep = (
        (n_case >= 2)
        & (n_control >= 2)
        & np.isfinite(effect)
        & np.isfinite(uncertainty)
        & (uncertainty > 0)
        & np.isfinite(p_value)
    )
    genes = np.asarray(genes, dtype=object)[keep].astype(str)
    symbols = np.asarray(symbols, dtype=object)[keep].astype(str)
    effect = effect[keep]
    uncertainty = uncertainty[keep]
    p_value = p_value[keep]
    if len(genes) == 0:
        raise ValueError(f"{eu_id}: no valid Hedges g rows")
    # One row per Entrez gene; highest all-sample mean was already used in raw-loader routes.
    order = np.argsort(uncertainty, kind="stable")
    unique = ~pd.Series(genes[order]).duplicated().to_numpy()
    idx = order[unique]
    genes, symbols = genes[idx], symbols[idx]
    effect, uncertainty, p_value = effect[idx], uncertainty[idx], p_value[idx]
    adjusted = bh_adjust(p_value)
    return RepairedEU(
        eu_id=eu_id,
        genes=genes,
        symbols=symbols,
        effect=effect,
        uncertainty=uncertainty,
        p_value=p_value,
        adjusted_p_value=adjusted,
        n_case=int(case.shape[1]),
        n_control=int(control.shape[1]),
        method=method,
        expression_transform=transform,
        source_path=source,
        selection_source=selection_source,
        rebuild_route=route,
    )


def auto_map_features(values: Any, legacy: ModuleType) -> tuple[np.ndarray, np.ndarray]:
    mapped: list[str] = []
    symbols: list[str] = []
    for value in map(str, values):
        raw = value.strip().replace(".0", "")
        gene = ""
        if raw in legacy.ORG["valid_entrez"]:
            gene = raw
        elif raw.upper().startswith("ENSG"):
            gene = legacy.ORG["ensembl"].get(raw.split(".")[0].upper(), "")
        else:
            gene = legacy.ORG["symbol"].get(raw.upper(), "")
        mapped.append(gene)
        symbols.append(legacy.ORG["entrez_symbol"].get(gene, "") if gene else "")
    return np.asarray(mapped, dtype=object), np.asarray(symbols, dtype=object)


NATIVE_PROBE_PLATFORMS = {
    "gse17048_1": "GPL6947",
    "gse48060_1": "GPL570",
    "gse137340_1": "GPL10558",
}


def map_native_features(values: Any, dataset_dir: str, legacy: ModuleType) -> tuple[np.ndarray, np.ndarray]:
    platform = NATIVE_PROBE_PLATFORMS.get(dataset_dir)
    if not platform:
        return auto_map_features(values, legacy)
    probe_map = legacy.PROBE_MAPS[platform]
    genes = np.asarray([probe_map.get(str(value), "") for value in values], dtype=object)
    symbols = np.asarray([legacy.ORG["entrez_symbol"].get(str(gene), "") if gene else "" for gene in genes], dtype=object)
    return genes, symbols


def standardized_legacy_builder(legacy: ModuleType):
    def build(
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
    ) -> RepairedEU:
        if feature_space != "mRNA_Entrez":
            raise ValueError(f"{eu_id}: non-mRNA feature space is outside the EU57 repair")
        if len(case_cols) < 2 or len(control_cols) < 2:
            raise ValueError(f"{eu_id}: insufficient sample groups")
        columns = case_cols + control_cols
        missing = [column for column in columns if column not in frame.columns]
        if missing:
            raise ValueError(f"{eu_id}: missing sample columns: {missing[:8]}")
        raw = frame[columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        features = frame[feature_col].astype(str).to_numpy(object)
        mapped = np.asarray([legacy.map_feature(value, feature_type, platform) for value in features], dtype=object)
        if symbol_col and symbol_col in frame.columns:
            symbols = frame[symbol_col].fillna("").astype(str).to_numpy(object)
        else:
            symbols = np.asarray([legacy.ORG["entrez_symbol"].get(str(gene), "") if gene else "" for gene in mapped], dtype=object)
        keep = np.asarray([gene is not None and str(gene) != "" for gene in mapped]) & np.isfinite(raw).any(axis=1)
        raw, mapped, symbols = raw[keep], mapped[keep].astype(str), symbols[keep]
        all_sample_mean = np.nanmean(raw, axis=1)
        chosen: dict[str, int] = {}
        for index, gene in enumerate(mapped):
            if gene not in chosen or all_sample_mean[index] > all_sample_mean[chosen[gene]]:
                chosen[gene] = index
        idx = np.asarray(sorted(chosen.values()), dtype=int)
        raw, mapped, symbols = raw[idx], mapped[idx], symbols[idx]
        expression, transform = preprocess_expression(raw, data_kind, source)
        return hedges_from_matrix(
            eu_id,
            mapped,
            symbols,
            expression[:, : len(case_cols)],
            expression[:, len(case_cols) :],
            "Hedges_g_from_sample_matrix; Welch_p; outcome_blind_feature_collapse",
            transform,
            source,
            "phase1b signed sample/selection manifest used by legacy loader",
            "direct_sample_matrix_phase1b_loader",
        )

    return build


def capture_native_main(script: Path) -> dict[str, Any]:
    module = runpy.run_path(str(script), run_name=f"capture_{script.parent.parent.name}")
    main = module.get("main")
    if main is None:
        raise ValueError(f"No main() in native script: {script}")
    captured: dict[str, Any] = {}

    def tracer(frame, event, arg):
        if event == "return" and frame.f_code is main.__code__:
            captured.update(frame.f_locals)
        return tracer

    original_trace = sys.gettrace()
    with patch.object(pd.DataFrame, "to_csv", autospec=True, return_value=None):
        sys.settrace(tracer)
        try:
            main()
        finally:
            sys.settrace(original_trace)
    if not captured:
        raise RuntimeError(f"Failed to capture native sample matrices: {script}")
    return captured


def repair_native_welch(row: pd.Series, legacy: ModuleType) -> RepairedEU:
    dataset_dir = DATASET / str(row["DatasetDirectory"])
    script = dataset_dir / "process/run_analysis.py"
    local = capture_native_main(script)
    case = local.get("A", local.get("a"))
    control = local.get("B", local.get("b"))
    frame = local.get("x")
    if case is None or control is None or not isinstance(frame, pd.DataFrame):
        raise ValueError(f"{row['EU_ID']}: native script did not expose case/control/x matrices")
    genes, symbols = map_native_features(frame.index, str(row["DatasetDirectory"]), legacy)
    keep = genes != ""
    source_files = sorted((dataset_dir / "raw").glob("*"))
    selection_files = sorted((dataset_dir / "process").glob("*manifest*.tsv"))
    return hedges_from_matrix(
        str(row["EU_ID"]),
        genes[keep],
        symbols[keep],
        np.asarray(case)[keep],
        np.asarray(control)[keep],
        "Hedges_g_from_native_sample_matrix; Welch_p",
        "native_dataset_script_preprocessing",
        ";".join(str(path) for path in source_files),
        ";".join(str(path) for path in selection_files),
        "direct_sample_matrix_dataset_native_script",
    )


def repair_gse186507(row: pd.Series, legacy: ModuleType) -> RepairedEU:
    """Rebuild GSE186507 while preserving the Ensembl identifier column.

    The pre-existing native script coerced the identifier column to numeric and
    consequently changed every gene identifier to zero.  This route fixes that
    defect before any effect calculation.
    """
    dataset_dir = DATASET / "gse186507_1"
    manifest_path = dataset_dir / "process/sample_manifest.tsv"
    raw_path = dataset_dir / "raw/GSE186507_MSCCR_Blood_counts.txt.gz"
    meta = pd.read_csv(manifest_path, sep="\t")
    meta = meta[meta.selected.astype(str).str.upper().eq("TRUE")]
    with pd.io.common.get_handle(raw_path, "rt", compression="gzip", encoding="latin1") as handle:
        names = [value.strip('"') for value in handle.handle.readline().split()]
    case_ids = meta.loc[meta.group.eq("CD"), "subject"].astype(str).tolist()
    control_ids = meta.loc[meta.group.eq("Healthy_Control"), "subject"].astype(str).tolist()
    selected_ids = case_ids + control_ids
    columns = ["gene_id"] + names
    table = pd.read_csv(
        raw_path,
        sep=r"\s+",
        header=None,
        skiprows=1,
        names=columns,
        usecols=["gene_id"] + selected_ids,
        compression="gzip",
    )
    matrix = table[selected_ids].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy(float)
    expression, transform = preprocess_expression(matrix, "counts", str(raw_path))
    genes = np.asarray([legacy.map_feature(value, "ensembl") or "" for value in table.gene_id], dtype=object)
    symbols = np.asarray([legacy.ORG["entrez_symbol"].get(str(gene), "") if gene else "" for gene in genes], dtype=object)
    keep = genes != ""
    genes, symbols, expression = genes[keep], symbols[keep], expression[keep]
    all_mean = np.nanmean(expression, axis=1)
    chosen: dict[str, int] = {}
    for index, gene in enumerate(genes.astype(str)):
        if gene not in chosen or all_mean[index] > all_mean[chosen[gene]]:
            chosen[gene] = index
    index = np.asarray(sorted(chosen.values()), dtype=int)
    genes, symbols, expression = genes[index], symbols[index], expression[index]
    return hedges_from_matrix(
        str(row["EU_ID"]),
        genes,
        symbols,
        expression[:, : len(case_ids)],
        expression[:, len(case_ids) :],
        "Hedges_g_from_sample_matrix; library_size_log2_CPM; Welch_p",
        transform,
        str(raw_path),
        str(manifest_path),
        "direct_sample_matrix_gse186507_identifier_repair",
    )


GSE112057_TARGETS = {
    "gse112057_1": "Crohn's Disease",
    "gse112057_2": "Ulcerative Colitis",
    "gse112057_3": "Systemic JIA",
}


def repair_gse112057(row: pd.Series, legacy: ModuleType) -> RepairedEU:
    dataset_dir = DATASET / str(row["DatasetDirectory"])
    target = GSE112057_TARGETS[str(row["DatasetDirectory"])]
    meta = pd.read_csv(dataset_dir / "process/sample_metadata.tsv", sep="\t", dtype=str)
    raw = pd.read_csv(dataset_dir / "raw/GSE112057_RawCounts_dataset.txt.gz", sep="\t")
    symbols_raw = raw.iloc[:, 0].astype(str).str.replace(r"_\d+$", "", regex=True)
    counts = raw.iloc[:, 1:].apply(pd.to_numeric, errors="raise").astype("int64")
    counts.columns = counts.columns.astype(str)
    counts.index = symbols_raw
    counts = counts.groupby(level=0, sort=False).sum()
    selected = meta[meta["diagnosis"].isin(["Control", target])].copy()
    selected["group"] = selected["diagnosis"].map({"Control": "control", target: "case"})
    selected["raw_sample"] = selected["title"].str.split("_", n=1).str[0]
    counts = counts.loc[:, selected["raw_sample"]]
    counts = counts.loc[counts.sum(axis=1) >= 10]
    library = counts.sum(axis=0).to_numpy(float)
    expression = np.log2((counts.to_numpy(float) / library[None, :]) * 1e6 + 0.5)
    genes, canonical_symbols = auto_map_features(counts.index, legacy)
    keep = genes != ""
    case_columns = np.flatnonzero(selected["group"].eq("case").to_numpy())
    control_columns = np.flatnonzero(selected["group"].eq("control").to_numpy())
    return hedges_from_matrix(
        str(row["EU_ID"]),
        genes[keep],
        canonical_symbols[keep],
        expression[keep][:, case_columns],
        expression[keep][:, control_columns],
        "Hedges_g_from_sample_matrix; library_size_log2_CPM; Welch_p",
        "library_size_log2_CPM_plus_0.5",
        str(dataset_dir / "raw/GSE112057_RawCounts_dataset.txt.gz"),
        str(dataset_dir / "process/sample_metadata.tsv"),
        "direct_sample_matrix_gse112057",
    )


def current_audit(path: Path) -> dict[str, Any]:
    frame = pd.read_csv(path, sep="\t", compression="infer")
    effect = pd.to_numeric(frame.get("effect"), errors="coerce")
    uncertainty = pd.to_numeric(frame.get("uncertainty"), errors="coerce")
    valid = np.isfinite(effect) & np.isfinite(uncertainty) & (uncertainty > 0)
    p = pd.to_numeric(frame.get("p_value"), errors="coerce") if "p_value" in frame else pd.Series(dtype=float)
    return {
        "current_rows": int(len(frame)),
        "current_valid_rows": int(valid.sum()),
        "current_u_unique": int(uncertainty[valid].nunique()),
        "current_p_all_one": bool(len(p.dropna()) > 0 and p.dropna().eq(1).all()),
        "current_effect_mad": float(np.median(np.abs(effect[valid] - np.median(effect[valid])))),
        "current_abs_effect_p95": float(np.quantile(np.abs(effect[valid]), 0.95)),
    }


def save_repaired(result: RepairedEU, output_root: Path, dataset_dir: str) -> Path:
    directory = output_root / dataset_dir
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{result.eu_id}_standardized_E_U_v1.tsv.gz"
    frame = pd.DataFrame(
        {
            "gene_id": result.genes,
            "canonical_symbol": result.symbols,
            "effect": result.effect,
            "uncertainty": result.uncertainty,
            "p_value": result.p_value,
            "adjusted_p_value": result.adjusted_p_value,
            "n_case": result.n_case,
            "n_control": result.n_control,
            "contrast_id": result.eu_id,
            "effect_definition": "Hedges_g_case_minus_control",
            "uncertainty_definition": "sampling_SE_of_Hedges_g",
            "analysis_method": result.method,
            "expression_transform": result.expression_transform,
            "source_path": result.source_path,
            "selection_source": result.selection_source,
            "rebuild_route": result.rebuild_route,
            "feature_space": "mRNA_Entrez",
        }
    ).sort_values("gene_id", kind="stable")
    frame.to_csv(path, sep="\t", index=False, compression="gzip")
    return path


def repaired_audit(result: RepairedEU, expected_case: int, expected_control: int) -> dict[str, Any]:
    finite = np.isfinite(result.effect) & np.isfinite(result.uncertainty) & (result.uncertainty > 0)
    p_ok = np.isfinite(result.p_value) & (result.p_value >= 0) & (result.p_value <= 1)
    q_ok = np.isfinite(result.adjusted_p_value) & (result.adjusted_p_value >= 0) & (result.adjusted_p_value <= 1)
    duplicate_genes = int(pd.Series(result.genes).duplicated().sum())
    checks = {
        "finite_positive_eu": bool(finite.all()),
        "p_in_range": bool(p_ok.all()),
        "q_in_range": bool(q_ok.all()),
        "unique_entrez": duplicate_genes == 0,
        "variable_u": int(pd.Series(result.uncertainty).nunique()) > 1,
        "sample_counts_match": result.n_case == expected_case and result.n_control == expected_control,
        "minimum_gene_rows": len(result.genes) >= 1000,
    }
    return {
        "repair_status": "PASS" if all(checks.values()) else "HOLD",
        "repair_checks": json.dumps(checks, sort_keys=True),
        "repaired_rows": int(len(result.genes)),
        "repaired_u_unique": int(pd.Series(result.uncertainty).nunique()),
        "repaired_effect_mad": float(np.median(np.abs(result.effect - np.median(result.effect)))),
        "repaired_abs_effect_p95": float(np.quantile(np.abs(result.effect), 0.95)),
        "repaired_abs_z_p95": float(np.quantile(np.abs(result.effect / result.uncertainty), 0.95)),
        "rebuild_route": result.rebuild_route,
        "expression_transform": result.expression_transform,
    }


def governance_fields(row: pd.Series, f0: pd.DataFrame, latest_phase1b: pd.DataFrame) -> dict[str, str]:
    eu_id = str(row["EU_ID"])
    accession = str(row["Accession"])
    source_state = str(row["source_state"])
    formal = f0[f0["contrast_id"].eq(eu_id)]
    phase = latest_phase1b[latest_phase1b["accession"].eq(accession)]
    formal_pass = bool(
        len(formal)
        and formal.iloc[-1]["lock_status"] == "PASS"
        and formal.iloc[-1]["inclusion_status"] == "primary_candidate"
    )
    phase_status = str(phase.iloc[-1]["current_status"]) if len(phase) else "not_in_phase1b"
    if source_state == "existing_locked" and formal_pass:
        eligible = "PASS"
        authority = "F0_formal_primary_lock"
        authority_status = "PASS_primary_candidate"
    elif source_state != "existing_locked" and phase_status == "EU_locked":
        eligible = "PASS"
        authority = "Phase1b_latest_accession_status"
        authority_status = "EU_locked"
    else:
        eligible = "HOLD"
        authority = "unresolved"
        authority_status = "no_applicable_current_pass"
    conflict = ""
    if source_state == "existing_locked" and formal_pass and phase_status not in {"EU_locked", "not_in_phase1b"}:
        conflict = (
            f"Phase1b latest status is {phase_status}; it belongs to a separate later candidate process and "
            "does not silently override the pre-existing formal F0 primary lock. Retained for audit."
        )
    return {
        "model_input_eligibility": eligible,
        "eligibility_authority": authority,
        "eligibility_authority_status": authority_status,
        "phase1b_latest_status": phase_status,
        "governance_conflict_note": conflict,
    }


def main() -> None:
    args = parse_args()
    eu_output = args.eu_output.resolve()
    report_output = args.report_output.resolve()
    for target in (eu_output, report_output):
        if target.exists():
            if not args.overwrite:
                raise FileExistsError(f"Refusing to overwrite existing output: {target}")
            shutil.rmtree(target)
    eu_output.mkdir(parents=True)
    report_output.mkdir(parents=True)

    legacy = import_legacy()
    status = pd.read_csv(STATUS, sep="\t")
    status = status[status["纳入说明"].eq("纳入57条mRNA共有结构分解")].copy()
    mapping = pd.read_csv(ORIGINAL_MAP, sep="\t", dtype=str)
    merged = status.merge(mapping, on=["EU_ID", "Accession"], how="inner", validate="one_to_one")
    if len(merged) != 57:
        raise ValueError(f"Expected 57 merged EU rows, got {len(merged)}")
    f0 = pd.read_csv(F0_LEDGER, sep="\t", dtype=str, keep_default_na=False)
    phase1b = pd.read_csv(
        PHASE1B_LEDGER,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        engine="python",
        on_bad_lines="skip",
    )
    latest_phase1b = phase1b.sort_values("last_updated").drop_duplicates("accession", keep="last")

    original_legacy_builder = legacy.build_eu
    legacy.build_eu = standardized_legacy_builder(legacy)
    audit_rows: list[dict[str, Any]] = []
    manifest_rows: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    try:
        for _, row in merged.sort_values("序号").iterrows():
            eu_id = str(row["EU_ID"])
            current_path = Path(str(row["DatasetPath"]))
            base_audit = current_audit(current_path)
            try:
                dataset_dir = str(row["DatasetDirectory"])
                script = DATASET / dataset_dir / "process/run_analysis.py"
                script_text = script.read_text(encoding="utf-8")
                if dataset_dir in GSE112057_TARGETS:
                    result = repair_gse112057(row, legacy)
                elif dataset_dir == "gse186507_1":
                    result = repair_gse186507(row, legacy)
                elif "No source result table found" not in script_text:
                    result = repair_native_welch(row, legacy)
                else:
                    result = legacy.load_missing_eu(row)
                    if not isinstance(result, RepairedEU):
                        raise TypeError(f"Legacy loader returned unexpected type for {eu_id}")
                expected_case = int(row["病例N"])
                expected_control = int(row["对照N"])
                repaired_metrics = repaired_audit(result, expected_case, expected_control)
                governance = governance_fields(row, f0, latest_phase1b)
                repaired_path = save_repaired(result, eu_output, dataset_dir)
                digest = sha256(repaired_path)
                audit_rows.append(
                    {
                        "序号": int(row["序号"]),
                        "EU_ID": eu_id,
                        "Accession": row["Accession"],
                        "DatasetDirectory": dataset_dir,
                        "n_case_expected": expected_case,
                        "n_control_expected": expected_control,
                        **base_audit,
                        **repaired_metrics,
                        **governance,
                        "repaired_path": str(repaired_path),
                        "sha256": digest,
                        "limitation": "training candidate only; no model validation performed",
                    }
                )
                manifest_rows.append(
                    {
                        "EU_ID": eu_id,
                        "Accession": row["Accession"],
                        "DatasetDirectory": dataset_dir,
                        "DatasetPath": str(repaired_path),
                        "SHA256": digest,
                        "repair_status": repaired_metrics["repair_status"],
                        "effect_definition": "Hedges_g_case_minus_control",
                        "uncertainty_definition": "sampling_SE_of_Hedges_g",
                        "rebuild_route": result.rebuild_route,
                        **governance,
                    }
                )
                print(f"[{len(audit_rows):02d}/57] {eu_id}: {repaired_metrics['repair_status']} ({len(result.genes)} genes)", flush=True)
            except Exception as exc:
                failures.append({"EU_ID": eu_id, "DatasetDirectory": str(row["DatasetDirectory"]), "error": repr(exc)})
                print(f"[FAIL] {eu_id}: {exc!r}", flush=True)
    finally:
        legacy.build_eu = original_legacy_builder

    audit = pd.DataFrame(audit_rows).sort_values("序号") if audit_rows else pd.DataFrame()
    manifest = pd.DataFrame(manifest_rows)
    failure_frame = pd.DataFrame(failures, columns=["EU_ID", "DatasetDirectory", "error"])
    audit.to_csv(report_output / "eu57_audit_repair_ledger_v1.tsv", sep="\t", index=False)
    manifest.to_csv(eu_output / "EU57_STANDARDIZED_REPAIRED_MANIFEST_v1.tsv", sep="\t", index=False)
    failure_frame.to_csv(report_output / "repair_failures_v1.tsv", sep="\t", index=False)

    pass_count = int((audit.get("repair_status") == "PASS").sum()) if len(audit) else 0
    hold_count = int((audit.get("repair_status") == "HOLD").sum()) if len(audit) else 0
    eligible_count = int((audit.get("model_input_eligibility") == "PASS").sum()) if len(audit) else 0
    current_mad_ratio = float(audit["current_effect_mad"].max() / audit["current_effect_mad"].min()) if len(audit) else math.nan
    repaired_mad_ratio = float(audit["repaired_effect_mad"].max() / audit["repaired_effect_mad"].min()) if len(audit) else math.nan
    current_p95_ratio = float(audit["current_abs_effect_p95"].max() / audit["current_abs_effect_p95"].min()) if len(audit) else math.nan
    repaired_p95_ratio = float(audit["repaired_abs_effect_p95"].max() / audit["repaired_abs_effect_p95"].min()) if len(audit) else math.nan
    summary = {
        "schema": "eu57_audit_repair_v1",
        "requested_eu": 57,
        "rebuilt_eu": int(len(audit)),
        "pass": pass_count,
        "hold": hold_count,
        "failed": int(len(failures)),
        "model_input_eligible": eligible_count,
        "effect_definition": "Hedges_g_case_minus_control",
        "uncertainty_definition": "sampling_SE_of_Hedges_g",
        "original_assets_overwritten": False,
        "model_validation_performed": False,
        "scale_audit": {
            "effect_mad_max_min_ratio_before": current_mad_ratio,
            "effect_mad_max_min_ratio_after": repaired_mad_ratio,
            "abs_effect_p95_max_min_ratio_before": current_p95_ratio,
            "abs_effect_p95_max_min_ratio_after": repaired_p95_ratio,
        },
    }
    (report_output / "summary_v1.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    report_lines = [
        "# EU57 audit and repair v1",
        "",
        f"- Rebuilt: {len(audit)}/57; PASS: {pass_count}; HOLD: {hold_count}; failed: {len(failures)}.",
        f"- Governance-reconciled model-input eligibility: {eligible_count}/57 PASS.",
        "- Existing dataset E/U and the original mapping manifest were not overwritten.",
        "- Repaired E is Hedges' g (case minus control); repaired U is its sampling standard error.",
        "- Each contrast was rebuilt from a sample-level matrix and signed sample selection evidence.",
        f"- Effect MAD scale ratio improved from {current_mad_ratio:.2f} to {repaired_mad_ratio:.2f}; absolute-effect p95 ratio improved from {current_p95_ratio:.2f} to {repaired_p95_ratio:.2f}.",
        "- This package is a training-input candidate only; validation was not performed.",
    ]
    (report_output / "REPORT.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    artifact_rows = []
    for path in sorted(eu_output.rglob("*")) + sorted(report_output.rglob("*")) + [Path(__file__)]:
        if path.is_file() and path.name != "artifact_manifest_v1.tsv":
            artifact_rows.append(
                {
                    "path": str(path),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    pd.DataFrame(artifact_rows).to_csv(report_output / "artifact_manifest_v1.tsv", sep="\t", index=False)
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
    if failures or hold_count:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
