"""Reproducible nested cross-disease validation for the 57 standardized E/U contrasts.

Primary analysis: purged leave-one-disease-class-out validation.  A held-out
class is removed together with every training contrast sharing its accession or
independence group.  Feature selection and rank selection are repeated using
training data only.  Test scores are estimated from anchor genes and evaluated
on deterministically masked genes.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import platform
import sys
import time
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import norm, rankdata
from sklearn.model_selection import GroupKFold


HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parents[1]
DEFAULT_MANIFEST = PACKAGE_ROOT / "data/input_manifest.tsv"
DEFAULT_AUDIT = PACKAGE_ROOT / "metadata/input_audit.tsv"
DEFAULT_CLASS_MAP = PACKAGE_ROOT / "metadata/disease_class_map.tsv"
DEFAULT_OUTPUT = PACKAGE_ROOT / "outputs/formal_results"
SEED = 20260828


@dataclass
class Model:
    mean: np.ndarray
    scores: np.ndarray
    loadings: np.ndarray
    fitted: np.ndarray
    history: list[float]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    p.add_argument("--class-map", type=Path, default=DEFAULT_CLASS_MAP)
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--genes", type=int, default=1000)
    p.add_argument("--coverage", type=float, default=0.90)
    p.add_argument("--max-k", type=int, default=6)
    p.add_argument("--inner-folds", type=int, default=5)
    p.add_argument("--inner-masks", type=int, default=2)
    p.add_argument("--outer-masks", type=int, default=5)
    p.add_argument("--mask-fraction", type=float, default=0.20)
    p.add_argument("--max-iterations", type=int, default=80)
    p.add_argument("--ridge", type=float, default=1e-3)
    p.add_argument("--bootstrap-repeats", type=int, default=1000)
    p.add_argument("--stability-repeats", type=int, default=200)
    p.add_argument("--quick", action="store_true", help="Small smoke test; not a reportable analysis.")
    return p.parse_args()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest().upper()


def package_relative(path: Path, base: Path = PACKAGE_ROOT) -> str:
    """Return a portable path for metadata written into public result files.

    Standard package inputs are stored relative to the package root. For an
    explicitly supplied external path, only the final name is retained so a
    public result never records a user's local directory structure.
    """
    resolved = path.resolve()
    try:
        return resolved.relative_to(base.resolve()).as_posix()
    except ValueError:
        return resolved.name


def slug(text: str) -> str:
    return "".join(c.lower() if c.isalnum() else "_" for c in text).strip("_")


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    good = np.isfinite(p)
    x = p[good]
    if not len(x):
        return out
    order = np.argsort(x)
    ranked = x[order]
    q = np.minimum.accumulate((ranked * len(x) / np.arange(1, len(x) + 1))[::-1])[::-1]
    q = np.clip(q, 0, 1)
    restored = np.empty_like(q)
    restored[order] = q
    out[good] = restored
    return out


def load_inputs(manifest_path: Path, audit_path: Path, class_map_path: Path):
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str)
    audit = pd.read_csv(audit_path, sep="\t", dtype=str)
    classes = pd.read_csv(class_map_path, sep="\t", dtype=str)
    if len(manifest) != 57 or manifest.EU_ID.nunique() != 57:
        raise ValueError("The manifest must contain exactly 57 unique EU_ID values")
    if set(manifest.EU_ID) != set(classes.EU_ID):
        raise ValueError("Class map EU_ID values do not exactly match the manifest")
    if classes.disease_class.nunique() != 12:
        raise ValueError("Class map must contain exactly 12 disease classes")
    manifest = manifest.merge(classes, on="EU_ID", validate="one_to_one")
    manifest = manifest.merge(
        audit[["EU_ID", "n_case_expected", "n_control_expected", "expression_transform"]],
        on="EU_ID", validate="one_to_one"
    )
    records = []
    all_genes: set[str] = set()
    symbols: dict[str, str] = {}
    for row in manifest.itertuples(index=False):
        path = Path(row.DatasetPath)
        if not path.is_absolute():
            path = manifest_path.resolve().parent / path
        if not path.is_file():
            raise FileNotFoundError(path)
        if sha256(path) != row.SHA256.upper():
            raise ValueError(f"SHA-256 mismatch: {path}")
        frame = pd.read_csv(path, sep="\t", compression="gzip", dtype={"gene_id": str})
        required = {"gene_id", "canonical_symbol", "effect", "uncertainty"}
        if not required.issubset(frame.columns):
            raise ValueError(f"Missing required columns in {path}")
        frame = frame.drop_duplicates("gene_id", keep=False).copy()
        frame.effect = pd.to_numeric(frame.effect, errors="coerce")
        frame.uncertainty = pd.to_numeric(frame.uncertainty, errors="coerce")
        frame = frame[np.isfinite(frame.effect) & np.isfinite(frame.uncertainty) & (frame.uncertainty > 0)]
        s_e = pd.Series(frame.effect.to_numpy(float), index=frame.gene_id.astype(str))
        s_u = pd.Series(frame.uncertainty.to_numpy(float), index=frame.gene_id.astype(str))
        records.append((s_e, s_u))
        all_genes.update(s_e.index)
        for g, s in zip(frame.gene_id.astype(str), frame.canonical_symbol.fillna("").astype(str)):
            if s and g not in symbols:
                symbols[g] = s
    genes = sorted(all_genes, key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x))
    pos = {g: j for j, g in enumerate(genes)}
    effect = np.full((len(records), len(genes)), np.nan)
    uncertainty = np.full_like(effect, np.nan)
    for i, (se, su) in enumerate(records):
        idx = np.array([pos[g] for g in se.index])
        effect[i, idx] = se.to_numpy(float)
        uncertainty[i, idx] = su.reindex(se.index).to_numpy(float)
    gene_table = pd.DataFrame({"gene_id": genes, "canonical_symbol": [symbols.get(g, "") for g in genes]})
    return manifest, gene_table, effect, uncertainty


def select_genes(effect: np.ndarray, uncertainty: np.ndarray, rows: np.ndarray,
                 coverage: float, n_genes: int) -> tuple[np.ndarray, pd.DataFrame]:
    e = effect[rows]
    u = uncertainty[rows]
    observed = np.isfinite(e) & np.isfinite(u) & (u > 0)
    counts = observed.sum(axis=0)
    minimum = math.ceil(coverage * len(rows))
    eligible = counts >= minimum
    active = counts > 0
    med = np.full(effect.shape[1], np.nan)
    mad = np.full(effect.shape[1], np.nan)
    med_u = np.full(effect.shape[1], np.nan)
    med[active] = np.nanmedian(np.where(observed[:, active], e[:, active], np.nan), axis=0)
    mad[active] = np.nanmedian(
        np.where(observed[:, active], np.abs(e[:, active] - med[active]), np.nan), axis=0
    ) * 1.4826
    med_u[active] = np.nanmedian(np.where(observed[:, active], u[:, active], np.nan), axis=0)
    score = np.divide(mad, med_u, out=np.zeros_like(mad), where=np.isfinite(med_u) & (med_u > 0))
    score *= np.sqrt(counts / len(rows))
    score[~eligible | ~np.isfinite(score)] = -np.inf
    candidates = np.flatnonzero(eligible)
    if len(candidates) < n_genes:
        raise ValueError(f"Only {len(candidates)} genes meet coverage in a training fold")
    order = np.lexsort((candidates, -counts[candidates], -score[candidates]))
    chosen = candidates[order[:n_genes]]
    audit = pd.DataFrame({
        "gene_index": np.arange(effect.shape[1]), "training_coverage_count": counts,
        "training_coverage_fraction": counts / len(rows), "training_selection_score": score,
        "eligible": eligible, "selected": False,
    })
    audit.loc[chosen, "selected"] = True
    audit.loc[chosen, "selection_rank"] = np.arange(1, len(chosen) + 1)
    return chosen, audit


def objective_weights(effect: np.ndarray, uncertainty: np.ndarray, groups: np.ndarray,
                      use_uncertainty: bool = True) -> np.ndarray:
    observed = np.isfinite(effect) & np.isfinite(uncertainty) & (uncertainty > 0)
    w = np.zeros_like(effect, dtype=float)
    for i in range(effect.shape[0]):
        m = observed[i]
        if not m.any():
            continue
        raw = 1.0 / np.square(uncertainty[i, m]) if use_uncertainty else np.ones(m.sum())
        lo, hi = np.quantile(raw, [0.05, 0.95])
        clipped = np.clip(raw, lo, hi)
        w[i, m] = clipped / clipped.mean()
    for group in np.unique(groups):
        r = groups == group
        total = w[r].sum()
        if total > 0:
            w[r] /= total
    if w.sum() > 0:
        w *= observed.sum() / w.sum()
    return w


def solve_wls(design: np.ndarray, target: np.ndarray, weights: np.ndarray,
              ridge: float, penalize_intercept: bool = True) -> np.ndarray:
    penalty = np.eye(design.shape[1]) * ridge
    if not penalize_intercept:
        penalty[0, 0] = 0
    gram = design.T @ (weights[:, None] * design) + penalty
    rhs = design.T @ (weights * target)
    return np.linalg.solve(gram, rhs)


def orient(scores: np.ndarray, loadings: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    u, s, vt = np.linalg.svd(scores @ loadings, full_matrices=False)
    rank = scores.shape[1]
    scores = u[:, :rank] * s[:rank]
    loadings = vt[:rank]
    for k in range(rank):
        j = np.argmax(np.abs(loadings[k]))
        if loadings[k, j] < 0:
            scores[:, k] *= -1
            loadings[k] *= -1
    return scores, loadings


def fit_model(effect: np.ndarray, uncertainty: np.ndarray, groups: np.ndarray, rank: int,
              use_uncertainty: bool, max_iterations: int, ridge: float) -> Model:
    weights = objective_weights(effect, uncertainty, groups, use_uncertainty)
    observed = weights > 0
    mean = np.divide((weights * np.nan_to_num(effect)).sum(0), weights.sum(0),
                     out=np.zeros(effect.shape[1]), where=weights.sum(0) > 0)
    filled = np.where(observed, effect - mean, 0.0)
    u, s, vt = np.linalg.svd(filled, full_matrices=False)
    scores, loadings = u[:, :rank] * s[:rank], vt[:rank].copy()
    history = []
    previous = np.inf
    for _ in range(max_iterations):
        for i in range(effect.shape[0]):
            m = observed[i]
            scores[i] = solve_wls(loadings[:, m].T, effect[i, m] - mean[m], weights[i, m], ridge)
        design = np.column_stack([np.ones(effect.shape[0]), scores])
        for j in range(effect.shape[1]):
            m = observed[:, j]
            coef = solve_wls(design[m], effect[m, j], weights[m, j], ridge, False)
            mean[j], loadings[:, j] = coef[0], coef[1:]
        scores, loadings = orient(scores, loadings)
        fitted = mean + scores @ loadings
        loss = float(np.sum(weights[observed] * np.square(effect[observed] - fitted[observed])))
        history.append(loss)
        if previous < np.inf and abs(previous - loss) <= 1e-7 * max(previous, 1.0):
            break
        previous = loss
    return Model(mean, scores, loadings, mean + scores @ loadings, history)


def deterministic_mask(observed: np.ndarray, fraction: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    hidden = np.zeros_like(observed, dtype=bool)
    for i in range(observed.shape[0]):
        idx = np.flatnonzero(observed[i])
        n = max(1, int(round(len(idx) * fraction)))
        hidden[i, rng.choice(idx, n, replace=False)] = True
    return hidden


def project(model: Model, effect: np.ndarray, uncertainty: np.ndarray, anchors: np.ndarray,
            use_uncertainty: bool, ridge: float) -> np.ndarray:
    scores = np.full((effect.shape[0], model.loadings.shape[0]), np.nan)
    for i in range(effect.shape[0]):
        m = anchors[i] & np.isfinite(effect[i]) & np.isfinite(uncertainty[i]) & (uncertainty[i] > 0)
        if m.sum() < model.loadings.shape[0] + 2:
            continue
        raw = 1.0 / np.square(uncertainty[i, m]) if use_uncertainty else np.ones(m.sum())
        lo, hi = np.quantile(raw, [0.05, 0.95])
        w = np.clip(raw, lo, hi)
        w /= w.mean()
        scores[i] = solve_wls(model.loadings[:, m].T, effect[i, m] - model.mean[m], w, ridge)
    return scores


def metric_components(y: np.ndarray, pred: np.ndarray, baseline: np.ndarray, u: np.ndarray,
                      hidden: np.ndarray, groups: np.ndarray) -> pd.DataFrame:
    weights = objective_weights(y, u, groups, True)
    rows = []
    for i in range(y.shape[0]):
        m = hidden[i] & np.isfinite(y[i]) & np.isfinite(pred[i]) & (weights[i] > 0)
        w = weights[i, m]
        rows.append({
            "row": i, "heldout_cells": int(m.sum()), "weight_sum": float(w.sum()),
            "sse": float(np.sum(w * np.square(y[i, m] - pred[i, m]))),
            "sst": float(np.sum(w * np.square(y[i, m] - baseline[m]))),
        })
    return pd.DataFrame(rows)


def inner_select_rank(effect: np.ndarray, uncertainty: np.ndarray, meta: pd.DataFrame,
                      outer_train: np.ndarray, args, outer_id: int, fold_dir: Path):
    groups_all = meta.independence_group.to_numpy()
    unique_groups = np.unique(groups_all[outer_train])
    n_splits = min(args.inner_folds, len(unique_groups))
    splitter = GroupKFold(n_splits=n_splits)
    detail = []
    base_rows = np.flatnonzero(outer_train)
    dummy = np.zeros(len(base_rows))
    for inner_fold, (tr_rel, va_rel) in enumerate(splitter.split(dummy, groups=groups_all[base_rows]), 1):
        tr, va = base_rows[tr_rel], base_rows[va_rel]
        selected, _ = select_genes(effect, uncertainty, tr, args.coverage, args.genes)
        e_tr, u_tr = effect[np.ix_(tr, selected)], uncertainty[np.ix_(tr, selected)]
        e_va, u_va = effect[np.ix_(va, selected)], uncertainty[np.ix_(va, selected)]
        observed = np.isfinite(e_va) & np.isfinite(u_va) & (u_va > 0)
        for k in range(1, args.max_k + 1):
            model = fit_model(e_tr, u_tr, groups_all[tr], k, True, args.max_iterations, args.ridge)
            for repeat in range(args.inner_masks):
                hidden = deterministic_mask(observed, args.mask_fraction,
                                            SEED + outer_id * 10000 + inner_fold * 100 + repeat)
                anchors = observed & ~hidden
                scores = project(model, e_va, u_va, anchors, True, args.ridge)
                pred = model.mean + scores @ model.loadings
                comp = metric_components(e_va, pred, model.mean, u_va, hidden, groups_all[va])
                detail.append({"inner_fold": inner_fold, "mask_repeat": repeat + 1, "K": k,
                               "sse": comp.sse.sum(), "sst": comp.sst.sum(),
                               "heldout_cells": comp.heldout_cells.sum(),
                               "weighted_rmse": math.sqrt(comp.sse.sum() / comp.weight_sum.sum())})
    d = pd.DataFrame(detail)
    agg = d.groupby("K", as_index=False).agg(mean_rmse=("weighted_rmse", "mean"),
                                               sd_rmse=("weighted_rmse", "std"),
                                               evaluations=("weighted_rmse", "size"))
    agg["sem_rmse"] = agg.sd_rmse / np.sqrt(agg.evaluations)
    best = agg.loc[agg.mean_rmse.idxmin()]
    eligible = agg[agg.mean_rmse <= best.mean_rmse + best.sem_rmse]
    chosen = int(eligible.K.min())
    d.to_csv(fold_dir / "inner_rank_detail.tsv", sep="\t", index=False)
    agg.to_csv(fold_dir / "inner_rank_summary.tsv", sep="\t", index=False)
    return chosen, agg


def extend_loadings(model: Model, effect: np.ndarray, uncertainty: np.ndarray,
                    groups: np.ndarray, ridge: float) -> tuple[np.ndarray, np.ndarray]:
    w = objective_weights(effect, uncertainty, groups, True)
    design = np.column_stack([np.ones(effect.shape[0]), model.scores])
    means = np.full(effect.shape[1], np.nan)
    loads = np.full((model.loadings.shape[0], effect.shape[1]), np.nan)
    for j in range(effect.shape[1]):
        m = w[:, j] > 0
        if m.sum() >= model.loadings.shape[0] + 2:
            coef = solve_wls(design[m], effect[m, j], w[m, j], ridge, False)
            means[j], loads[:, j] = coef[0], coef[1:]
    return means, loads


def run_outer_cv(manifest, genes, effect, uncertainty, args, out):
    all_predictions, fold_summaries, outer_ranks, crossfit = [], [], [], []
    groups = manifest.independence_group.to_numpy()
    for outer_id, disease_class in enumerate(sorted(manifest.disease_class.unique()), 1):
        fold_dir = out / "outer_folds" / f"{outer_id:02d}_{slug(disease_class)}"
        fold_dir.mkdir(parents=True, exist_ok=True)
        test = manifest.disease_class.eq(disease_class).to_numpy()
        blocked_accessions = set(manifest.loc[test, "Accession"])
        blocked_groups = set(manifest.loc[test, "independence_group"])
        purged = (~test) & (manifest.Accession.isin(blocked_accessions) | manifest.independence_group.isin(blocked_groups)).to_numpy()
        train = ~(test | purged)
        chosen_k, rank_table = inner_select_rank(effect, uncertainty, manifest, train, args, outer_id, fold_dir)
        selected, selection_audit = select_genes(effect, uncertainty, np.flatnonzero(train), args.coverage, args.genes)
        selection_audit = pd.concat([genes, selection_audit], axis=1)
        selection_audit.to_csv(fold_dir / "gene_selection.tsv.gz", sep="\t", index=False, compression="gzip")
        tr, te = np.flatnonzero(train), np.flatnonzero(test)
        e_tr, u_tr = effect[np.ix_(tr, selected)], uncertainty[np.ix_(tr, selected)]
        e_te, u_te = effect[np.ix_(te, selected)], uncertainty[np.ix_(te, selected)]
        weighted = fit_model(e_tr, u_tr, groups[tr], chosen_k, True, args.max_iterations, args.ridge)
        unweighted = fit_model(e_tr, u_tr, groups[tr], chosen_k, False, args.max_iterations, args.ridge)
        observed = np.isfinite(e_te) & np.isfinite(u_te) & (u_te > 0)
        for repeat in range(args.outer_masks):
            hidden = deterministic_mask(observed, args.mask_fraction, SEED + outer_id * 100 + repeat)
            anchors = observed & ~hidden
            for method, model, weighted_projection in [
                ("uncertainty_weighted_low_rank", weighted, True),
                ("group_balanced_unweighted_low_rank", unweighted, False),
            ]:
                scores = project(model, e_te, u_te, anchors, weighted_projection, args.ridge)
                pred = model.mean + scores @ model.loadings
                comp = metric_components(e_te, pred, weighted.mean, u_te, hidden, groups[te])
                for local_i, global_i in enumerate(te):
                    row = comp.iloc[local_i].to_dict()
                    row.update({"disease_class": disease_class, "mask_repeat": repeat + 1,
                                "method": method, "EU_ID": manifest.EU_ID.iloc[global_i],
                                "Accession": manifest.Accession.iloc[global_i],
                                "independence_group": groups[global_i], "K": chosen_k})
                    all_predictions.append(row)
            base_comp = metric_components(e_te, np.broadcast_to(weighted.mean, e_te.shape), weighted.mean,
                                          u_te, hidden, groups[te])
            for local_i, global_i in enumerate(te):
                row = base_comp.iloc[local_i].to_dict()
                row.update({"disease_class": disease_class, "mask_repeat": repeat + 1,
                            "method": "training_gene_mean", "EU_ID": manifest.EU_ID.iloc[global_i],
                            "Accession": manifest.Accession.iloc[global_i],
                            "independence_group": groups[global_i], "K": 0})
                all_predictions.append(row)
        # Cross-fitted residuals over every gene meeting the outer-training coverage gate.
        coverage_count = np.isfinite(effect[tr]).sum(0)
        eligible = coverage_count >= math.ceil(args.coverage * len(tr))
        residual_anchor_idx = selected[: min(500, len(selected) // 2)]
        eligible[residual_anchor_idx] = False
        ext_idx = np.flatnonzero(eligible)
        ext_mean, ext_load = extend_loadings(weighted, effect[np.ix_(tr, ext_idx)],
                                             uncertainty[np.ix_(tr, ext_idx)], groups[tr], args.ridge)
        residual_anchors = np.zeros_like(observed)
        residual_anchors[:, :len(residual_anchor_idx)] = observed[:, :len(residual_anchor_idx)]
        test_scores = project(weighted, e_te, u_te, residual_anchors, True, args.ridge)
        ext_pred = ext_mean + test_scores @ ext_load
        for local_i, global_i in enumerate(te):
            m = np.isfinite(effect[global_i, ext_idx]) & np.isfinite(ext_pred[local_i])
            for j, pred in zip(ext_idx[m], ext_pred[local_i, m]):
                crossfit.append({"EU_ID": manifest.EU_ID.iloc[global_i], "disease_class": disease_class,
                                 "independence_group": groups[global_i], "gene_index": j,
                                 "effect": effect[global_i, j], "uncertainty": uncertainty[global_i, j],
                                 "shared_prediction": pred, "residual": effect[global_i, j] - pred})
        split = manifest[["EU_ID", "Accession", "independence_group", "disease_class"]].copy()
        split["role"] = np.where(test, "test", np.where(purged, "purged", "train"))
        split.to_csv(fold_dir / "split_manifest.tsv", sep="\t", index=False)
        np.savez_compressed(fold_dir / "model_parameters.npz", selected_gene_indices=selected,
                            gene_mean=weighted.mean, loadings=weighted.loadings, scores=weighted.scores)
        outer_ranks.append({"disease_class": disease_class, "selected_K": chosen_k})
        fold_summaries.append({"disease_class": disease_class, "train_contrasts": int(train.sum()),
                               "test_contrasts": int(test.sum()), "purged_contrasts": int(purged.sum()),
                               "train_groups": int(manifest.loc[train, "independence_group"].nunique()),
                               "test_groups": int(manifest.loc[test, "independence_group"].nunique()),
                               "selected_K": chosen_k})
        print(f"[outer {outer_id:02d}/12] {disease_class}: K={chosen_k}, train={train.sum()}, test={test.sum()}, purged={purged.sum()}", flush=True)
    pred = pd.DataFrame(all_predictions)
    pred.to_csv(out / "outer_prediction_components.tsv", sep="\t", index=False)
    pd.DataFrame(fold_summaries).to_csv(out / "outer_fold_summary.tsv", sep="\t", index=False)
    pd.DataFrame(outer_ranks).to_csv(out / "outer_selected_ranks.tsv", sep="\t", index=False)
    pd.DataFrame(crossfit).to_csv(out / "crossfit_residuals.tsv.gz", sep="\t", index=False, compression="gzip")
    return pred, pd.DataFrame(outer_ranks), pd.DataFrame(crossfit)


def summarize_predictions(pred: pd.DataFrame, repeats: int, bootstrap_repeats: int, out: Path):
    metrics = []
    for (method, repeat), x in pred.groupby(["method", "mask_repeat"]):
        sse, sst, w = x.sse.sum(), x.sst.sum(), x.weight_sum.sum()
        metrics.append({"method": method, "mask_repeat": repeat, "pooled_cv_r2": 1 - sse / sst,
                        "weighted_rmse": math.sqrt(sse / w), "heldout_cells": int(x.heldout_cells.sum())})
    repeat_metrics = pd.DataFrame(metrics)
    repeat_metrics.to_csv(out / "outer_repeat_metrics.tsv", sep="\t", index=False)
    class_metrics = []
    for (method, cls), x in pred.groupby(["method", "disease_class"]):
        class_metrics.append({"method": method, "disease_class": cls,
                              "cv_r2": 1 - x.sse.sum() / x.sst.sum(),
                              "weighted_rmse": math.sqrt(x.sse.sum() / x.weight_sum.sum()),
                              "contrasts": x.EU_ID.nunique(), "groups": x.independence_group.nunique()})
    class_metrics = pd.DataFrame(class_metrics)
    class_metrics.to_csv(out / "outer_class_metrics.tsv", sep="\t", index=False)
    rng = np.random.default_rng(SEED + 900000)
    primary = pred[pred.method.eq("uncertainty_weighted_low_rank")]
    group_comp = primary.groupby(["independence_group", "mask_repeat"], as_index=False)[["sse", "sst"]].sum()
    group_ids = group_comp.independence_group.unique()
    boots = []
    for b in range(bootstrap_repeats):
        sampled = rng.choice(group_ids, len(group_ids), replace=True)
        sse = sst = 0.0
        for g in sampled:
            z = group_comp[group_comp.independence_group.eq(g)]
            sse += z.sse.sum(); sst += z.sst.sum()
        boots.append(1 - sse / sst)
    primary_r2 = 1 - primary.sse.sum() / primary.sst.sum()
    summary = {
        "primary_metric": "purged_leave_one_disease_class_out_pooled_weighted_cv_r2",
        "primary_r2": primary_r2,
        "cluster_bootstrap_ci95": [float(np.quantile(boots, .025)), float(np.quantile(boots, .975))],
        "bootstrap_repeats": bootstrap_repeats,
        "outer_mask_repeats": repeats,
    }
    (out / "primary_result.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return repeat_metrics, class_metrics, summary


def run_accession_sensitivity(manifest, genes, effect, uncertainty, rank, args, out):
    groups = manifest.independence_group.to_numpy()
    rows = []
    accessions = sorted(manifest.Accession.unique())
    if args.quick:
        accessions = accessions[:3]
    for fold, accession in enumerate(accessions, 1):
        test = manifest.Accession.eq(accession).to_numpy()
        blocked = set(manifest.loc[test, "independence_group"])
        train = ~manifest.independence_group.isin(blocked).to_numpy()
        tr, te = np.flatnonzero(train), np.flatnonzero(test)
        selected, _ = select_genes(effect, uncertainty, tr, args.coverage, args.genes)
        e_tr, u_tr = effect[np.ix_(tr, selected)], uncertainty[np.ix_(tr, selected)]
        e_te, u_te = effect[np.ix_(te, selected)], uncertainty[np.ix_(te, selected)]
        model = fit_model(e_tr, u_tr, groups[tr], rank, True, args.max_iterations, args.ridge)
        observed = np.isfinite(e_te) & np.isfinite(u_te) & (u_te > 0)
        for repeat in range(3 if not args.quick else 1):
            hidden = deterministic_mask(observed, args.mask_fraction, SEED + 500000 + fold * 10 + repeat)
            scores = project(model, e_te, u_te, observed & ~hidden, True, args.ridge)
            pred = model.mean + scores @ model.loadings
            comp = metric_components(e_te, pred, model.mean, u_te, hidden, groups[te])
            rows.append({"Accession": accession, "mask_repeat": repeat + 1, "sse": comp.sse.sum(),
                         "sst": comp.sst.sum(), "weight_sum": comp.weight_sum.sum(),
                         "heldout_cells": comp.heldout_cells.sum()})
        print(f"[accession {fold:02d}/{len(accessions)}] {accession}", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(out / "leave_one_accession_out_components.tsv", sep="\t", index=False)
    summary = {"rank": rank, "accessions": len(accessions),
               "pooled_cv_r2": float(1 - d.sse.sum() / d.sst.sum()),
               "weighted_rmse": float(math.sqrt(d.sse.sum() / d.weight_sum.sum()))}
    (out / "leave_one_accession_out_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def align_loadings(reference: np.ndarray, candidate: np.ndarray) -> np.ndarray:
    cross = candidate @ reference.T
    u, _, vt = np.linalg.svd(cross)
    return (u @ vt).T @ candidate


def fit_final_and_stability(manifest, genes, effect, uncertainty, rank, args, out):
    rows = np.arange(len(manifest))
    selected, audit = select_genes(effect, uncertainty, rows, args.coverage, args.genes)
    audit = pd.concat([genes, audit], axis=1)
    audit.to_csv(out / "final_gene_selection.tsv.gz", sep="\t", index=False, compression="gzip")
    groups = manifest.independence_group.to_numpy()
    model = fit_model(effect[:, selected], uncertainty[:, selected], groups, rank, True,
                      args.max_iterations, args.ridge)
    np.savez_compressed(out / "final_model_parameters.npz", selected_gene_indices=selected,
                        gene_mean=model.mean, scores=model.scores, loadings=model.loadings,
                        objective_history=np.asarray(model.history))
    scores = manifest[["EU_ID", "Accession", "independence_group", "disease_class"]].copy()
    for k in range(rank): scores[f"factor_{k+1}"] = model.scores[:, k]
    scores.to_csv(out / "final_contrast_factor_scores.tsv", sep="\t", index=False)
    ext_idx = np.flatnonzero(audit.eligible.to_numpy())
    ext_mean, ext_load = extend_loadings(model, effect[:, ext_idx], uncertainty[:, ext_idx], groups, args.ridge)
    load_table = genes.iloc[ext_idx].copy()
    load_table["gene_mean"] = ext_mean
    load_table["coverage_count"] = np.isfinite(effect[:, ext_idx]).sum(0)
    for k in range(rank): load_table[f"factor_{k+1}_loading"] = ext_load[k]
    load_table.to_csv(out / "extended_gene_loadings.tsv.gz", sep="\t", index=False, compression="gzip")
    rng = np.random.default_rng(SEED + 700000)
    unique_groups = np.unique(groups)
    stability = []
    reps = args.stability_repeats
    for b in range(reps):
        sampled = rng.choice(unique_groups, len(unique_groups), replace=True)
        idx, boot_groups = [], []
        for instance, g in enumerate(sampled):
            take = np.flatnonzero(groups == g)
            idx.extend(take.tolist()); boot_groups.extend([f"b{instance}"] * len(take))
        boot = fit_model(effect[np.ix_(idx, selected)], uncertainty[np.ix_(idx, selected)],
                         np.asarray(boot_groups), rank, True, args.max_iterations, args.ridge)
        aligned = align_loadings(model.loadings, boot.loadings)
        for k in range(rank):
            cosine = np.dot(model.loadings[k], aligned[k]) / (np.linalg.norm(model.loadings[k]) * np.linalg.norm(aligned[k]))
            stability.append({"bootstrap_repeat": b + 1, "factor": k + 1, "cosine_similarity": cosine})
    stability = pd.DataFrame(stability)
    stability.to_csv(out / "factor_stability_bootstrap.tsv", sep="\t", index=False)
    return model, load_table, stability


def random_effect_meta(values: np.ndarray, se: np.ndarray) -> tuple[float, float, float, int]:
    good = np.isfinite(values) & np.isfinite(se) & (se > 0)
    y, s = values[good], se[good]
    if len(y) < 2: return np.nan, np.nan, np.nan, len(y)
    w = 1 / s**2; fixed = np.sum(w*y)/np.sum(w); q = np.sum(w*(y-fixed)**2)
    c = np.sum(w) - np.sum(w**2)/np.sum(w)
    tau2 = max(0.0, (q-(len(y)-1))/c) if c > 0 else 0.0
    wr = 1/(s**2+tau2); est = np.sum(wr*y)/np.sum(wr); out_se = math.sqrt(1/np.sum(wr))
    return est, out_se, tau2, len(y)


def residual_meta(crossfit: pd.DataFrame, genes: pd.DataFrame, out: Path):
    # First collapse non-independent contrasts inside each independence group.
    z = crossfit.copy()
    z["ivw"] = 1 / np.square(z.uncertainty)
    z["weighted_residual"] = z.ivw * z.residual
    grouped = z.groupby(
        ["disease_class", "independence_group", "gene_index"], as_index=False
    ).agg(ivw=("ivw", "sum"), weighted_residual=("weighted_residual", "sum"))
    grouped["residual"] = grouped.weighted_residual / grouped.ivw
    grouped["se"] = np.sqrt(1 / grouped.ivw)
    grouped["w"] = 1 / np.square(grouped.se)
    grouped["wy"] = grouped.w * grouped.residual
    grouped["wyy"] = grouped.w * np.square(grouped.residual)
    grouped["w2"] = np.square(grouped.w)
    keys = ["disease_class", "gene_index"]
    first = grouped.groupby(keys, as_index=False).agg(
        sum_w=("w", "sum"), sum_wy=("wy", "sum"), sum_wyy=("wyy", "sum"),
        sum_w2=("w2", "sum"), independence_groups=("residual", "size")
    )
    first["fixed"] = first.sum_wy / first.sum_w
    first["q"] = first.sum_wyy - np.square(first.sum_wy) / first.sum_w
    first["c"] = first.sum_w - first.sum_w2 / first.sum_w
    first["tau2"] = np.maximum(
        0.0,
        np.divide(first.q - (first.independence_groups - 1), first.c,
                  out=np.zeros(len(first)), where=first.c.to_numpy() > 0),
    )
    grouped = grouped.merge(first[keys + ["tau2"]], on=keys, validate="many_to_one")
    grouped["wr"] = 1 / (np.square(grouped.se) + grouped.tau2)
    grouped["wry"] = grouped.wr * grouped.residual
    second = grouped.groupby(keys, as_index=False).agg(sum_wr=("wr", "sum"), sum_wry=("wry", "sum"))
    result = first.merge(second, on=keys, validate="one_to_one")
    result["residual_effect"] = result.sum_wry / result.sum_wr
    result["meta_se"] = np.sqrt(1 / result.sum_wr)
    result["p_value"] = 2 * norm.sf(np.abs(result.residual_effect / result.meta_se))
    result.loc[result.independence_groups.lt(2), ["residual_effect", "meta_se", "tau2", "p_value"]] = np.nan
    result = result[keys + ["residual_effect", "meta_se", "tau2", "independence_groups", "p_value"]]
    result["adjusted_p_value"] = result.groupby("disease_class").p_value.transform(lambda x: bh(x.to_numpy()))
    result = result.merge(genes.reset_index().rename(columns={"index": "gene_index"}), on="gene_index", how="left")
    result.to_csv(out / "disease_class_crossfit_residual_meta.tsv.gz", sep="\t", index=False, compression="gzip")
    return result


def fetch_gmt(name: str, path: Path) -> dict:
    url = f"https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={name}"
    if not path.exists():
        req = urllib.request.Request(url, headers={"User-Agent": "EU57-reproducible-analysis/1.0"})
        with urllib.request.urlopen(req, timeout=120) as response:
            data = response.read()
        path.write_bytes(data)
    return {
        "library": name,
        "url": url,
        # The manifest is stored beside the GMT files, so a file name is
        # sufficient and remains valid after the package is moved.
        "path": path.name,
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
    }


def parse_gmt(path: Path) -> dict[str, set[str]]:
    sets = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3:
            sets[parts[0]] = {x.upper() for x in parts[2:] if x}
    return sets


def rank_enrichment(table: pd.DataFrame, score_col: str, gene_sets: dict[str, set[str]],
                    analysis: str, min_size=10, max_size=500) -> pd.DataFrame:
    x = table[["canonical_symbol", score_col]].dropna().drop_duplicates("canonical_symbol")
    x = x[x.canonical_symbol.ne("")]
    symbols = x.canonical_symbol.str.upper().to_numpy(); scores = x[score_col].to_numpy(float)
    symbol_index = {symbol: i for i, symbol in enumerate(symbols)}
    ranks = rankdata(scores, method="average")
    _, tie_counts = np.unique(scores, return_counts=True)
    tie_term = np.sum(tie_counts**3 - tie_counts) / (len(scores) * max(len(scores) - 1, 1))
    rows = []
    for term, members in gene_sets.items():
        idx = np.fromiter((symbol_index[g] for g in members if g in symbol_index), dtype=int)
        n, other = len(idx), len(scores) - len(idx)
        if n < min_size or n > max_size or other < 2: continue
        stat = ranks[idx].sum() - n * (n + 1) / 2
        variance = n * other / 12 * ((len(scores) + 1) - tie_term)
        z = (stat - n * other / 2) / math.sqrt(variance) if variance > 0 else 0.0
        p = 2 * norm.sf(abs(z))
        auc = stat / (n * other)
        rows.append({"analysis": analysis, "term": term, "set_size_in_universe": n,
                     "rank_auc": auc, "signed_enrichment": 2*(auc-.5), "p_value": p})
    out = pd.DataFrame(rows)
    if len(out): out["adjusted_p_value"] = bh(out.p_value.to_numpy())
    return out.sort_values(["adjusted_p_value", "p_value"]) if len(out) else out


def run_enrichment(loadings, residuals, rank, args, out):
    annotation_dir = out / "annotation_sources"; annotation_dir.mkdir(exist_ok=True)
    manifests, libraries = [], {}
    for name in ["MSigDB_Hallmark_2020", "Reactome_2022"]:
        path = annotation_dir / f"{name}.gmt"
        manifests.append(fetch_gmt(name, path)); libraries[name] = parse_gmt(path)
    (annotation_dir / "annotation_manifest.json").write_text(json.dumps(manifests, indent=2), encoding="utf-8")
    enriched = []
    for k in range(rank):
        for lib, sets in libraries.items():
            enriched.append(rank_enrichment(loadings, f"factor_{k+1}_loading", sets, f"factor_{k+1}|{lib}"))
    for cls, x in residuals.groupby("disease_class"):
        for lib, sets in libraries.items():
            enriched.append(rank_enrichment(x, "residual_effect", sets, f"residual:{cls}|{lib}"))
    result = pd.concat(enriched, ignore_index=True)
    result.to_csv(out / "pathway_rank_enrichment.tsv", sep="\t", index=False)
    return result


def plots(class_metrics, repeat_metrics, ranks, stability, enrichment, out):
    figdir = out / "figures"; figdir.mkdir(exist_ok=True)
    primary = class_metrics[class_metrics.method.eq("uncertainty_weighted_low_rank")].sort_values("cv_r2")
    fig, ax = plt.subplots(figsize=(9, 6)); ax.barh(primary.disease_class, primary.cv_r2)
    ax.axvline(0, color="black", lw=.8); ax.set_xlabel("Weighted cross-validated R²")
    fig.tight_layout(); fig.savefig(figdir / "class_cv_r2.png", dpi=180); plt.close(fig)
    order = ["training_gene_mean", "group_balanced_unweighted_low_rank", "uncertainty_weighted_low_rank"]
    data = [repeat_metrics.loc[repeat_metrics.method.eq(m), "pooled_cv_r2"] for m in order]
    fig, ax = plt.subplots(figsize=(8, 5)); ax.boxplot(data, tick_labels=["Gene mean", "Unweighted", "Uncertainty-weighted"])
    ax.axhline(0, color="black", lw=.8); ax.set_ylabel("Pooled weighted CV R²")
    fig.tight_layout(); fig.savefig(figdir / "method_comparison.png", dpi=180); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4)); counts = ranks.selected_K.value_counts().sort_index()
    ax.bar(counts.index.astype(str), counts.values); ax.set_xlabel("Nested-CV selected K"); ax.set_ylabel("Outer folds")
    fig.tight_layout(); fig.savefig(figdir / "selected_rank_distribution.png", dpi=180); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4)); stability.boxplot(column="cosine_similarity", by="factor", ax=ax)
    ax.set_title(""); fig.suptitle(""); ax.set_xlabel("Factor"); ax.set_ylabel("Bootstrap loading cosine similarity")
    fig.tight_layout(); fig.savefig(figdir / "factor_stability.png", dpi=180); plt.close(fig)


def artifact_manifest(out: Path):
    rows = []
    for path in sorted(p for p in out.rglob("*") if p.is_file() and p.name != "artifact_manifest.tsv"):
        rows.append({
            "path": path.relative_to(out).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    pd.DataFrame(rows).to_csv(out / "artifact_manifest.tsv", sep="\t", index=False)


def main():
    args = parse_args()
    if args.quick:
        args.inner_folds, args.inner_masks, args.outer_masks = 2, 1, 1
        args.max_k, args.max_iterations = 2, 12
        args.bootstrap_repeats, args.stability_repeats = 20, 5
    out = args.output.resolve()
    smoke_marker = out / "SMOKE_TEST_ONLY.txt"
    if not args.quick and smoke_marker.exists():
        raise RuntimeError(f"Refusing to write formal results into a smoke-test directory: {out}")
    if args.quick and out == DEFAULT_OUTPUT.resolve():
        raise RuntimeError("Quick mode must use a separate --output directory; formal results are protected")
    out.mkdir(parents=True, exist_ok=True)
    if args.quick:
        smoke_marker.write_text(
            "SMOKE TEST ONLY. These artifacts are incomplete, use reduced settings, and must never be reported as formal results.\n",
            encoding="utf-8",
        )
    started = time.time()
    manifest, genes, effect, uncertainty = load_inputs(args.manifest.resolve(), args.audit.resolve(), args.class_map.resolve())
    manifest.to_csv(out / "analysis_manifest.tsv", sep="\t", index=False)
    genes.to_csv(out / "global_gene_index.tsv.gz", sep="\t", index=False, compression="gzip")
    np.savez_compressed(out / "global_eu_matrices.npz", effect=effect, uncertainty=uncertainty)
    pred, ranks, crossfit = run_outer_cv(manifest, genes, effect, uncertainty, args, out)
    repeat_metrics, class_metrics, primary = summarize_predictions(pred, args.outer_masks, args.bootstrap_repeats, out)
    final_rank_dir = out / "final_rank_selection"
    final_rank_dir.mkdir(exist_ok=True)
    final_rank, final_rank_table = inner_select_rank(
        effect, uncertainty, manifest, np.ones(len(manifest), dtype=bool), args, 99, final_rank_dir
    )
    (final_rank_dir / "selected_rank.json").write_text(
        json.dumps({"selected_K": final_rank, "selection_rule": "one_standard_error"}, indent=2),
        encoding="utf-8",
    )
    accession = run_accession_sensitivity(manifest, genes, effect, uncertainty, final_rank, args, out)
    model, loadings, stability = fit_final_and_stability(manifest, genes, effect, uncertainty, final_rank, args, out)
    residuals = residual_meta(crossfit, genes, out)
    enrichment = run_enrichment(loadings, residuals, final_rank, args, out)
    plots(class_metrics, repeat_metrics, ranks, stability, enrichment, out)
    summary = {
        "schema": "eu57_cross_disease_cv_20260828_v1", "seed": SEED,
        "n_contrasts": len(manifest), "n_accessions": int(manifest.Accession.nunique()),
        "n_independence_groups": int(manifest.independence_group.nunique()),
        "n_disease_classes": int(manifest.disease_class.nunique()), "final_rank": final_rank,
        "primary_result": primary, "accession_sensitivity": accession,
        "median_factor_stability": stability.groupby("factor").cosine_similarity.median().to_dict(),
        "runtime_seconds": time.time() - started, "quick_mode": args.quick,
        "software": {"python": sys.version, "platform": platform.platform(), "numpy": np.__version__,
                     "pandas": pd.__version__},
        "parameters": vars(args),
    }
    # Store logical package locations rather than machine-specific absolute
    # paths. The output is represented by "." because the summary lives there.
    summary["parameters"] = {
        k: ("." if k == "output" else package_relative(v)) if isinstance(v, Path) else v
        for k, v in summary["parameters"].items()
    }
    (out / "experiment_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    artifact_manifest(out)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
