"""Training-fold-only 500/1000/2000-gene sensitivity at the selected final rank."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
CORE_PATH = HERE / "run_experiment.py"
SPEC = importlib.util.spec_from_file_location("eu57_core", CORE_PATH)
core = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = core
SPEC.loader.exec_module(core)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Repeat disease-domain validation across three training-selected gene panels."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=core.DEFAULT_OUTPUT,
        help="Formal result directory containing final_rank_selection/selected_rank.json.",
    )
    args = parser.parse_args()
    out = args.output.resolve()
    if (out / "SMOKE_TEST_ONLY.txt").exists():
        raise RuntimeError("Panel sensitivity for formal analysis cannot use a smoke-test directory")
    manifest, genes, effect, uncertainty = core.load_inputs(
        core.DEFAULT_MANIFEST, core.DEFAULT_AUDIT, core.DEFAULT_CLASS_MAP
    )
    groups = manifest.independence_group.to_numpy()
    rank_info = json.loads((out / "final_rank_selection" / "selected_rank.json").read_text(encoding="utf-8"))
    fixed_rank = int(rank_info["selected_K"])
    components = []
    for panel_size in (500, 1000, 2000):
        for fold, disease_class in enumerate(sorted(manifest.disease_class.unique()), 1):
            test = manifest.disease_class.eq(disease_class).to_numpy()
            blocked_a = set(manifest.loc[test, "Accession"])
            blocked_g = set(manifest.loc[test, "independence_group"])
            purge = (~test) & (
                manifest.Accession.isin(blocked_a) | manifest.independence_group.isin(blocked_g)
            ).to_numpy()
            train = ~(test | purge)
            tr, te = np.flatnonzero(train), np.flatnonzero(test)
            selected, _ = core.select_genes(effect, uncertainty, tr, 0.90, panel_size)
            e_tr, u_tr = effect[np.ix_(tr, selected)], uncertainty[np.ix_(tr, selected)]
            e_te, u_te = effect[np.ix_(te, selected)], uncertainty[np.ix_(te, selected)]
            model = core.fit_model(e_tr, u_tr, groups[tr], fixed_rank, True, 80, 1e-3)
            observed = np.isfinite(e_te) & np.isfinite(u_te) & (u_te > 0)
            for repeat in range(5):
                hidden = core.deterministic_mask(
                    observed, 0.20, core.SEED + panel_size * 1000 + fold * 100 + repeat
                )
                scores = core.project(model, e_te, u_te, observed & ~hidden, True, 1e-3)
                pred = model.mean + scores @ model.loadings
                comp = core.metric_components(e_te, pred, model.mean, u_te, hidden, groups[te])
                components.append({
                    "panel_size": panel_size, "disease_class": disease_class,
                    "mask_repeat": repeat + 1, "sse": comp.sse.sum(), "sst": comp.sst.sum(),
                    "weight_sum": comp.weight_sum.sum(), "heldout_cells": comp.heldout_cells.sum(),
                })
            print(f"[panel={panel_size} fold={fold:02d}/12] {disease_class}", flush=True)
    detail = pd.DataFrame(components)
    detail.to_csv(out / "panel_size_sensitivity_components.tsv", sep="\t", index=False)
    rows = []
    for panel_size, x in detail.groupby("panel_size"):
        rows.append({
            "panel_size": panel_size, "fixed_rank": fixed_rank,
            "pooled_weighted_cv_r2": 1 - x.sse.sum() / x.sst.sum(),
            "weighted_rmse": math.sqrt(x.sse.sum() / x.weight_sum.sum()),
            "heldout_cells": int(x.heldout_cells.sum()),
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "panel_size_sensitivity.tsv", sep="\t", index=False)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
