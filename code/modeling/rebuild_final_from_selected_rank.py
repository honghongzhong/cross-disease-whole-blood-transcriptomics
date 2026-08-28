"""Rebuild post-validation artifacts after the all-data grouped CV rank is frozen."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("eu57_core", HERE / "run_experiment.py")
core = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = core
SPEC.loader.exec_module(core)
OUT = core.DEFAULT_OUTPUT


def main() -> None:
    selected = json.loads((OUT / "final_rank_selection/selected_rank.json").read_text(encoding="utf-8"))
    rank = int(selected["selected_K"])
    args = SimpleNamespace(
        coverage=0.90, genes=1000, max_iterations=80, ridge=1e-3,
        stability_repeats=200, quick=False, mask_fraction=0.20,
    )
    manifest, genes, effect, uncertainty = core.load_inputs(
        core.DEFAULT_MANIFEST, core.DEFAULT_AUDIT, core.DEFAULT_CLASS_MAP
    )
    accession = core.run_accession_sensitivity(manifest, genes, effect, uncertainty, rank, args, OUT)
    _, loadings, stability = core.fit_final_and_stability(
        manifest, genes, effect, uncertainty, rank, args, OUT
    )
    residuals = pd.read_csv(OUT / "disease_class_crossfit_residual_meta.tsv.gz", sep="\t")
    enrichment = core.run_enrichment(loadings, residuals, rank, args, OUT)
    class_metrics = pd.read_csv(OUT / "outer_class_metrics.tsv", sep="\t")
    repeat_metrics = pd.read_csv(OUT / "outer_repeat_metrics.tsv", sep="\t")
    ranks = pd.read_csv(OUT / "outer_selected_ranks.tsv", sep="\t")
    core.plots(class_metrics, repeat_metrics, ranks, stability, enrichment, OUT)
    summary_path = OUT / "experiment_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["final_rank"] = rank
    summary["final_rank_selection"] = selected
    summary["accession_sensitivity"] = accession
    summary["median_factor_stability"] = {
        str(int(k)): float(v) for k, v in stability.groupby("factor").cosine_similarity.median().items()
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    core.artifact_manifest(OUT)
    print(json.dumps({"final_rank": rank, "accession_sensitivity": accession,
                      "median_factor_stability": summary["median_factor_stability"]}, indent=2))


if __name__ == "__main__":
    main()
