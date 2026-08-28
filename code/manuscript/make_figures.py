"""Generate publication-quality figures from the formal EU57 result tables."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parents[1]
RESULTS = PACKAGE_ROOT / "reference_results"
FIGURES = PACKAGE_ROOT / "manuscript_figures"
FIGURES.mkdir(parents=True, exist_ok=True)

BLUE = "#2878B5"
ORANGE = "#E07A3F"
GREEN = "#3A9D77"
GREY = "#73808C"
LIGHT = "#EEF4F8"

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def box(ax, xy, width, height, title, detail, color=BLUE):
    x, y = xy
    patch = FancyBboxPatch(
        (x, y), width, height,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=1.3, edgecolor=color, facecolor="white",
    )
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height * 0.70, title, ha="center", va="center", weight="bold", color=color, fontsize=7.7, linespacing=1.05)
    ax.text(x + width / 2, y + height * 0.28, detail, ha="center", va="center", color="#26343F", fontsize=6.3, linespacing=1.18)


def arrow(ax, start, end):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=12, lw=1.2, color=GREY))


def figure1_design() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    box(ax, (0.03, 0.73), 0.25, 0.18, "Standardized\ninputs", "57 case-control contrasts\n49 GEO accessions\n48 independent sources")
    box(ax, (0.375, 0.73), 0.25, 0.18, "Disease-aware\ntesting", "Hold out one disease domain\nRemove all source-linked contrasts\nfrom model training")
    box(ax, (0.72, 0.73), 0.25, 0.18, "Training-only\nmodel selection", "Select covered, variable genes\nChoose 1–6 factors by grouped\ntraining-data validation")
    arrow(ax, (0.285, 0.82), (0.37, 0.82))
    arrow(ax, (0.63, 0.82), (0.715, 0.82))

    box(ax, (0.12, 0.42), 0.30, 0.18, "Fit shared\ntranscriptional structure", "Uncertainty-weighted low-rank model\nEqual total weight for each\nindependent source")
    box(ax, (0.58, 0.42), 0.30, 0.18, "Predict unseen\ngene effects", "Estimate test scores from 80%\nPredict the hidden 20%\nRepeat five fixed masks")
    arrow(ax, (0.84, 0.72), (0.42, 0.605))
    arrow(ax, (0.425, 0.51), (0.575, 0.51))

    box(ax, (0.03, 0.09), 0.28, 0.19, "Primary\nevaluation", "Pool hidden-gene errors\nCompare with the training mean\nBootstrap independent sources")
    box(ax, (0.36, 0.09), 0.28, 0.19, "Robustness\nanalyses", "Leave one accession out\nRepeat with three panel sizes\nResample sources for stability")
    box(ax, (0.69, 0.09), 0.28, 0.19, "Biological\ninterpretation", "Fit the final five-factor model\nRank-based pathway enrichment\nAnalyse cross-fitted residuals")
    arrow(ax, (0.72, 0.415), (0.18, 0.285))
    arrow(ax, (0.72, 0.415), (0.50, 0.285))
    arrow(ax, (0.72, 0.415), (0.83, 0.285))
    save(fig, "Figure_1_study_design")


def figure2_performance() -> None:
    repeats = pd.read_csv(RESULTS / "outer_repeat_metrics.tsv", sep="\t")
    classes = pd.read_csv(RESULTS / "outer_class_metrics.tsv", sep="\t")
    names = ["Gene mean", "Unweighted\nlow-rank", "Uncertainty-weighted\nlow-rank"]
    methods = ["training_gene_mean", "group_balanced_unweighted_low_rank", "uncertainty_weighted_low_rank"]
    colors = [GREY, ORANGE, BLUE]
    pooled = []
    for method in methods:
        x = repeats[repeats.method.eq(method)]
        pooled.append((x.pooled_cv_r2.mean(), x.weighted_rmse.mean(), x.pooled_cv_r2.to_numpy()))

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.3), gridspec_kw={"height_ratios": [0.72, 1.55]})
    ax = axes[0]
    ypos = np.arange(3)[::-1]
    ax.barh(ypos, [x[0] for x in pooled], color=colors, height=0.58, zorder=2)
    for i, (_, _, points) in enumerate(pooled):
        jitter = np.linspace(-0.10, 0.10, len(points))
        ax.scatter(points, np.full(len(points), ypos[i]) + jitter, color="white", edgecolor="#26343F", s=22, zorder=3, linewidth=0.7)
        label_x = max(points.max() + 0.018, 0.018)
        ax.text(label_x, ypos[i], f"{pooled[i][0]:.3f}", va="center", fontsize=7.2)
    ax.axvline(0, color="#26343F", lw=0.8)
    ax.set_yticks(ypos, [x.replace("\n", " ") for x in names])
    ax.set_xlabel("Pooled weighted cross-validated $R^2$")
    ax.set_xlim(-0.03, 0.64)
    ax.set_title("Model comparison")
    ax.text(-0.08, 1.06, "a", transform=ax.transAxes, fontsize=12, weight="bold")

    ax = axes[1]
    x = classes[classes.method.eq("uncertainty_weighted_low_rank")].sort_values("cv_r2")
    replacements = {
        "Cardiometabolic and vascular disease": "Cardiometabolic & vascular",
        "Chronic respiratory disease": "Chronic respiratory",
        "Endocrine autoimmunity": "Endocrine autoimmunity",
        "Genetic and congenital disease": "Genetic & congenital",
        "Hematologic neoplasm": "Hematologic neoplasm",
        "Infectious disease": "Infectious",
        "Inflammatory bowel disease": "Inflammatory bowel",
        "Neurologic disease": "Neurologic",
        "Primary immunodeficiency": "Primary immunodeficiency",
        "Solid cancer": "Solid cancer",
        "Systemic autoimmune and rheumatic disease": "Systemic autoimmune & rheumatic",
        "Transplantation and acute tissue injury": "Transplantation & acute injury",
    }
    labels = x.disease_class.map(replacements)
    y = np.arange(len(x))
    ax.barh(y, x.cv_r2, color=BLUE, alpha=0.9)
    ax.set_yticks(y, labels)
    ax.set_xlim(0, 0.9)
    ax.set_xlabel("Weighted cross-validated $R^2$")
    ax.set_title("Performance by held-out disease domain")
    for yi, value, groups in zip(y, x.cv_r2, x.groups):
        ax.text(value + 0.015, yi, f"{value:.3f}  (g={int(groups)})", va="center", fontsize=7.2)
    ax.text(-0.08, 1.03, "b", transform=ax.transAxes, fontsize=12, weight="bold")
    fig.subplots_adjust(left=0.34, right=0.96, bottom=0.08, top=0.95, hspace=0.48)
    save(fig, "Figure_2_predictive_performance")


def figure3_robustness() -> None:
    ranks = pd.read_csv(RESULTS / "final_rank_selection/inner_rank_summary.tsv", sep="\t")
    panel = pd.read_csv(RESULTS / "panel_size_sensitivity.tsv", sep="\t")
    primary = json.loads((RESULTS / "primary_result.json").read_text(encoding="utf-8"))
    accession = json.loads((RESULTS / "leave_one_accession_out_summary.json").read_text(encoding="utf-8"))

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.5), gridspec_kw={"height_ratios": [1.0, 1.15]})
    ax = axes[0]
    ax.errorbar(ranks.K, ranks.mean_rmse, yerr=ranks.sem_rmse, color=BLUE, marker="o", capsize=3, lw=1.5)
    best = ranks.loc[ranks.mean_rmse.idxmin()]
    threshold = best.mean_rmse + best.sem_rmse
    ax.axhspan(best.mean_rmse, threshold, color=GREEN, alpha=0.15, label="One-SE range")
    ax.scatter([5], [ranks.loc[ranks.K.eq(5), "mean_rmse"].iloc[0]], s=70, facecolor="white", edgecolor=ORANGE, lw=2, zorder=4, label="Selected K = 5")
    ax.set_xlabel("Number of factors (K)")
    ax.set_ylabel("Grouped validation RMSE")
    ax.set_xticks(ranks.K)
    ax.set_title("Final model rank selection")
    ax.legend(frameon=False, loc="upper right")
    ax.text(-0.08, 1.06, "a", transform=ax.transAxes, fontsize=12, weight="bold")

    ax = axes[1]
    labels = ["Primary disease-domain test", "Leave-one-accession-out"] + [f"Fixed K=5, {int(n):,} genes" for n in panel.panel_size]
    values = [primary["primary_r2"], accession["pooled_cv_r2"], *panel.pooled_weighted_cv_r2]
    y = np.arange(len(labels))[::-1]
    ax.hlines(y, 0, 0.70, color="#E1E7EC", lw=0.8, zorder=0)
    ax.scatter(values, y, color=[BLUE, GREEN, ORANGE, ORANGE, ORANGE], s=42, zorder=3)
    ci = primary["cluster_bootstrap_ci95"]
    ax.errorbar(values[0], y[0], xerr=[[values[0] - ci[0]], [ci[1] - values[0]]], fmt="none", ecolor="#26343F", capsize=3, lw=1.2)
    ax.text(ci[0], y[0] - 0.24, "Source-bootstrap 95% CI", ha="left", va="top", fontsize=6.8, color="#26343F")
    ax.axvline(0, color="#26343F", lw=0.8)
    ax.set_yticks(y, labels)
    ax.set_xlim(0, 0.72)
    ax.set_xlabel("Pooled weighted cross-validated $R^2$")
    ax.set_title("Robustness of reconstruction performance")
    for yi, value in zip(y, values):
        ax.text(value + 0.015, yi, f"{value:.3f}", va="center", fontsize=7.5)
    ax.text(-0.08, 1.04, "b", transform=ax.transAxes, fontsize=12, weight="bold")
    fig.subplots_adjust(left=0.34, right=0.96, bottom=0.09, top=0.94, hspace=0.55)
    save(fig, "Figure_3_model_robustness")


def hallmark_matrix(enrichment: pd.DataFrame, analyses: list[str], terms: list[str]) -> tuple[np.ndarray, np.ndarray]:
    values = np.full((len(terms), len(analyses)), np.nan)
    significant = np.zeros_like(values, dtype=bool)
    for i, term in enumerate(terms):
        for j, analysis in enumerate(analyses):
            x = enrichment[(enrichment.analysis == analysis) & (enrichment.term.str.lower() == term.lower())]
            if len(x):
                values[i, j] = x.signed_enrichment.iloc[0]
                significant[i, j] = x.adjusted_p_value.iloc[0] <= 0.05
    return values, significant


def figure4_factors() -> None:
    stability = pd.read_csv(RESULTS / "factor_stability_bootstrap.tsv", sep="\t")
    enrichment = pd.read_csv(RESULTS / "pathway_rank_enrichment.tsv", sep="\t")
    terms = [
        "heme Metabolism", "TNF-alpha Signaling via NF-kB", "Interferon Gamma Response",
        "Interferon Alpha Response", "MYC Targets V1", "Myogenesis", "Hypoxia",
        "Inflammatory Response", "Oxidative Phosphorylation", "Adipogenesis",
    ]
    analyses = [f"factor_{i}|MSigDB_Hallmark_2020" for i in range(1, 6)]
    values, significant = hallmark_matrix(enrichment, analyses, terms)

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.0), gridspec_kw={"height_ratios": [0.85, 1.35]})
    ax = axes[0]
    data = [stability.loc[stability.factor.eq(i), "cosine_similarity"] for i in range(1, 6)]
    bp = ax.boxplot(data, patch_artist=True, widths=0.6, showfliers=False, medianprops={"color": "white", "lw": 1.4})
    for patch in bp["boxes"]:
        patch.set_facecolor(BLUE)
        patch.set_alpha(0.85)
    ax.set_xticks(range(1, 6), [f"F{i}" for i in range(1, 6)])
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("Aligned loading cosine similarity")
    ax.set_xlabel("Factor")
    ax.set_title("Source-bootstrap stability")
    ax.text(-0.08, 1.06, "a", transform=ax.transAxes, fontsize=12, weight="bold")

    ax = axes[1]
    image = ax.imshow(values, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
    ax.set_xticks(range(5), [f"F{i}" for i in range(1, 6)])
    display_terms = {
        "heme Metabolism": "Heme metabolism",
        "TNF-alpha Signaling via NF-kB": "TNF-alpha signalling via NF-κB",
        "Interferon Gamma Response": "Interferon-gamma response",
        "Interferon Alpha Response": "Interferon-alpha response",
        "MYC Targets V1": "MYC targets V1",
        "Myogenesis": "Myogenesis",
        "Hypoxia": "Hypoxia",
        "Inflammatory Response": "Inflammatory response",
        "Oxidative Phosphorylation": "Oxidative phosphorylation",
        "Adipogenesis": "Adipogenesis",
    }
    clean_terms = [display_terms[term] for term in terms]
    ax.set_yticks(range(len(terms)), clean_terms)
    for i in range(len(terms)):
        for j in range(5):
            if significant[i, j]:
                ax.scatter(j, i, s=18, facecolor="white", edgecolor="#111111", linewidth=0.7, zorder=3)
    ax.set_title("Hallmark enrichment of factor loadings")
    ax.text(-0.08, 1.04, "b", transform=ax.transAxes, fontsize=12, weight="bold")
    cbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("Signed rank enrichment")
    ax.text(0.0, -0.13, "Outlined dot: FDR ≤ 0.05", transform=ax.transAxes, fontsize=7.5)
    fig.subplots_adjust(left=0.35, right=0.94, bottom=0.10, top=0.95, hspace=0.52)
    save(fig, "Figure_4_factor_stability_and_pathways")


def figure5_residuals() -> None:
    enrichment = pd.read_csv(RESULTS / "pathway_rank_enrichment.tsv", sep="\t")
    residual = enrichment[
        enrichment.analysis.str.startswith("residual:")
        & enrichment.analysis.str.endswith("|MSigDB_Hallmark_2020")
    ].copy()
    residual["disease_class"] = residual.analysis.str.extract(r"^residual:(.*)\|MSigDB")
    significant = residual[residual.adjusted_p_value <= 0.05]
    top_terms = (
        significant.assign(abs_enrichment=significant.signed_enrichment.abs())
        .groupby("term", as_index=False).abs_enrichment.max()
        .nlargest(12, "abs_enrichment").term.tolist()
    )
    classes = sorted(residual.disease_class.unique())
    values = np.full((len(top_terms), len(classes)), np.nan)
    sig = np.zeros_like(values, dtype=bool)
    for i, term in enumerate(top_terms):
        for j, disease in enumerate(classes):
            x = residual[(residual.term == term) & (residual.disease_class == disease)]
            if len(x):
                values[i, j] = x.signed_enrichment.iloc[0]
                sig[i, j] = x.adjusted_p_value.iloc[0] <= 0.05

    short = {
        "Cardiometabolic and vascular disease": "Cardiometabolic/vascular",
        "Chronic respiratory disease": "Chronic respiratory",
        "Endocrine autoimmunity": "Endocrine autoimmunity",
        "Genetic and congenital disease": "Genetic/congenital",
        "Hematologic neoplasm": "Hematologic neoplasm",
        "Infectious disease": "Infectious",
        "Inflammatory bowel disease": "Inflammatory bowel",
        "Neurologic disease": "Neurologic",
        "Primary immunodeficiency": "Primary immunodeficiency",
        "Solid cancer": "Solid cancer",
        "Systemic autoimmune and rheumatic disease": "Systemic autoimmune/rheumatic",
        "Transplantation and acute tissue injury": "Transplantation/acute injury",
    }
    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    image = ax.imshow(values, cmap="RdBu_r", vmin=-0.75, vmax=0.75, aspect="auto")
    labels = [short[x] for x in classes]
    ax.set_xticks(range(len(classes)), labels, rotation=50, ha="right", rotation_mode="anchor")
    pathway_labels = {
        "Interferon Alpha Response": "Interferon-alpha response",
        "Heme Metabolism": "Heme metabolism",
        "Interferon Gamma Response": "Interferon-gamma response",
        "Tgf-Beta Signaling": "TGF-beta signalling",
        "Wnt-Beta Catenin Signaling": "WNT/beta-catenin signalling",
        "Protein Secretion": "Protein secretion",
        "Uv Response Dn": "UV response down",
        "Myc Targets V2": "MYC targets V2",
        "Coagulation": "Coagulation",
        "Oxidative Phosphorylation": "Oxidative phosphorylation",
        "Il-6/Jak/Stat3 Signaling": "IL-6/JAK/STAT3 signalling",
        "E2F Targets": "E2F targets",
    }
    formatted_terms = [pathway_labels.get(x.title(), x.replace("Signaling", "signalling")) for x in top_terms]
    ax.set_yticks(range(len(top_terms)), formatted_terms)
    for i in range(len(top_terms)):
        for j in range(len(classes)):
            if sig[i, j]:
                ax.scatter(j, i, s=17, facecolor="white", edgecolor="#111111", linewidth=0.7, zorder=3)
    cbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Signed rank enrichment")
    ax.set_title("Residual pathway signals after shared-structure reconstruction\nOutlined dots indicate FDR ≤ 0.05 within each disease domain")
    fig.subplots_adjust(left=0.30, right=0.95, bottom=0.32, top=0.89)
    save(fig, "Figure_5_residual_pathways")


def main() -> None:
    figure1_design()
    figure2_performance()
    figure3_robustness()
    figure4_factors()
    figure5_residuals()
    print(f"Wrote figures to {FIGURES}")


if __name__ == "__main__":
    main()
