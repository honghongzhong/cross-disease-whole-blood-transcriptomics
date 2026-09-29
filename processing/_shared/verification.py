"""Compare regenerated matrices and effect tables with published references."""

from pathlib import Path

import numpy as np
import pandas as pd

VIEWS = ("hedges_g", "log2fc", "shrunken_log2fc")


def check_preparation(
    target: Path, original: Path, sample_manifest: Path
) -> tuple[bool | None, bool, list[int]]:
    """Check sample identity and prepared expression within text precision."""
    expression = pd.read_csv(target / "input_expression.tsv.gz", sep="\t", index_col=0)
    comparison = None
    reference_file = original / "input_expression.tsv.gz"
    if reference_file.exists():
        reference = pd.read_csv(reference_file, sep="\t", index_col=0)
        comparison = (
            expression.index.astype(str).equals(reference.index.astype(str))
            and expression.columns.equals(reference.columns)
            and expression.shape == reference.shape
            and np.allclose(
                expression, reference, rtol=1e-9, atol=1e-10, equal_nan=True
            )
        )
    samples = pd.read_csv(target / "samples.tsv", sep="\t").equals(
        pd.read_csv(sample_manifest, sep="\t")
    )
    return comparison, samples, list(expression.shape)


def compare_effect_tables(target: Path, source: Path, contrast_id: str) -> list[dict]:
    """Compare all effect and uncertainty columns when a reference is present."""
    comparisons = []
    for view in VIEWS:
        fresh = pd.read_csv(
            target / f"{view}_E_U.tsv.gz", sep="\t", dtype={"gene_id": str}
        )
        expected = (
            source
            / "effect-uncertainty/contrasts"
            / contrast_id
            / f"{contrast_id}__{view}__E_U.tsv.gz"
        )
        if not expected.exists():
            comparisons.append(
                {
                    "view": view,
                    "rows": len(fresh),
                    "match": None,
                    "note": "Generated from source inputs; original comparison table not provided.",
                }
            )
            continue
        reference = pd.read_csv(expected, sep="\t", dtype={"gene_id": str})
        same_ids = fresh.gene_id.equals(reference.gene_id)
        same_shape = fresh.shape == reference.shape
        columns = ["effect", "uncertainty"]
        difference = (
            float(np.max(np.abs(fresh[columns].values - reference[columns].values)))
            if same_shape and same_ids
            else None
        )
        match = (
            same_ids
            and same_shape
            and np.allclose(fresh[columns], reference[columns], rtol=1e-7, atol=1e-9)
            and fresh[["uncertainty_type", "contrast_id"]].equals(
                reference[["uncertainty_type", "contrast_id"]]
            )
        )
        comparisons.append(
            {
                "view": view,
                "rows": len(fresh),
                "match": bool(match),
                "max_absolute_difference": difference,
            }
        )
    return comparisons


def reproduction_state(
    prepared: bool | None, samples: bool, comparisons: list[dict]
) -> str:
    """Distinguish a matching result from one lacking reference data."""
    if prepared is None or any(item["match"] is None for item in comparisons):
        return "GENERATED_REFERENCE_NOT_PROVIDED"
    if prepared and samples and all(item["match"] for item in comparisons):
        return "PASS"
    return "MISMATCH"
