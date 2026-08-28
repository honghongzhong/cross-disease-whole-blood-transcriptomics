"""Compare manuscript-defining outputs from a new run with reference results."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
KEY_FILES = [
    "primary_result.json",
    "outer_prediction_components.tsv",
    "outer_selected_ranks.tsv",
    "outer_class_metrics.tsv",
    "leave_one_accession_out_summary.json",
    "panel_size_sensitivity.tsv",
    "final_model_parameters.npz",
    "factor_stability_bootstrap.tsv",
    "pathway_rank_enrichment.tsv",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", type=Path, required=True)
    args = parser.parse_args()
    candidate = args.candidate.resolve()
    reference = ROOT / "reference_results"
    comparisons = []
    for name in KEY_FILES:
        expected, observed = reference / name, candidate / name
        comparisons.append(
            {
                "file": name,
                "reference_exists": expected.is_file(),
                "candidate_exists": observed.is_file(),
                "sha256_identical": expected.is_file() and observed.is_file() and sha256(expected) == sha256(observed),
            }
        )
    result = {"schema": "eu57_key_result_comparison_v1", "comparisons": comparisons}
    output = candidate / "reference_comparison.json"
    output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    if not all(row["sha256_identical"] for row in comparisons):
        raise SystemExit(1)


if __name__ == "__main__":
    main()

