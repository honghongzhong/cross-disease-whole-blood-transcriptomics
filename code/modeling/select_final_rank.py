"""Select the deployment/final-model rank using all 57 contrasts and grouped inner CV."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("eu57_core", HERE / "run_experiment.py")
core = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = core
SPEC.loader.exec_module(core)


def main() -> None:
    # Keep helper-script outputs in the same public output tree used by run_all.py.
    out = core.DEFAULT_OUTPUT / "final_rank_selection"
    out.mkdir(parents=True, exist_ok=True)
    manifest, _, effect, uncertainty = core.load_inputs(
        core.DEFAULT_MANIFEST, core.DEFAULT_AUDIT, core.DEFAULT_CLASS_MAP
    )
    args = SimpleNamespace(
        inner_folds=5, inner_masks=2, max_k=6, coverage=0.90, genes=1000,
        mask_fraction=0.20, max_iterations=80, ridge=1e-3,
    )
    chosen, table = core.inner_select_rank(
        effect, uncertainty, manifest, np.ones(len(manifest), dtype=bool), args, 99, out
    )
    (out / "selected_rank.json").write_text(
        json.dumps({"selected_K": chosen, "selection_rule": "one_standard_error"}, indent=2),
        encoding="utf-8",
    )
    print(table.to_string(index=False))
    print(f"selected_K={chosen}")


if __name__ == "__main__":
    main()
