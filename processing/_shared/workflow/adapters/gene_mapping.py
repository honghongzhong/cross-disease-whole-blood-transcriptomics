"""Map deposited feature identifiers to canonical Entrez gene IDs."""

from types import ModuleType
from typing import Iterable

import numpy as np


def auto_map_features(
    values: Iterable[object], loaders: ModuleType
) -> tuple[np.ndarray, np.ndarray]:
    """Resolve Entrez IDs, Ensembl IDs, or symbols using the frozen maps."""
    mapped: list[str] = []
    symbols: list[str] = []
    for value in map(str, values):
        raw = value.strip().replace(".0", "")
        if raw in loaders.ORG["valid_entrez"]:
            gene = raw
        elif raw.upper().startswith("ENSG"):
            gene = loaders.ORG["ensembl"].get(raw.split(".")[0].upper(), "")
        else:
            gene = loaders.ORG["symbol"].get(raw.upper(), "")
        mapped.append(gene)
        symbols.append(loaders.ORG["entrez_symbol"].get(gene, "") if gene else "")
    return np.asarray(mapped, dtype=object), np.asarray(symbols, dtype=object)
