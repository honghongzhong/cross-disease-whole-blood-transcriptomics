from pathlib import Path
from collections import defaultdict
import gzip
import numpy as np
import pandas as pd
from scipy.stats import t

HERE = Path(__file__).resolve().parent
RAW = HERE.parent / "raw"
OUT = HERE.parent / "result"


def main():
    meta = pd.read_csv(HERE / "sample_manifest.tsv", sep="\t", dtype=str)
    groups = dict(zip(meta["sample_id"], meta["primary_role"]))
    case = {s for s, g in groups.items() if g == "case"}
    control = {s for s, g in groups.items() if g == "control"}
    assert len(case) == 304 and len(control) == 96

    counts_path = RAW / "GSE190125_Counts_for_GEO.txt.gz"
    library_sizes = defaultdict(float)
    with gzip.open(counts_path, "rt", encoding="utf-8-sig", errors="replace") as f:
        next(f)
        for line in f:
            sample, _ens, _gene, _chr, _type, raw = line.rstrip("\n").split("\t")
            library_sizes[sample] += float(raw)
    assert set(library_sizes) == set(groups)

    # Streaming summary avoids materializing the approximately 24-million-row long file.
    stats = {}
    with gzip.open(counts_path, "rt", encoding="utf-8-sig", errors="replace") as f:
        next(f)
        for line in f:
            sample, _ens, gene, _chr, _type, raw = line.rstrip("\n").split("\t")
            if not gene or sample not in groups:
                continue
            value = np.log2(float(raw) / library_sizes[sample] * 1e6 + 1.0)
            if gene not in stats:
                stats[gene] = [0.0, 0.0, 0, 0.0, 0.0, 0]
            s = stats[gene]
            if sample in case:
                s[0] += value; s[1] += value * value; s[2] += 1
            else:
                s[3] += value; s[4] += value * value; s[5] += 1

    rows = []
    for gene, s in stats.items():
        if s[2] < 2 or s[5] < 2:
            continue
        ma, mb = s[0] / s[2], s[3] / s[5]
        va = max((s[1] - s[0] * s[0] / s[2]) / (s[2] - 1), 0.0)
        vb = max((s[4] - s[3] * s[3] / s[5]) / (s[5] - 1), 0.0)
        u2 = va / s[2] + vb / s[5]
        if not np.isfinite(u2) or u2 <= 0:
            continue
        e = ma - mb
        df = u2 * u2 / ((va / s[2]) ** 2 / (s[2] - 1) + (vb / s[5]) ** 2 / (s[5] - 1))
        p = 2.0 * t.sf(abs(e / np.sqrt(u2)), df)
        rows.append((gene, e, np.sqrt(u2), p))
    out = pd.DataFrame(rows, columns=["gene_id", "E", "U", "p_value"]).sort_values("gene_id")
    p = out["p_value"].to_numpy()
    order = np.argsort(p); q = np.empty_like(p)
    q[order] = np.clip(np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1], 0, 1)
    out["q_value"] = q
    OUT.mkdir(exist_ok=True)
    target = OUT / "GSE190125_trisomy21_vs_euploid_control_E_U.tsv"
    out.to_csv(target, sep="\t", index=False, float_format="%.10g")
    print(f"Wrote {len(out)} genes to {target}")


if __name__ == "__main__":
    main()
