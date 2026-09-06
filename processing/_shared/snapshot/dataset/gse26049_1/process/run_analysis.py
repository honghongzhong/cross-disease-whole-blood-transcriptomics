from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent
RESULT = HERE.parent / "result"
def main():
    files = [p for p in RESULT.iterdir() if p.is_file() and "processed_E_U" not in p.name]
    if not files:
        raise FileNotFoundError("No source result table found")
    src = sorted(files, key=lambda p: p.stat().st_size, reverse=True)[0]
    d = pd.read_csv(src, sep="	", compression="infer")
    def pick(names):
        for name in names:
            if name in d.columns:
                return pd.to_numeric(d[name], errors="coerce")
        return pd.Series(np.nan, index=d.index)
    gid = next((d[name].astype(str) for name in ["gene_id", "entrez_id", "probe_id", "canonical_symbol"] if name in d.columns), pd.Series(d.index.astype(str)))
    out = pd.DataFrame({"gene_id": gid, "E": pick(["E", "effect", "effect_raw"]), "U": pick(["U", "uncertainty", "se", "se_raw"]), "p_value": pick(["p_value", "pvalue"]), "q_value": pick(["q_value", "adjusted_p_value", "padj"])})
    out = out[np.isfinite(out.E) & np.isfinite(out.U) & (out.U > 0)].drop_duplicates("gene_id").sort_values("gene_id")
    out.to_csv(RESULT / "processed_E_U.tsv", sep="	", index=False, float_format="%.10g")
    print(f"Wrote {len(out)} rows to {RESULT / 'processed_E_U.tsv'}")
if __name__ == "__main__":
    main()
