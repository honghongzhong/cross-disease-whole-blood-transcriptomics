"""Reproduce the GSE112057 Oligoarticular JIA E/U result with PyDESeq2.

Run from the project root with: python dataset/<dataset_id>/process/run_analysis.py
"""
from pathlib import Path
import re
import sqlite3
import sys
import argparse

import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


PROCESS = Path(__file__).resolve().parent
DATASET = PROCESS.parent
RAW = DATASET / "raw"
LOCK = PROCESS
RESULT = DATASET / "result"

parser = argparse.ArgumentParser()
parser.add_argument("--target", default="Systemic JIA")
parser.add_argument("--case", type=int, default=26)
parser.add_argument("--control", type=int, default=12)
parser.add_argument("--contrast", default="GSE112057_Systemic_JIA_vs_HC")
parser.add_argument("--output", default="GSE112057_Systemic_JIA_vs_HC_E_U.tsv")
args = parser.parse_args()

meta = pd.read_csv(LOCK / "sample_metadata.tsv", sep="\t", dtype=str)
raw = pd.read_csv(RAW / "GSE112057_RawCounts_dataset.txt.gz", sep="\t")
symbols = raw.iloc[:, 0].astype(str).str.replace(r"_\d+$", "", regex=True)
counts = raw.iloc[:, 1:].copy()
counts.columns = counts.columns.astype(str)
counts = counts.apply(pd.to_numeric, errors="raise").astype("int64")
counts.index = symbols
counts = counts.groupby(level=0, sort=False).sum()

sel = meta[meta["diagnosis"].isin(["Control", args.target])].copy()
sel["group"] = sel["diagnosis"].map({"Control": "Healthy", args.target: "Case"})
assert len(sel) == args.case + args.control and (sel["group"] == "Case").sum() == args.case and (sel["group"] == "Healthy").sum() == args.control
sel["raw_sample"] = sel["title"].str.split("_", n=1).str[0]
assert set(sel["raw_sample"]) <= set(counts.columns)
sel = sel.set_index("raw_sample").loc[:, ["gsm", "group"]]
sel = sel.rename_axis("raw_sample")
counts = counts.loc[:, sel.index].T
counts.index = sel["gsm"].values
sel = sel.set_index("gsm").loc[:, ["group"]]
counts = counts.loc[:, counts.sum(axis=0) > 0]
keep = counts.sum(axis=0) >= 10
counts = counts.loc[:, keep]

dds = DeseqDataSet(counts=counts, metadata=sel, design_factors="group", refit_cooks=True, n_cpus=1)
dds.deseq2()
stats = DeseqStats(dds, contrast=["group", "Case", "Healthy"], independent_filter=False, n_cpus=1)
stats.summary()
res = stats.results_df.reset_index()
res = res.rename(columns={res.columns[0]: "gene_symbol"})

db = Path(__file__).resolve().parents[3] / "tools" / "r461_bioc323_gse153315" / "library" / "org.Hs.eg.db" / "extdata" / "org.Hs.eg.sqlite"
if not db.exists():
    raise FileNotFoundError(f"Missing Entrez annotation database: {db}")
with sqlite3.connect(db) as con:
    mapping = pd.read_sql_query(
        "select i.symbol as gene_symbol, cast(g.gene_id as text) as gene_id "
        "from gene_info i join genes g on g._id=i._id where i.symbol is not null", con
    )
mapping = mapping.drop_duplicates("gene_symbol")
res = res.merge(mapping, on="gene_symbol", how="inner").drop_duplicates("gene_id")
res = res.rename(columns={"log2FoldChange": "effect", "lfcSE": "se", "pvalue": "p_value", "padj": "adjusted_p_value"})
res = res[["gene_id", "effect", "se", "p_value", "adjusted_p_value"]]
res = res[res["effect"].notna() & res["se"].notna() & (res["se"] > 0)].copy()
res["n_case"] = args.case
res["n_control"] = args.control
res["contrast_id"] = args.contrast
res["analysis_method"] = "PyDESeq2 0.5.4; design=~group; raw counts; Entrez mapping from org.Hs.eg.db"
out = RESULT / args.output
res.to_csv(out, sep="\t", index=False)
sel.reset_index().to_csv(PROCESS / "sample_selection.tsv", sep="\t", index=False)
print(f"Wrote {len(res)} genes to {out}")

