from pathlib import Path
import csv
import gzip
import re

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind_from_stats


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RAW = ROOT / "raw"
RESULT = ROOT / "result"

CONTROL_IDS = {"927", "928", "938", "941", "942", "943", "944", "947", "967", "970"}


def bh_fdr(p_values):
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(ranked)
    q[order] = np.minimum(ranked, 1.0)
    return q


def read_sample_metadata(path):
    meta = {}
    with gzip.open(path, "rt", errors="replace", newline="") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                if line.startswith("!series_matrix_table_begin"):
                    break
                continue
            fields = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            meta[fields[0][8:]] = fields[1:]
    return meta


def subject_id(title):
    m = re.search(r"patient ([^ ]+)", title, flags=re.I)
    if m:
        return m.group(1)
    m = re.search(r"healthy control ([^ ]+)$", title, flags=re.I)
    return m.group(1) if m else ""


def main():
    RESULT.mkdir(parents=True, exist_ok=True)
    matrix_path = RAW / "GSE69063_series_matrix.txt.gz"
    meta = read_sample_metadata(matrix_path)
    titles = meta["title"]
    accessions = meta["geo_accession"]

    case_idx = [
        i for i, title in enumerate(titles)
        if "anaphylaxis patient" in title.lower() and "time t0" in title.lower()
    ]
    control_idx = [
        i for i, title in enumerate(titles)
        if re.search(r"healthy control (\d+)$", title, flags=re.I)
        and re.search(r"healthy control (\d+)$", title, flags=re.I).group(1) in CONTROL_IDS
    ]
    if len(case_idx) != 16 or len(control_idx) != 10:
        raise RuntimeError(f"Unexpected locked sample counts: case={len(case_idx)}, control={len(control_idx)}")

    manifest = []
    for i, (gsm, title) in enumerate(zip(accessions, titles)):
        is_case = i in case_idx
        is_control = i in control_idx
        included = is_case or is_control
        if is_case:
            role, cohort, tp, reason = "case", "anaphylaxis", "T0", ""
        elif is_control:
            role, cohort, tp, reason = "control", "anaphylaxis_matched_healthy", "baseline", ""
        else:
            role, cohort, tp = "excluded", "other_study_phenotype", ""
            reason = "Different phenotype, non-baseline timepoint, or healthy control from another cohort"
        manifest.append({
            "sample_id": gsm,
            "title": title,
            "role": role,
            "cohort": cohort,
            "timepoint": tp,
            "subject_id": subject_id(title),
            "included": int(included),
            "exclusion_reason": reason,
        })
    pd.DataFrame(manifest).to_csv(HERE / "sample_manifest.tsv", sep="\t", index=False)

    matrix = pd.read_csv(matrix_path, sep="\t", comment="!", compression="gzip")
    matrix = matrix.rename(columns={matrix.columns[0]: "probe_id"})
    case_cols = [accessions[i] for i in case_idx]
    ctrl_cols = [accessions[i] for i in control_idx]
    case = matrix[case_cols].apply(pd.to_numeric, errors="coerce")
    ctrl = matrix[ctrl_cols].apply(pd.to_numeric, errors="coerce")
    n_case = case.notna().sum(axis=1)
    n_ctrl = ctrl.notna().sum(axis=1)
    mean_case = case.mean(axis=1)
    mean_ctrl = ctrl.mean(axis=1)
    var_case = case.var(axis=1, ddof=1)
    var_ctrl = ctrl.var(axis=1, ddof=1)
    u = np.sqrt(var_case / n_case + var_ctrl / n_ctrl)
    e = mean_case - mean_ctrl
    p = ttest_ind_from_stats(
        mean1=mean_case, std1=np.sqrt(var_case), nobs1=n_case,
        mean2=mean_ctrl, std2=np.sqrt(var_ctrl), nobs2=n_ctrl,
        equal_var=False,
    ).pvalue
    q = bh_fdr(np.nan_to_num(p, nan=1.0))
    out = pd.DataFrame({
        "probe_id": matrix["probe_id"].astype(str),
        "E": e,
        "U": u,
        "p_value": p,
        "q_value": q,
        "n_case": n_case,
        "n_control": n_ctrl,
    })
    out = out[np.isfinite(out["E"]) & np.isfinite(out["U"]) & (out["U"] > 0)]
    out = out.drop_duplicates("probe_id").sort_values("probe_id")
    out.to_csv(RESULT / "GSE69063_Anaphylaxis_vs_HC_E_U.tsv", sep="\t", index=False, float_format="%.10g")
    print(f"case={len(case_idx)} control={len(control_idx)} probes={len(out)}")


if __name__ == "__main__":
    main()
