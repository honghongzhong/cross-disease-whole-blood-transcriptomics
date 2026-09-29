"""Study-specific input routes used when the standard reader cannot load a GEO deposit."""

import re, gzip, json
import numpy as np
import pandas as pd


def extra(row, api):
    root = api.ROOT
    d = row["directory"]
    acc = row["accession"]
    base = root / "dataset" / d
    raw = base / "raw"
    cid = row["contrast_id"]
    read = api.read
    freeze = (
        root
        / "census_step0/candidate_extensions_2026-08-22_v1"
        / (acc + "_lock_2026-08-22_v1")
    )
    cand = root / "phase1b_EU_coordination_v1_2026-08-23/01_candidates" / acc

    def finish(f, ca, co, kind="normalized", typ="probe", platform=None, source=""):
        return api.capture_selection(
            cid, f, list(ca), list(co), f.columns[0], typ, kind, str(source), platform
        )

    def series(p):
        f, meta = api.source_loaders.read_geo_series_matrix(p)
        plats = meta.get("!Sample_platform_id", [[]])[0]
        return f, meta, plats[0] if plats else None

    if d in [
        "gse201332_1",
        "gse44314_1",
        "gse69063_1",
        "gse94648_1",
    ]:
        m = read(base / "process/sample_manifest.tsv")
        p = next(raw.glob("*series_matrix*"))
        f, meta, pl = series(p)
        if d == "gse44314_1":
            ca = m.loc[
                m.title.str.contains(
                    "Fulminant" if "fulminant" in cid else "Classical 1A", case=False
                ),
                "sample_id",
            ]
            co = m.loc[m.role.eq("control"), "sample_id"]
        elif d == "gse94648_1":
            m = m[~m.title.str.contains("bis", case=False)]
            ca = m.loc[m.title.str.contains("CD" if "_CD_" in cid else "UC"), "gsm"]
            co = m.loc[m.group.eq("healthy_control"), "gsm"]
        elif d == "gse69063_1":
            ca = m.loc[m.role.eq("case") & m.included.eq(1), "sample_id"]
            co = m.loc[m.role.eq("control") & m.included.eq(1), "sample_id"]
        else:
            ca = m.loc[m.group.eq("case"), "gsm"]
            co = m.loc[m.group.eq("control"), "gsm"]
        return finish(f, ca, co, platform=pl, source=p)
    if d in ["gse136371_1", "gse185855_1"]:
        if d == "gse136371_1":
            m = read(freeze / (acc + "_sample_metadata.tsv"))
            ca = m.loc[m.diagnosis.ne("healthy control"), "gsm"]
            co = m.loc[m.diagnosis.eq("healthy control"), "gsm"]
            p = raw / (acc + "_counts_matrix.tsv.gz")
            typ = "symbol"
        elif d == "gse185855_1":
            m = read(freeze / (acc + "_baseline_sample_lock.tsv"))
            ca = m.loc[m.group.eq("MDD"), "gsm"]
            co = m.loc[m.group.ne("MDD"), "gsm"]
            p = raw / (acc + "_baseline_counts_matrix.tsv.gz")
            typ = "ensembl"
        f = pd.read_csv(p, sep="\t")
        f.columns = f.columns.str.strip()
        if d == "gse136371_1":
            f.iloc[:, 0] = f.iloc[:, 0].str.replace(r"_\d+$", "", regex=True)
        return finish(f, ca, co, "counts", typ, source=p)
    if d == "gse18781_3":
        m = read(freeze / (acc + "_sample_subject_data_lock_v1.tsv"))
        p = next(raw.glob("*series_matrix*"))
        f, meta, pl = series(p)
        m = m[m.contrast_id.eq(cid)]
        ca = m.loc[m.disease_group.eq("case"), "sample_id"]
        co = m.loc[m.disease_group.eq("healthy_control"), "sample_id"]
        return finish(f, ca, co, platform=pl, source=p)
    if d == "gse177044_1":
        m = read(cand / "sample_manifest.tsv")
        fs = [
            pd.read_csv(p)
            .set_index("Geneid")
            .drop(columns=["gene_name"], errors="ignore")
            for p in sorted(raw.glob("*.csv.gz"))
        ]
        f = pd.concat(fs, axis=1)
        assert not f.columns.duplicated().any()
        return finish(
            f.reset_index(),
            m.loc[m.disease.eq("UC"), "title"],
            m.loc[m.disease.eq("Control"), "title"],
            "counts",
            "ensembl",
            source=";".join(str(p) for p in raw.glob("*.csv.gz")),
        )
    if d == "gse190125_1":
        p = raw / "GSE190125_Counts_for_GEO.txt.gz"
        m = read(base / "process/sample_manifest.tsv")
        arrays = {}
        for chunk in pd.read_csv(p, sep="\t", chunksize=200000):
            for sid, part in chunk.groupby(chunk.columns[0], sort=False):
                ser = pd.Series(
                    pd.to_numeric(part.iloc[:, -1]).values,
                    index=part.iloc[:, 1].astype(str),
                )
                arrays.setdefault(sid, []).append(ser)
        f = pd.DataFrame(
            {
                sid: pd.concat(parts).groupby(level=0).sum(min_count=1)
                for sid, parts in arrays.items()
            }
        ).reset_index(names="gene_id")
        return finish(
            f,
            m.loc[m.primary_role.eq("case"), "sample_id"],
            m.loc[m.primary_role.eq("control"), "sample_id"],
            "counts",
            "ensembl",
            source=p,
        )
    if d == "gse28750_1":
        p = api.OUT / "cel_normalized" / d / "rma.tsv.gz"
        f = read(p)
        m = read(base / "process/sample_manifest.tsv")
        m = m[m.include_locked.astype(str).str.upper().eq("TRUE")]
        key = "sample_id"
        return finish(
            f,
            m.loc[m.disease_group.eq("case"), key],
            m.loc[m.disease_group.eq("healthy_control"), key],
            platform="GPL570",
            source="Affymetrix RMA from " + str(raw),
        )
    if d == "gse184876_1":
        p = raw / "counts.txt.gz"
        m = read(base / "process/sample_manifest.tsv")
        f = read(p)
        return finish(
            f,
            m.loc[m.group.eq("case"), "count_column"],
            m.loc[m.group.eq("control"), "count_column"],
            "counts",
            "ensembl",
            source=p,
        )
    return None
