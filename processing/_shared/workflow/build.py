"""Prepare study expression matrices for the shared effect estimator."""

from pathlib import Path
import importlib
import json, sys, gzip, re, hashlib, socket
import os

socket.setdefaulttimeout(30)
import numpy as np
import pandas as pd
from adapters import gene_mapping as repair
from adapters import study_loaders as source_loaders

ROOT = Path(os.environ["BIOLOGY_WORK_ROOT"])
HERE = Path(__file__).parent
OUT = ROOT / "outputs/stage1_three_effects_20260903_v1"
OUT.mkdir(parents=True, exist_ok=True)

AUD = ROOT / "audit"


EVID = {
    x["directory"]: x
    for x in json.loads((AUD / "source_index.json").read_text(encoding="utf8"))
}


def read(p, **kw):
    return pd.read_csv(p, sep="\t", **kw)


def sha(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda: f.read(2**20), b""):
            h.update(b)
    return h.hexdigest()


def ensure_platform(platform, directory=None):
    if not platform or source_loaders.PROBE_MAPS.get(platform):
        return
    annout = OUT / "annotations"
    annout.mkdir(exist_ok=True)
    candidates = list((ROOT / "dataset").glob("*/raw/" + platform + "*")) + list(
        (ROOT / "tmp/eu58_full_model_20260827/annotations").glob(platform + "*")
    )
    if directory:
        candidates += list(
            (ROOT / "dataset" / directory / "raw").glob("*family.soft.gz")
        )
    table = None
    for p in candidates:
        try:
            table = source_loaders.read_annotation(p)
            break
        except Exception:
            pass
    if table is None:
        import urllib.request

        p = annout / (platform + "_platform_only.txt")
        if not p.exists():
            url = (
                "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
                + platform
                + "&targ=self&form=text&view=full"
            )
            with urllib.request.urlopen(url, timeout=30) as response:
                data = response.read(100 * 1024 * 1024 + 1)
            if len(data) > 100 * 1024 * 1024:
                raise ValueError("Platform-only response exceeds size limit")
            p.write_bytes(data)
        table = source_loaders.read_annotation(p)
    cols = {c.lower().replace("_", " ").strip(): c for c in table.columns}
    ic = cols.get("id", table.columns[0])
    gc = next(
        (
            cols[k]
            for k in [
                "entrez gene id",
                "entrez gene",
                "gene id",
                "entrezgeneid",
                "entrez id",
                "locuslink id",
            ]
            if k in cols
        ),
        None,
    )
    sc = next(
        (
            cols[k]
            for k in [
                "gene symbol",
                "symbol",
                "gene name",
                "genename",
                "gene assignment",
            ]
            if k in cols
        ),
        None,
    )
    out = {}
    valid = set(source_loaders.ORG["valid_entrez"])
    accession_genes = {}
    if "gb list" in cols or "gb acc" in cols:
        import sqlite3

        with sqlite3.connect(source_loaders.ORG_DB) as con:
            for accession, gene in con.execute(
                "select a.accession,g.gene_id from accessions a join genes g on a._id=g._id"
            ):
                accession_genes.setdefault(str(accession).split(".")[0], set()).add(
                    str(gene)
                )
        if platform == "GPL4302":
            with sqlite3.connect(ROOT / "external/hgu133plus2.sqlite") as con:
                for accession, gene in con.execute(
                    "select a.accession,p.gene_id from accessions a join probes p on a.probe_id=p.probe_id where p.gene_id is not null"
                ):
                    if gene in valid:
                        accession_genes.setdefault(
                            str(accession).split(".")[0], set()
                        ).add(str(gene))
    sequence_genes = {}
    if platform == "GPL32193":
        # Custom Agilent layout has no gene annotation. Transfer only exact probe sequences
        # with a unique Entrez assignment from two deposited human Agilent platforms.
        for pl, ap in [
            (
                "GPL6480",
                ROOT / "tmp/eu58_full_model_20260827/annotations/GPL6480.annot.gz",
            ),
            ("GPL21185", ROOT / "dataset/gse124272_1/raw/GPL21185_platform.txt"),
        ]:
            if pl not in source_loaders.PROBE_MAPS:
                mp = read(annout / (pl + "_mapping.tsv"), dtype=str)
                source_loaders.PROBE_MAPS[pl] = dict(zip(mp.probe_id, mp.gene_id))
            at = source_loaders.read_annotation(ap)
            sq = next(
                c for c in at.columns if c.upper() in ["SEQUENCE", "PLATFORM_SEQUENCE"]
            )
            for record in at.to_dict("records"):
                gene = source_loaders.PROBE_MAPS[pl].get(str(record["ID"]))
                seq = str(record[sq]).strip().upper()
                if gene and len(seq) >= 40 and re.fullmatch("[ACGT]+", seq):
                    sequence_genes.setdefault(seq, set()).add(gene)
    for rr in table.to_dict("records"):
        gs = set(re.findall(r"\d+", str(rr[gc]))) if gc else set()
        gs &= valid
        if not gs and "gene assignment" in cols:
            for part in str(rr[cols["gene assignment"]]).split(" /// "):
                token = part.split(" // ")[-1].strip()
                if token in valid:
                    gs.add(token)
        if not gs and sc:
            syms = re.split(r"\s*///\s*|\s*//\s*|[;,]", str(rr[sc]))
            gs = {
                source_loaders.ORG["symbol"].get(s.strip().upper(), "") for s in syms
            } - {""}
        if not gs and accession_genes:
            for key in ["gb list", "gb acc"]:
                if key in cols:
                    for acc in re.split(r"[;,\s/]+", str(rr[cols[key]])):
                        gs.update(accession_genes.get(acc.split(".")[0], set()))
        if not gs and sequence_genes:
            gs = sequence_genes.get(
                str(rr.get(cols.get("sequence", ""), "")).strip().upper(), set()
            )
        if len(gs) == 1:
            out[str(rr[ic])] = next(iter(gs))
    source_loaders.PROBE_MAPS[platform] = out
    pd.DataFrame(list(out.items()), columns=["probe_id", "gene_id"]).to_csv(
        annout / (platform + "_mapping.tsv"), sep="\t", index=False
    )


def native(row):
    d = ROOT / "dataset" / row["directory"]
    adapter = HERE / "adapters/studies" / f"{d.name}.py"
    if not adapter.is_file():
        raise ValueError(f"No direct study adapter for {d.name}")
    study = importlib.import_module(f"adapters.studies.{d.name}")
    local = study.extract()
    code = adapter.read_text(encoding="utf8")
    a = local["__case_frame"]
    b = local["__control_frame"]
    if not isinstance(a, pd.DataFrame) or not isinstance(b, pd.DataFrame):
        raise ValueError("captured selection is not a DataFrame")
    x = pd.concat([a, b], axis=1)
    kind = "normalized"
    src = "native log2/normalized expression"
    raw = local.get("__raw_counts")
    if isinstance(raw, pd.DataFrame) and set(x.columns) <= set(raw.columns):
        raw = raw[x.columns].copy()
        if isinstance(raw.index, pd.RangeIndex):
            gene = local.get("gene")
            if gene is None or len(gene) != len(raw):
                raise ValueError("raw counts lost feature IDs")
            raw.index = np.asarray(gene)
        x = raw
        kind = "counts"
        src = "native raw counts intercepted before normalization"
    if d.name == "gse34404_1":
        src = "non_normalized Illumina signal; native sample selection"
    platform = {
        "gse17048_1": "GPL6947",
        "gse48060_1": "GPL570",
        "gse137340_1": "GPL10558",
        "gse51404_1": "GPL6947",
        "gse51405_2": "GPL6947",
        "gse54248_1": "GPL10558",
        "gse94648_1": "GPL19109",
    }.get(row["directory"])
    if not platform:
        alltext = (EVID[row["directory"]].get("platform") or "") + " " + code
        plats = re.findall(r"GPL\d+", alltext)
        platform = plats[0] if plats else None
    return dict(
        x=x,
        cases=list(a.columns),
        controls=list(b.columns),
        kind=kind,
        platform=platform,
        feature_type="auto",
        source=src,
        adapter_id=adapter.stem,
    )


capture_selection = source_loaders.capture_selection


def special(row, sample_selection=None):
    d = row["directory"]
    base = ROOT / "dataset" / d
    acc = row["accession"]
    if d in ["gse44314_1", "gse94648_1"]:
        from extra_loaders import extra

        return extra(row, sys.modules[__name__])
    if row["contrast_id"] == "GSE69683_severe_asthma_vs_HC":
        p = ROOT / "dataset" / d / "raw/GSE69683_series_matrix.txt.gz"
        f, _ = source_loaders.read_geo_series_matrix(p)
        selection = read(sample_selection, dtype=str, keep_default_na=False)
        cases = selection.loc[selection.group.eq("case"), "sample_id"].tolist()
        controls = selection.loc[selection.group.eq("control"), "sample_id"].tolist()
        return capture_selection(
            "",
            f,
            cases,
            controls,
            "feature_id",
            "probe",
            "normalized",
            str(p),
            "GPL13158",
        )
    if d.startswith("gse112057_"):
        m = read(base / "process/sample_metadata.tsv")
        target = {
            1: "Crohn's Disease",
            2: "Ulcerative Colitis",
            3: "Systemic JIA",
            4: "Oligoarticular JIA",
            5: "Polyarticular JIA",
        }[int(d[-1])]
        m["col"] = m.title.str.split("_", n=1).str[0]
        p = base / "raw/GSE112057_RawCounts_dataset.txt.gz"
        f = read(p)
        f.iloc[:, 0] = f.iloc[:, 0].astype(str).str.replace(r"_\d+$", "", regex=True)
        return capture_selection(
            "",
            f,
            m.loc[m.diagnosis.eq(target), "col"].tolist(),
            m.loc[m.diagnosis.eq("Control"), "col"].tolist(),
            f.columns[0],
            "symbol",
            "counts",
            str(p),
        )
    if d.startswith("gse186507_"):
        base = ROOT / "dataset/gse186507_1"
        m = read(base / "process/sample_manifest.tsv")
        p = base / "raw/GSE186507_MSCCR_Blood_counts.txt.gz"
        with gzip.open(p, "rt") as f:
            names = [t.strip('"') for t in f.readline().split()]
        cases = m.loc[
            m.group.eq("CD" if d.endswith("_1") else "UC"), "subject"
        ].tolist()
        controls = m.loc[m.group.eq("Healthy_Control"), "subject"].tolist()
        f = pd.read_csv(
            p,
            sep=r"\s+",
            header=None,
            skiprows=1,
            names=["gene_id"] + names,
            usecols=["gene_id"] + cases + controls,
        )
        return capture_selection(
            "", f, cases, controls, "gene_id", "ensembl", "counts", str(p)
        )
    if d == "gse181228_1":
        m = read(base / "process/sample_manifest.tsv")
        m = m[m.inclusion.eq("include")]
        cases = m.loc[
            m.group.eq("case")
            & m.characteristics.str.contains("visit: Baseline", regex=False),
            "gsm",
        ].tolist()
        ctrl = m.loc[m.group.eq("healthy_control"), "gsm"].tolist()
        p = base / "raw/counts.tsv.gz"
        f = read(p)
        return capture_selection(
            "", f, cases, ctrl, f.columns[0], "entrez", "counts", str(p)
        )
    return None


def validate_sample_selection(row, cases, controls, manifest):
    """Require the source reader to reproduce the declared sample selection."""
    expected = read(manifest, dtype=str, keep_default_na=False)
    columns = ["sample_id", "group", "accession", "study_group"]
    if list(expected.columns) != columns:
        raise ValueError(f"Invalid sample selection columns: {manifest}")
    if expected.sample_id.duplicated().any():
        raise ValueError(f"Duplicate sample IDs in {manifest}")
    actual = pd.DataFrame(
        {
            "sample_id": cases + controls,
            "group": ["case"] * len(cases) + ["control"] * len(controls),
            "accession": row["accession"],
            "study_group": row["group"],
        }
    )
    if not actual.equals(expected):
        raise ValueError(f"Source sample IDs, groups, or order differ from {manifest}")


def prepare(row, frozen_mapping=None, sample_selection=None):
    """Build one expression matrix using the published feature mapping when supplied."""
    d = row["directory"]
    target = OUT / "contrasts" / row["contrast_id"]
    target.mkdir(parents=True, exist_ok=True)
    if (target / "input.json").exists():
        if sample_selection is not None:
            samples = read(target / "samples.tsv", dtype=str, keep_default_na=False)
            validate_sample_selection(
                row,
                samples.loc[samples.group.eq("case"), "sample_id"].tolist(),
                samples.loc[samples.group.eq("control"), "sample_id"].tolist(),
                sample_selection,
            )
        info = json.loads((target / "input.json").read_text(encoding="utf8"))
        info.update(row)
        (target / "input.json").write_text(
            json.dumps(info, ensure_ascii=False, indent=2), encoding="utf8"
        )
        return info
    r = special(row, sample_selection)
    if r is None:
        try:
            r = native(row)
        except Exception:
            r = (
                source_loaders.load_declared_source(row, sample_selection)
                if sample_selection is not None
                else None
            )
            if r is None:
                from extra_loaders import extra

                r = extra(row, sys.modules[__name__])
                if r is None:
                    raise
    x = r["x"].apply(pd.to_numeric, errors="coerce")
    cases = list(map(str, r["cases"]))
    ctrl = list(map(str, r["controls"]))
    x.columns = x.columns.astype(str)
    assert (
        len(cases) >= 2
        and len(ctrl) >= 2
        and not set(cases) & set(ctrl)
        and len(set(cases + ctrl)) == len(cases + ctrl)
    ), "invalid/overlapping samples"
    if sample_selection is not None:
        validate_sample_selection(row, cases, ctrl, sample_selection)
    for k, n in [("n_case_audit", len(cases)), ("n_control_audit", len(ctrl))]:
        if row.get(k) is not None:
            assert int(row[k]) == n, f"{k}: {n} != {row[k]}"
    ids = x.index.astype(str)
    if frozen_mapping is not None:
        # Lock feature order and gene IDs to the original analysis.
        frozen = read(frozen_mapping, dtype=str, keep_default_na=False)
        assert list(ids) == frozen.original_id.tolist(), (
            "Source feature IDs/order differ from frozen annotation mapping"
        )
        mapped = frozen.gene_id.to_numpy()
    else:
        typ = r["feature_type"]
        if typ == "auto":
            mapped = repair.auto_map_features(ids, source_loaders)[0]
            if (mapped != "").mean() < 0.5 and r.get("platform"):
                ensure_platform(r["platform"], d)
                mapped = np.array(
                    [
                        source_loaders.map_feature(v, "probe", r["platform"]) or ""
                        for v in ids
                    ]
                )
        else:
            if typ == "probe":
                ensure_platform(r.get("platform"), d)
            mapped = np.array(
                [
                    source_loaders.map_feature(v, typ, r.get("platform")) or ""
                    for v in ids
                ]
            )
    mapt = pd.DataFrame({"original_id": ids, "gene_id": mapped})
    mapt.to_csv(target / "feature_mapping.tsv.gz", sep="\t", index=False)
    x.index = mapped
    x = x[x.index != ""]
    pre = len(x)
    # All feature collapse precedes effects. Counts sum; normalized probes average.
    x = (
        x.groupby(level=0).sum(min_count=1)
        if r["kind"] == "counts"
        else x.groupby(level=0).mean()
    )
    x.index.name = "gene_id"
    x = x[cases + ctrl]
    missing = int(x.isna().sum().sum())
    x = x.replace([np.inf, -np.inf], np.nan).dropna()
    if r["kind"] == "counts":
        assert (x.to_numpy() >= 0).all(), "negative counts"
    assert len(x) >= 1000, f"Only {len(x)} mapped complete genes"
    x.to_csv(target / "input_expression.tsv.gz", sep="\t", float_format="%.10g")
    sm = pd.DataFrame(
        {
            "sample_id": cases + ctrl,
            "group": ["case"] * len(cases) + ["control"] * len(ctrl),
            "accession": row["accession"],
            "study_group": row["group"],
        }
    )
    sm.to_csv(target / "samples.tsv", sep="\t", index=False)
    info = {
        **row,
        "n_case": len(cases),
        "n_control": len(ctrl),
        "input_kind": r["kind"],
        "source": r["source"],
        "platform": r.get("platform"),
        "input_genes": len(x),
        "mapped_rows_before_collapse": pre,
        "missing_values_removed_with_genes": missing,
        "input_sha256": sha(target / "input_expression.tsv.gz"),
        "source_adapter": r.get("adapter_id"),
        "scope": "core" if row["review_status"] in "AB" else "supplemental",
    }
    (target / "input.json").write_text(
        json.dumps(info, ensure_ascii=False, indent=2), encoding="utf8"
    )
    return info
