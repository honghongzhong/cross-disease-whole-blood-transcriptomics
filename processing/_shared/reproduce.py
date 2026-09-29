"""Rebuild from source-level data in an isolated tree, then compare all E/U rows."""

from pathlib import Path
import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time
import traceback

# Per-comparison launchers execute this file with runpy.run_path.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from staging import stage
from workflow_setup import install_workflow
from verification import check_preparation, compare_effect_tables, reproduction_state
import pandas as pd

BASE = Path(__file__).resolve().parent.parent


def run(args):
    source = Path(args.source_root).resolve()
    work = Path(args.work_root).resolve()
    if args.stage_only:
        stage(source, work)
        return
    if not (work / "STAGED.json").exists():
        stage(source, work)
    cid = args.contrast
    pkg = BASE / cid
    metadata = json.loads((pkg / "meta.json").read_text(encoding="utf8"))
    info = metadata["input"]
    audit = {
        "contrast_id": cid,
        "code_copy": "COMPLETE",
        "reproduction": "RUNNING",
        "started": time.time(),
        "limitations": [
            "GEO-deposited sample-level expression is not necessarily instrument-level raw data.",
            "Unadjusted case-minus-control model; uncertainty assumes independent samples.",
            "Hedges g sampling SE is an approximation; ashr uses normal likelihood.",
        ],
    }
    output = work / "outputs/stage1_three_effects_20260903_v1"
    target = output / "contrasts" / cid
    reads = set()
    checking = False

    def hook(event, a):
        if event == "open" and a and isinstance(a[0], (str, bytes)):
            p = str(a[0]).replace("\\", "/")
            if not checking and (
                str(source).replace("\\", "/") + "/effect-uncertainty/" in p
                or str(source).replace("\\", "/")
                + "/outputs/stage1_three_effects_20260903_v1/contrasts/"
                in p
            ):
                raise RuntimeError("Attempted reuse of original computed output: " + p)
            if not checking and (str(work).replace("\\", "/") in p):
                reads.add(p)

    sys.addaudithook(hook)
    try:
        missing_inputs = [
            item["relative"]
            for item in info["source_inputs"]
            if not item["bundled"] and not (work / item["relative"]).is_file()
        ]
        if missing_inputs:
            raise FileNotFoundError(
                "Required source inputs are missing under --source-root: "
                + ", ".join(missing_inputs)
            )
        qold = metadata["qc"]
        historical = (
            "shrinkage_weight_prior" in qold or "noncurrent_entrez_excluded" in qold
        )
        code = install_workflow(cid, work)
        audit["annotation_mapping"] = (
            "Frozen per-contrast feature-to-gene mapping; source IDs/order checked exactly"
        )
        if historical:
            audit["historical_serialized_refit"] = True
        # Respect the original platform mapping, including the sequence-transfer maps.
        sys.path.insert(0, str(code))
        os.environ["BIOLOGY_ANNOTATION_CACHE"] = str(
            BASE / "_shared/annotation_cache.json.gz"
        )
        os.environ["BIOLOGY_WORK_ROOT"] = str(work)
        spec = importlib.util.spec_from_file_location("build", code / "build.py")
        b = importlib.util.module_from_spec(spec)
        sys.modules["build"] = b
        spec.loader.exec_module(b)
        if target.exists():
            if (target / "qc.json").exists():
                raise RuntimeError("Completed output exists; use a fresh --work-root")
            # Preserve a failed attempt before restarting preparation from source data.
            backup = target.with_name(cid + "_failed_" + str(time.time_ns()))
            assert target.resolve().is_relative_to(
                work
            ) and backup.resolve().is_relative_to(work)
            target.rename(backup)
        if info["directory"] in ["gse28750_1", "gse55201_1"]:
            import tarfile

            raw = work / "dataset" / info["directory"] / "raw"
            wanted = set(
                pd.read_csv(pkg / "sample_selection.tsv", sep="\t").sample_id.astype(
                    str
                )
            )
            available = {
                p.name.split("_")[0].split(".")[0]
                for p in raw.iterdir()
                if p.name.lower().endswith((".cel", ".cel.gz"))
            }
            if wanted - available:
                for archive in raw.iterdir():
                    if not archive.name.lower().endswith((".tar", ".tar.gz", ".tgz")):
                        continue
                    with tarfile.open(archive) as tf:
                        for member in tf:
                            name = Path(member.name).name
                            sid = name.split("_")[0].split(".")[0]
                            if (
                                not member.isfile()
                                or sid not in wanted - available
                                or not name.lower().endswith((".cel", ".cel.gz"))
                            ):
                                continue
                            dest = raw / name
                            assert dest.resolve().is_relative_to(work)
                            with tf.extractfile(member) as src, dest.open("wb") as dst:
                                shutil.copyfileobj(src, dst)
                            available.add(sid)
                assert not wanted - available, "Missing selected CEL files: " + str(
                    sorted(wanted - available)
                )
            s = (
                (code / "cel.R")
                .read_text(encoding="utf8")
                .replace("c('gse28750_1','gse55201_1')", repr(info["directory"]))
            )
            s = (
                ".libPaths(c('"
                + (source / "tools/r461_bioc323_gse153315/library").as_posix()
                + "',.libPaths()))\n"
                + s
            )
            (code / "cel_one.R").write_text(s, encoding="utf8")
            subprocess.run([args.rscript, str(code / "cel_one.R")], check=True)
        b.prepare(
            info,
            pkg / "original_feature_mapping.tsv.gz",
            pkg / "sample_selection.tsv",
        )
        # Compare preparation before estimating: preserve numerical tolerance for text serialization.
        checking = True
        original = source / "outputs/stage1_three_effects_20260903_v1/contrasts" / cid
        prep, samples, shape = check_preparation(
            target, original, pkg / "sample_selection.tsv"
        )
        audit.update(
            input_reproduced=None if prep is None else bool(prep),
            samples_reproduced=bool(samples),
            input_shape=shape,
        )
        checking = False
        with (target / "effects.log").open("w") as log:
            subprocess.run(
                [args.rscript, str(code / "effects.R"), cid],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
                env={
                    **os.environ,
                    "BIOLOGY_SERIALIZED_REFIT": "1" if historical else "0",
                },
            )
        checking = True
        comparisons = compare_effect_tables(target, source, cid)
        state = reproduction_state(prep, samples, comparisons)
        audit.update(comparisons=comparisons, reproduction=state)
    except Exception as e:
        checking = True
        audit.update(
            reproduction="ERROR", error=str(e), traceback=traceback.format_exc()
        )
        print(traceback.format_exc(), flush=True)
    audit.update(
        finished=time.time(), read_paths=sorted(reads), work_directory=str(target)
    )
    target.mkdir(parents=True, exist_ok=True)
    (target / "processing_log.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf8"
    )
    print(cid, audit["reproduction"], flush=True)
    if audit["reproduction"] in ["ERROR", "MISMATCH"]:
        raise SystemExit(1)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--contrast")
    p.add_argument("--source-root", required=True)
    p.add_argument("--work-root", required=True)
    p.add_argument("--rscript", default="Rscript")
    p.add_argument("--stage-only", action="store_true")
    run(p.parse_args())
