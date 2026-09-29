"""Bounded concurrent R workers; one contrast owns its output directory."""

from pathlib import Path
from config import work_root, rscript_command
import concurrent.futures as cf, json, subprocess, time, sys

HERE = Path(__file__).parent
ROOT = work_root() / "outputs/stage1_three_effects_20260903_v1"
RS = rscript_command()
skip100 = "--skip100150" in sys.argv


def run(d):
    with (d / "effect_run.log").open("w") as log:
        p = subprocess.run(
            [RS, str(HERE / "effects.R"), d.name], stdout=log, stderr=subprocess.STDOUT
        )
    return d.name, p.returncode, (d / "qc.json").exists()


seen = set()
results = []
with cf.ThreadPoolExecutor(max_workers=3) as pool:
    jobs = {}
    while True:
        for p in sorted((ROOT / "contrasts").glob("*/input.json")):
            d = p.parent
            if (
                d.name in seen
                or (d / "qc.json").exists()
                or (skip100 and d.name.startswith("GSE100150_"))
            ):
                continue
            if len(jobs) >= 3:
                break
            seen.add(d.name)
            jobs[pool.submit(run, d)] = d
        if not jobs:
            break
        done, _ = cf.wait(jobs, timeout=30, return_when=cf.FIRST_COMPLETED)
        for j in done:
            d = jobs.pop(j)
            r = j.result()
            results.append(r)
            print(r, flush=True)
            (ROOT / "effect_worker_results.json").write_text(
                json.dumps(results, indent=2)
            )
