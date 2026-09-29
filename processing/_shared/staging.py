"""Stage bundled resources and external inputs in an isolated work tree."""

from pathlib import Path
import json
import re
import shutil
import time

BASE = Path(__file__).resolve().parent.parent
SOURCE_PREFIX = re.compile(
    r"(?i)\b[A-Z]:[\\/]+[^\\/]+[\\/]+"
    r"(?=(?:dataset|outputs|census_step0|phase1b_[^\\/]+|tmp|tools|"
    r"effect-uncertainty|modeling|GPL\d+)[\\/])"
)


def rewrite(s, root):
    """Map historical Windows source paths into the isolated work tree."""
    return SOURCE_PREFIX.sub(root.as_posix() + "/", s)


def stage(source, work):
    """Copy bundled source data and link available external inputs into work."""
    work.mkdir(parents=True, exist_ok=True)
    rows = json.loads((BASE / "_shared/inputs.json").read_text(encoding="utf8"))
    for r in rows:
        dest = work / r["relative"]
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists() and not r["bundled"]:
            continue
        if r["bundled"]:
            src = BASE / "_shared/resources" / r["bundled_path"]
            if src.suffix in [
                ".py",
                ".R",
                ".r",
                ".json",
                ".tsv",
                ".csv",
                ".yaml",
                ".yml",
                ".md",
                ".ps1",
            ]:
                try:
                    dest.write_text(
                        rewrite(src.read_text(encoding="utf-8-sig"), work),
                        encoding="utf8",
                    )
                except UnicodeError:
                    shutil.copy2(src, dest)
            else:
                shutil.copy2(src, dest)
        else:
            src = source / r["relative"]
            if src.is_file():
                dest.symlink_to(src)
    # Libraries are environment dependencies, not source expression data.
    for rel in [
        "tools/stage1_three_effects_20260903/library",
        "tools/r461_bioc323_gse153315/library",
    ]:
        lib = work / rel
        origin = source / rel
        if not lib.exists() and origin.exists():
            lib.parent.mkdir(parents=True, exist_ok=True)
            lib.symlink_to(origin, target_is_directory=True)
        elif origin.exists():
            for p in origin.iterdir():
                if not (lib / p.name).exists():
                    (lib / p.name).symlink_to(p, target_is_directory=p.is_dir())
    (work / "STAGED.json").write_text(
        json.dumps({"source": str(source), "created": time.time()})
    )
