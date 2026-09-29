"""Download public platform resources into the source layout used by preprocessing."""

from pathlib import Path
import argparse, csv, urllib.request

ap = argparse.ArgumentParser()
ap.add_argument("--source-root", type=Path, required=True)
args = ap.parse_args()
for row in csv.DictReader(
    (Path(__file__).parent / "public_annotations.csv").open(encoding="utf8")
):
    dest = args.source_root / row["relative"]
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        print("Downloading", row["platform"], dest.name, flush=True)
        urllib.request.urlretrieve(row["url"], dest)
