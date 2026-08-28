# GSE22229_renal_tolerance_vs_HC

## Data selected

The dataset-specific sample selection and inclusion decisions are recorded in `process/sample_manifest.tsv` when available. The contrast represented here is `GSE22229_renal_tolerance_vs_HC`.

## Original data

Original files are stored in `raw/`. NCBI sources: [NCBI GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22229) and [GEO supplementary files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE222nnn/GSE22229/suppl/).

## Processing

`run_analysis.py` reads the source gene-level table already stored in `result/`, standardizes its effect and uncertainty fields as `E` and `U`, removes invalid rows, and writes `result/processed_E_U.tsv`. All paths are relative to this dataset directory.

## Run

```text
python dataset/gse22229_1/process/run_analysis.py
```

The generated table is `result/processed_E_U.tsv`.
