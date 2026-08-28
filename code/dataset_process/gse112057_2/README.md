# GSE112057: Ulcerative Colitis vs Healthy Controls

## Purpose

This folder contains the files used to process the GSE112057 blood-expression dataset for people with ulcerative colitis and healthy controls.

## Input files

The `raw` folder contains the original GEO files:

- `GSE112057_RawCounts_dataset.txt.gz` — deposited gene-count table.
- `GSE112057_family.soft.gz` — GEO sample and study metadata.

NCBI download links:

- Series page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112057
- Count table: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112057/suppl/GSE112057_RawCounts_dataset.txt.gz
- GEO metadata: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112057/soft/GSE112057_family.soft.gz

`sample_metadata.tsv` records the sample labels used for selection. `dataset_protocol.md` records the study-specific selection rules.

## Selection and processing

The analysis selects 15 samples labelled Ulcerative Colitis and 12 samples labelled Control. All selected samples are whole blood. Sample names are matched to count-table columns before analysis.

The script groups repeated gene symbols, removes genes with fewer than three samples having at least 10 counts, fits a two-group negative-binomial model, calculates log2 fold change, standard error, p value and false-discovery adjusted p value, and maps gene symbols to unique Entrez identifiers using the project annotation database.

## Run

Install Python 3.10 or newer and the packages listed in `requirements.txt`, then run from the project root:

```powershell
python -m pip install -r .\dataset\gse112057_2\process\requirements.txt
python .\dataset\gse112057_2\process\run_analysis.py
```

The script reads files relative to its own location and writes `result/GSE112057_Ulcerative_Colitis_vs_HC_E_U.tsv`.

## Output columns

The output contains `gene_id`, `effect`, `se`, `p_value`, `adjusted_p_value`, `n_case`, `n_control`, `contrast_id`, and `analysis_method`.

## Reproducibility notes

The expected group sizes are 15 case samples and 12 control samples. The earlier result is retained beside the current output for comparison. The current script uses PyDESeq2 and the annotation database already included in the project tools.
