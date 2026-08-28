# GSE124548: cystic fibrosis baseline versus healthy controls

## Data selected

The analysis uses 20 baseline cystic-fibrosis samples and 20 non-CF healthy controls listed in `sample_manifest.tsv`. Post-treatment samples are excluded.

## Original data

`../raw/GSE124548_raw_baseline_counts.tsv.gz` is the baseline count table extracted from the GEO supplementary data, and `../raw/GSE124548_AllData_170308_RNAseq_Kopp_Results.xlsx` is the original supplementary workbook. NCBI sources: [GSE124548](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124548) and [supplementary files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124548/suppl/).

## Processing

`run_analysis.py` selects baseline samples, converts counts to log2 counts-per-million, and calculates the CF-minus-healthy mean difference, standard-error estimate, Welch t-test p-value, and Benjamini–Hochberg q-value.

## Run

```text
python dataset/gse124548_1/process/run_analysis.py
```

The result is written to `result/GSE124548_CF_baseline_vs_HC_E_U.tsv`.
