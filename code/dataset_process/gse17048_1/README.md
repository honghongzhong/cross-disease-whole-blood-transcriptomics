# GSE17048: multiple sclerosis versus healthy controls

## Data selected

The analysis uses the 99 multiple-sclerosis samples and 45 healthy controls with `selected=TRUE` in `sample_manifest.tsv`.

## Original data

`../raw/GSE17048_series_matrix.txt.gz` is the NCBI GEO series matrix. Sources: [GSE17048](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17048) and [GEO supplementary files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE017nnn/GSE17048/suppl/).

## Processing

`run_analysis.py` reads the selected GEO columns, averages duplicate probe identifiers, and calculates the multiple-sclerosis-minus-healthy mean difference, standard-error estimate, Welch t-test p-value, and Benjamini–Hochberg q-value.

## Run

```text
python dataset/gse17048_1/process/run_analysis.py
```

The result is written to `result/GSE17048_MS_vs_HC_E_U.tsv`.
