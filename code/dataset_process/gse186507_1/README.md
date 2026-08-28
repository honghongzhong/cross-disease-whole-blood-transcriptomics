# GSE186507: Crohn's disease versus healthy controls

## Data selected

The analysis uses the 432 Crohn's disease samples and 209 healthy controls marked `selected=TRUE` in `sample_manifest.tsv`. Ulcerative-colitis samples are excluded.

## Original data

`../raw/GSE186507_MSCCR_Blood_counts.txt.gz` is the GEO supplementary blood count matrix. Sources: [GSE186507](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186507) and [GEO supplementary files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186507/suppl/).

## Processing

`run_analysis.py` reads the selected count columns, converts counts to log2 counts-per-million, and calculates the Crohn's-disease-minus-healthy mean difference, standard-error estimate, Welch t-test p-value, and Benjamini–Hochberg q-value.

## Run

```text
python dataset/gse186507_1/process/run_analysis.py
```

The result is written to `result/GSE186507_CD_vs_HC_E_U.tsv`.
