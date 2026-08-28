# GSE165082: Parkinson disease versus comparison controls

## Data selected

The analysis uses 12 samples labelled `PD` and 14 samples labelled `CC` in `sample_manifest.tsv`. All selected samples are whole-blood profiles.

## Original data

`../raw/GSE165082_PD-CC.counts.txt.gz` is the GEO count matrix. NCBI sources: [GSE165082](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165082) and [supplementary files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165082/suppl/).

## Processing

`run_analysis.py` converts counts to log2 counts-per-million and calculates the PD-minus-control mean difference, standard-error estimate, Welch t-test p-value, and Benjamini–Hochberg q-value.

## Run

```text
python dataset/gse165082_1/process/run_analysis.py
```

The result is written to `result/GSE165082_PD_vs_HC_E_U.tsv`.
