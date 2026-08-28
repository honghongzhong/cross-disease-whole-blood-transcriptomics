# GSE196822: moderate COVID-19 versus healthy controls

## Data selected

The analysis uses 10 moderate COVID-19 samples and 9 healthy samples. The sample definition is recorded in `sample_manifest.tsv`; all other conditions in the series are excluded.

## Original data

The raw count matrix is stored in `../raw/`. Download source: [GSE196822](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196822) and [GEO supplementary files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196822/suppl/).

`../raw/GSE196822_Raw_counts_matrix.csv.gz` is the GEO supplementary count matrix. Sources: [GSE196822](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196822) and [GEO supplementary files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE196nnn/GSE196822/suppl/).

## Processing

`run_analysis.py` selects the moderate and healthy samples, converts counts to log2 counts-per-million, and computes the case-minus-healthy mean difference, standard-error estimate, Welch t-test p-value, and Benjamini–Hochberg q-value.

## Run

```text
python dataset/gse196822_1/process/run_analysis.py
```

The result is written to `result/GSE196822_Moderate_vs_HC_E_U.tsv`.
