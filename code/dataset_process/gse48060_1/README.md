# GSE48060: acute myocardial infarction versus healthy controls

## Data selected

The analysis uses the 31 patient samples labelled `case` and 21 normal controls labelled `healthy_control` in `sample_manifest.tsv`. Recurrent/non-recurrent subgroup labels are not used in this contrast.

## Original data

`../raw/GSE48060_series_matrix.txt.gz` is the NCBI GEO series matrix. Sources: [GSE48060](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48060) and [GEO supplementary files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE048nnn/GSE48060/suppl/).

## Processing

`run_analysis.py` reads the selected GEO columns, averages duplicate probe identifiers, and calculates the case-minus-control mean difference, standard-error estimate, Welch t-test p-value, and Benjamini–Hochberg q-value.

## Run

```text
python dataset/gse48060_1/process/run_analysis.py
```

The result is written to `result/GSE48060_AMI_vs_HC_E_U.tsv`.
