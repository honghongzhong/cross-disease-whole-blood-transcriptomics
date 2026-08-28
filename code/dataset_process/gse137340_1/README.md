# GSE137340: sepsis day 1 versus healthy controls

## Data selected

The analysis uses the 21 sepsis day-1 samples (`analysis_group=Sepsis_D1`) and 12 healthy controls (`analysis_group=Healthy_Control`) with `selected=TRUE` in `sample_manifest.tsv`. Later time points are excluded.

## Original data

`../raw/GSE137340_series_matrix.txt.gz` is the NCBI GEO series matrix. Sources: [GSE137340](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137340) and [GEO supplementary files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137340/suppl/).

## Processing

`run_analysis.py` reads the selected GEO columns, averages duplicate probe identifiers, and calculates the sepsis-minus-healthy mean difference, standard-error estimate, Welch t-test p-value, and Benjamini–Hochberg q-value.

## Run

```text
python dataset/gse137340_1/process/run_analysis.py
```

The result is written to `result/GSE137340_Sepsis_Day1_vs_HC_E_U.tsv`.
