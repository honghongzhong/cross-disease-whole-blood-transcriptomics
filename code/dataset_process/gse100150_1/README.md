# GSE100150: HIV samples versus healthy controls

## Data selected

The analysis uses the quality-controlled samples listed in `sample_manifest.tsv`: 28 HIV samples and 35 healthy controls. Samples from the other conditions in the series are not included.

## Original data

Files downloaded from or derived from NCBI GEO are in `../raw/`: `GSE100150_series_matrix.txt.gz` and `GPL6884.annot.gz`.

NCBI sources: [GSE100150](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100150), [GEO series files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100150/), and [GPL6884](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6884).

## Processing

`run_analysis.py` reads the selected samples, joins probe annotation, removes probes without a valid gene symbol, averages probes mapping to the same gene, and calculates the mean expression difference (HIV minus healthy), a standard-error estimate, Welch t-test p-values, and Benjamini–Hochberg q-values.

## Reproduce the result

From the project root, install the packages in `process/requirements.txt` if needed and run:

```text
python dataset/gse100150_1/process/run_analysis.py
```

The output is `../result/GSE100150_HIV_vs_HC_E_U.tsv`.
