# GSE100150 contrast

## Data selected

The quality-controlled samples in `sample_manifest.tsv` define this contrast: rows marked `Case` are the condition named by this dataset directory, and rows marked `Healthy` are controls. Other conditions from this GEO series are excluded.

## Original data

`../raw/GSE100150_series_matrix.txt.gz` is the GEO series matrix and `../raw/GPL6884.annot.gz` is the GPL6884 probe annotation. NCBI sources: [GSE100150](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100150), [GEO series files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100150/), and [GPL6884](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6884).

## Processing

`run_analysis.py` reads the selected samples, joins probe annotation, removes probes without valid gene symbols, averages probes for each gene, and computes the mean expression difference (case minus healthy), standard-error estimate, Welch t-test p-value, and Benjamini–Hochberg q-value.

## Run

From the project root:

```text
python dataset/gse100150_2/process/run_analysis.py
```

The script writes the result table to this dataset's `result` directory.
