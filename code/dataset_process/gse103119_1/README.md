# GSE103119: pneumonia versus healthy controls

## Data selected

The analysis uses the 152 samples labelled `case` and 20 samples labelled `healthy_control` in `sample_manifest.tsv`. All selected samples are whole-blood samples from this GEO series.

## Original data

`../raw/GSE103119_non-normalized.txt.gz` contains the GEO expression signals and detection p-values. `../raw/GPL10558.annot.gz` contains the Illumina HumanHT-12 v4 probe annotation. NCBI sources: [GSE103119](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103119), [GEO supplementary files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103119/), and [GPL10558](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558).

## Processing

`run_analysis.py` selects the manifest samples, filters probes by detection p-value, maps probes to gene symbols, averages probes per gene, and calculates the case-minus-control mean difference, its standard-error estimate, Welch t-test p-value, and Benjamini–Hochberg q-value.

## Run

From the project root:

```text
python dataset/gse103119_1/process/run_analysis.py
```

The result is written to `result/GSE103119_Pneumonia_vs_HC_E_U.tsv`.
