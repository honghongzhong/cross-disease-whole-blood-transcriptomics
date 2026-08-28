# GSE48000: high-risk venous thromboembolism versus healthy controls

## Data selected

The analysis uses the 40 high-risk VTE samples and 25 healthy controls listed in `sample_manifest.tsv`.

## Original data

The raw signal matrix is stored in `../raw/`. Download source: [GSE48000](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48000) and [GEO supplementary files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48000/suppl/).

`../raw/GSE48000_non_normalized_set1.txt.gz` contains the expression signals and `../raw/GPL10558.annot.gz` contains probe annotation. Sources: [GSE48000](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48000), [GEO supplementary files](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE048nnn/GSE48000/suppl/), and [GPL10558](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558).

## Processing

`run_analysis.py` log2-transforms and quantile-normalizes the signals, maps probes to gene symbols, averages duplicate probes, and calculates the VTE-minus-healthy mean difference, standard-error estimate, Welch t-test p-value, and Benjamini–Hochberg q-value.

## Run

```text
python dataset/gse48000_1/process/run_analysis.py
```

The result is written to `result/GSE48000_high_risk_VTE_vs_HC_E_U.tsv`.
