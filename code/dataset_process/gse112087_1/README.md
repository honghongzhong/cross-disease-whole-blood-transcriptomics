# GSE112087: systemic lupus erythematosus versus healthy controls

## Data selected

The analysis uses the 31 case donors and 29 healthy-control donors with `include_locked=TRUE` in `sample_manifest.tsv`. Each donor may have two sequencing lanes; the lanes listed in `lane_gsm_ids` are combined for that donor.

## Original data

`../raw/GSE112087_RAW.tar` is the NCBI GEO supplementary archive containing the original lane-level count files. Source: [GSE112087](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112087) and [GEO supplementary archive](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112087/suppl/GSE112087_RAW.tar).

## Processing

`run_analysis.py` reads the selected lane files directly from the archive, sums lanes belonging to each donor, converts counts to log2 counts-per-million, and calculates the case-minus-control mean difference, standard-error estimate, Welch t-test p-value, and Benjamini–Hochberg q-value.

## Run

From the project root:

```text
python dataset/gse112087_1/process/run_analysis.py
```

The result is written to `result/GSE112087_SLE_vs_HC_E_U.tsv`.
