# EU57 cross-disease transcriptomic model

Reproducibility package for *Uncertainty-weighted low-rank modelling reveals
shared transcriptional structure across 57 disease-control blood contrasts*.

The package contains the standardized gene-level inputs used in the paper,
portable analysis code, fixed pathway libraries, reference results, manuscript
figure code, and accession-specific processing records. No script relies on a
developer-specific directory. Paths stored in public manifests and generated
metadata are relative to the package root or to their containing result folder.

## What can be reproduced

The package supports two levels of verification:

1. **Direct reproduction of the reported analysis.** The 57 packaged
   standardized contrasts are sufficient to repeat disease-domain validation,
   rank selection, accession sensitivity, final-model fitting, stability,
   residual analysis, pathway analysis, and gene-panel sensitivity.
2. **Audit of upstream processing.** `code/dataset_process/` contains the
   available processing script and sample-selection record for each source
   dataset directory. Raw and normalized GEO downloads are not redistributed;
   rerunning this upstream layer therefore requires the corresponding GEO files.

The reproducible analysis unit is a disease-control contrast, not an individual
patient. Each input reports Hedges' g for case minus control and its sampling
standard error. Missing genes remain missing and are never replaced with zero.

## Package layout

```text
code/modeling/                 Main model and sensitivity analyses
code/dataset_process/          Accession-specific processing and selection records
code/manuscript/               Figure-generation code
data/standardized_eu/          The 57 standardized input files
data/input_manifest.tsv        Portable input paths, hashes and source groups
metadata/input_audit.tsv       Sample sizes and input-processing metadata
metadata/disease_class_map.tsv Disease labels and the 12 evaluation domains
resources/annotation_sources/  Fixed Hallmark and Reactome 2022 GMT files
reference_results/             Formal results reported in the manuscript
validation/smoke_results/      Packaged reduced-run evidence; not reportable
manuscript_figures/            Regenerated manuscript figures
outputs/                       New user-generated runs; absent until execution
```

`reference_results/` is an immutable comparison target. Do not write a new run
into that directory.

## Requirements

- Python 3.13.5
- Approximately 4 GB of free memory
- Approximately 1 GB of free disk space in addition to the extracted package
- The exact Python package versions in `requirements-lock.txt`

The reported formal run used Windows 11 and completed in approximately 42
minutes. Runtime depends on processor, BLAS implementation and operating system.

### Conda

From the extracted package root:

```bash
conda env create -f environment.yml
conda activate eu57-paper
```

### Python virtual environment

Linux or macOS:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-lock.txt
```

Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements-lock.txt
```

## Quick verification

First validate the extracted files:

```bash
python validate_package.py
```

A valid release reports 57 unique inputs, 57 dataset-processing directories,
successful Python compilation, no reference hash errors,
`"absolute_path_files": []`, and `"passed": true`.

Then run the reduced executability check:

```bash
python run_all.py --mode smoke
```

This writes to `outputs/smoke_results/`. It uses fewer folds, masks and
iterations and creates `SMOKE_TEST_ONLY.txt`. Its numerical estimates must not
be reported as scientific results. Every run requires a new or empty output
directory; for a repeated check, supply another location:

```bash
python run_all.py --mode smoke --output outputs/smoke_results_2
```

## Reproduce the formal analysis

Run the manuscript settings from the package root:

```bash
python run_all.py --mode formal
```

The wrapper performs the following operations in order:

1. verifies every standardized input against `data/input_manifest.tsv`;
2. repeats purged leave-one-disease-domain-out validation with training-only
   gene and rank selection;
3. selects the final rank by grouped internal validation and the
   one-standard-error rule;
4. runs leave-one-accession-out and 500/1,000/2,000-gene sensitivity analyses;
5. fits the final model and performs source-group bootstrap stability analysis;
6. calculates cross-fitted residual and Hallmark/Reactome enrichment results.

New formal results are written only to `outputs/formal_results/`. The fixed GMT
files are copied into the output before enrichment, so rerunning the analysis
does not depend on a later online library version.

## Compare a formal run with the reference

After the formal workflow finishes:

```bash
python compare_results.py --candidate outputs/formal_results
```

The comparison checks the manuscript-defining result files and writes
`outputs/formal_results/reference_comparison.json`. With the pinned environment
and matching numerical libraries, all entries are expected to report
`"sha256_identical": true`. Different operating systems or BLAS builds can
produce very small floating-point differences even when scientific conclusions
are unchanged; such a run should be reviewed numerically rather than described
as byte-identical.

## Regenerate manuscript figures

The five figures can be regenerated directly from the packaged formal results:

```bash
python code/manuscript/make_figures.py
```

The command writes PNG and vector PDF files to `manuscript_figures/`. It does
not read the reduced smoke-test results.

## Reported settings

- Random seed: `20260828`
- Disease-control contrasts: 57
- GEO accessions: 49
- Independent source groups: 48
- Disease domains: 12
- Training-selected genes: 1,000 with at least 90% training coverage
- Candidate ranks: 1-6
- Held-out target genes: 20%, repeated five times per disease-domain test
- Ridge penalty: `0.001`
- Source-group bootstrap repetitions: 1,000
- Final-model stability repetitions: 200
- Selected final rank: 5

## Result map

| Manuscript component | Reproduction artifact |
|---|---|
| Primary pooled cross-validated result | `reference_results/primary_result.json` |
| Disease-domain performance | `reference_results/outer_class_metrics.tsv` |
| Rank selected in each outer fold | `reference_results/outer_selected_ranks.tsv` |
| Final all-data rank selection | `reference_results/final_rank_selection/` |
| Accession sensitivity | `reference_results/leave_one_accession_out_summary.json` |
| Gene-panel sensitivity | `reference_results/panel_size_sensitivity.tsv` |
| Final scores and loadings | `reference_results/final_contrast_factor_scores.tsv`, `reference_results/extended_gene_loadings.tsv.gz` |
| Factor stability | `reference_results/factor_stability_bootstrap.tsv` |
| Hallmark and Reactome results | `reference_results/pathway_rank_enrichment.tsv` |
| Cross-fitted residual effects | `reference_results/disease_class_crossfit_residual_meta.tsv.gz` |

## Integrity and provenance

- `PACKAGE_MANIFEST.tsv` gives the byte size and SHA-256 digest of every
  released file except the manifest itself.
- `data/input_manifest.tsv` gives the SHA-256 digest of each standardized input.
- `reference_results/artifact_manifest.tsv` covers all formal reference results.
- `metadata/dataset_process_code_manifest.tsv` covers the packaged upstream
  processing records.

The full formal analysis starts from the standardized contrasts. The upstream
GEO audit materials document how those contrasts were produced but do not make
the package a completely offline raw-GEO-to-paper reconstruction.

## Licensing and citation

Author-written code is released under the MIT License (`LICENSE`).
Author-created tables, results and figures are released under CC BY 4.0
(`LICENSE_DATA.md`). Fixed third-party gene-set libraries retain their original
terms; see `THIRD_PARTY_NOTICES.md`.

Please cite the accompanying manuscript and the original GEO studies when using
the packaged data products. Version and author metadata are provided in
`CITATION.cff`; a DOI can be added by the archival service when this exact
release is deposited.
