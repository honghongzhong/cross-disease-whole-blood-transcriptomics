# Study-specific expression processing

These scripts rebuild the 48 comparisons from the GEO-deposited sample-level inputs named in `GEO_downloads.csv`. The file lists the GEO page, input filename and expected path below your source directory. Download those inputs from their GEO records and place them at the listed paths. Some deposited inputs are normalized expression or counts; the workflow does not claim to reprocess every study from FASTQ.

Install Python dependencies in `_shared/requirements.txt` and R with limma, edgeR, ashr, affy, jsonlite and data.table. The original R package versions are in `_shared/R_sessionInfo_original.txt`.

Run `python download_annotations.py --source-root SOURCE` to obtain the public platform resources. `public_annotations.csv` also gives their direct URLs. The compact feature mappings and cached identifiers retained here preserve the gene mapping used in the paper; large public platform dumps and expression measurements are not included.

To process one comparison, run:

```text
python GSE34404_malaria_vs_HC/run.py --source-root SOURCE --work-root NEW_OUTPUT --rscript PATH_TO_RSCRIPT
```

Use a fresh work directory. Results are written under
`NEW_OUTPUT/outputs/stage1_three_effects_20260903_v1/contrasts/<comparison>/`.
The single shared processing workflow is in `_shared/workflow/`; comparison
folders supply metadata, frozen feature mappings, and a `sample_selection.tsv`
that declares the case and control samples in order. The shared workflow
checks this declaration against the sample IDs selected from source data before
effect estimation.

The shared code is organized as follows:

- `_shared/reproduce.py`: command-line orchestration and audit log.
- `_shared/staging.py`: stage bundled resources and link external inputs.
- `_shared/workflow_setup.py`: copy the shared workflow into an isolated work tree.
- `_shared/verification.py`: compare regenerated inputs and effect tables.
- `_shared/workflow/`: expression preparation, R effect estimation, and output utilities.
- `_shared/workflow/adapters/`: direct study modules, source-format readers, and gene-ID mapping.
- `_shared/resources/samples/`: original study metadata and sample evidence used by some source readers. Every comparison's final sample selection is declared in its own `sample_selection.tsv`.
- `_shared/resources/annotations/`: bundled feature and platform mappings.
- `_shared/resources/audit/`: compact source index and review inventory needed by the workflow. These JSON files use English field names and relative paths; the historical full audit snapshot is not bundled.
- `_shared/inputs.json`: maps each bundled resource to its required runtime location and lists external inputs.

The release contains one shared implementation under `_shared/workflow/`.
Historical source directory names in `inputs.json` describe the input layout
expected by the readers; they are created only in the separate work directory.

Each comparison has a `meta.json` with `input` (study and input details, including relative source paths and required external files) and `qc` (the original processing quality-control summary). The metadata contains only ASCII text and no local absolute paths. Historical review reasons and notes are translated into English. The GEO accession and relative evidence paths identify the underlying sources.

The original source-layout names are retained so existing sample loaders can read the correct inputs. They are folder names, not additional experiments required by this paper. The primary analysis package can be run directly from the supplied derived effect tables without repeating this stage.
