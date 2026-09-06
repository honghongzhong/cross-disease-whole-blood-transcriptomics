# Study-specific expression processing

These scripts rebuild the 48 comparisons from the GEO-deposited sample-level inputs named in `GEO_downloads.csv`. The file lists the GEO page, input filename and expected path below your source directory. Download those inputs from their GEO records and place them at the listed paths. Some deposited inputs are normalized expression or counts; the workflow does not claim to reprocess every study from FASTQ.

Install Python dependencies in `_shared/requirements.txt` and R with limma, edgeR, ashr, affy, jsonlite and data.table. The original R package versions are in `_shared/R_sessionInfo_original.txt`.

Run `python download_annotations.py --source-root SOURCE` to obtain the public platform resources. `public_annotations.csv` also gives their direct URLs. The compact feature mappings and cached identifiers retained here preserve the gene mapping used in the paper; large public platform dumps and expression measurements are not included.

To process one comparison, run:

```text
python GSE34404_malaria_vs_HC/run.py --source-root SOURCE --work-root NEW_OUTPUT --rscript PATH_TO_RSCRIPT
```

Use a fresh work directory. The output is under `NEW_OUTPUT/outputs/stage1_three_effects_20260903_v1/contrasts/<comparison>/`, with separate Hedges g, log2FC and shrunken-log2FC effect tables. The scripts preserve the study-specific transformations and serialization used to create the manuscript's effects. Shared pipeline code is stored once in `_shared/pipelines`; each comparison specifies its pipeline and sample selection.

The original source-layout names are retained so existing sample loaders can read the correct inputs. They are folder names, not additional experiments required by this paper. The primary analysis package can be run directly from the supplied derived effect tables without repeating this stage.
