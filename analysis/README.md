# Shared whole-blood structure across diseases

This code runs the analyses in the BMC Genomics manuscript: uncentred SVD, training-based dimension selection, study and disease-category reconstruction, gene-permuted references, Hallmark enrichment, group resampling, heme sensitivity and pathway residuals.

## Run

1. Use Python 3.11–3.13. Install `pip install -r requirements.txt`. PyTorch uses a supported CUDA GPU if available and otherwise uses the CPU. The original numerical calculations used double precision and PyTorch 2.6.0 with CUDA 12.4.
2. The repository includes `../data/contrasts/` with all derived effect tables.
3. Run `python run_analysis.py --data ../data --output ../new_results`.

The program downloads human Hallmark **2025.1.Hs, Entrez identifiers** from the Broad Institute. To use a downloaded copy, add `--hallmark path/to/h.all.v2025.1.Hs.entrez.gmt`.

The main results are `results/analysis/selected_metrics.tsv`, `selected_domain_transfer.tsv`, `dimension_decision.tsv`, `selected_enrichment.tsv`, `selected_sensitivity.tsv` and `pathway_residuals.tsv`. The five publication figures are in `results/figures` as vector PDF and 450-dpi PNG. To redraw them: `python make_figures.py --results results/analysis --output figures`.

The native-scale selected dimensions should be 5, 4 and 3. Median study-disjoint R² should be approximately 0.429014, 0.420562 and 0.368671; category-excluded medians approximately 0.428970, 0.399688 and 0.343879. Random permutations can differ between CPU and GPU backends even with the same seed, particularly near a rank threshold. Deterministic SVD/projection results should agree to numerical precision.

The first stage retains the historical K=5 calculations needed by subsequent selection and sensitivity steps. Use the `selected_*` tables for the manuscript's representation-specific dimensions. The scripts use the fixed 36/12 study lists; they do not select a new cohort. Test scores are fitted from each observed test vector, so these results concern representation, not unmeasured-gene prediction.

## Source data

The 144 supplied tables are derived analysis inputs, not raw GEO measurements. `GEO_downloads.csv` lists each source study, its GEO page and the input filename used in processing. Raw measurements remain at GEO. The separate processing-code archive contains study-specific sample assignments, transformations and scripts for rebuilding these effects. GEO studies can provide processed matrices rather than raw reads; use the input named for each study.

Hallmark download: https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/h.all.v2025.1.Hs.entrez.gmt

Code: MIT License. Author-created input lists and figures: CC BY 4.0. Public data and annotations retain their original terms.
