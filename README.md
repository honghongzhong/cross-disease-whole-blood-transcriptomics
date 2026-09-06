# Shared low-dimensional structure across diseases in whole-blood transcriptomes

Version 2.0.0 contains the code, derived data and results for the strengthened whole-blood study: 48 disease–healthy control comparisons from 43 GEO studies, 36 training comparisons, 12 evaluation comparisons and a 6144-gene panel.

## Reproduce the paper

Download this release's source archive and extract it. From its root, use Python 3.11–3.13:

```sh
python -m pip install -r analysis/requirements.txt
python analysis/run_analysis.py --data data --output new_results
```

The entry point downloads the specified Hallmark 2025.1.Hs annotation, runs dimension selection, reconstruction across studies and excluded disease categories, permutation references, enrichment, resampling, sensitivity and residual analyses, and draws the five figures. It uses a supported CUDA GPU when available and otherwise the CPU. See `analysis/README.md` for expected results and the CPU/GPU random-number distinction.

## Contents

| Directory | Contents |
|---|---|
| `analysis/` | Main analysis and figure scripts, requirements and fixed input lists |
| `data/contrasts/` | All 144 derived effect-and-uncertainty tables for the 48 comparisons |
| `results/` | Machine-readable study lists and reference results |
| `figures/` | Five figures in vector PDF and high-resolution PNG |
| `processing/` | Shared upstream code, study-specific sample selections and GEO download instructions |

The native-scale selected dimensions are 5, 4 and 3 for standardized, log2 and shrunken-log2 effects. Median study reconstruction R² is approximately 0.429, 0.421 and 0.369. Category-excluded medians are approximately 0.429, 0.400 and 0.344. These are representations of observed study-level effects; the combination scores are fitted to each test profile.

The archived `results/` tables are reference outputs, not computation inputs. New runs write to a separate directory. The full primary pipeline was rerun for this release's analysis and agreed with the manuscript results. Upstream sample-level processing was additionally rerun for one microarray and one RNA-seq comparison; all 48 upstream comparisons were not rerun during packaging.

Public GEO measurements and large public annotation files are not duplicated. `processing/GEO_downloads.csv` names the precise inputs and source locations; `processing/public_annotations.csv` and `download_annotations.py` provide public platform resources. The compact mappings retained in the package preserve the gene mapping used in the study.

## Version and citation

Use this version for the manuscript “Shared low-dimensional structure across diseases in whole-blood transcriptomes”. Version 1.0.0 remains available under its original Git tag and archived release and corresponds to the earlier EU57 study; its input set and evaluation design differ.

Version-specific Zenodo archive: [10.5281/zenodo.22446550](https://doi.org/10.5281/zenodo.22446550). Cite this record for v2.0.0. Author: Hongzhong Hong, South China University of Technology; ORCID 0009-0005-2114-143X.

Author-written code: MIT. Author-created tables and figures: CC BY 4.0. Third-party data and annotations retain their applicable terms.
