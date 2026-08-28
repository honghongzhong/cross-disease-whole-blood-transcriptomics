# GSE112057 Protocol and E/U Lock v1

## Status

**Passed for five candidate-extension E/U contrasts using the accession's explicit Control group:** Crohn's Disease, Ulcerative Colitis, Polyarticular JIA, Oligoarticular JIA, and Systemic JIA, each versus Control.

## Frozen scope

- Accession: GSE112057; 202 samples, all `Whole Blood` and GPL11154.
- Sample-level diagnosis counts: Crohn's Disease 60, Ulcerative Colitis 15, Polyarticular JIA 46, Oligoarticular JIA 43, Systemic JIA 26, Control 12.
- Sample titles are unique and the family metadata contains no repeated timepoint field for these samples; the accession is treated as one independence group.
- The 12 Control samples are used as the same shared-control group across the five contrasts.

## Data and model lock

- Input: `GSE112057_RawCounts_dataset.txt.gz`, 26,305 gene-symbol rows and 202 integer sample columns.
- Filter: count >=10 in at least 3 samples within each binary contrast.
- Model: DESeq2 1.52.0, group-only `~ disease_group`; no age/sex adjustment was added because the Control samples have `NA` for those fields and no validated balanced covariate lock exists.
- E: named disease subgroup minus Control log2 fold change. U: DESeq2 standard error for the same coefficient.
- Gene harmonization: one-to-one numeric Entrez IDs through `org.Hs.eg.db` after symbol-level duplicate collapse.

## Output inventory

The five E/U files and `GSE112057_run_summary.tsv` record the following validated row counts: Crohn's 13,000; Ulcerative Colitis 12,821; Polyarticular JIA 12,992; Oligoarticular JIA 13,061; Systemic JIA 12,895. All retained E/U values are finite with positive standard errors. Sample QC records library size and detected-gene count for all 202 samples; no outcome-informed sample was removed.

## Current verification and hashes

The five E/U files were re-read on 2026-08-22 and each passed unique Entrez, finite E/U, and U>0 checks. Current artifact SHA-256 values are: `GSE112057_family.soft.gz` `958E2983136CBBFC9E6ECA96F18C82EA9DD97A628554B259030DF7611B4CBFBE`; `GSE112057_RawCounts_dataset.txt.gz` `9DFE0B66E71FBB47233148BCA8F73021A54CE8F8742AFFCEB16FB246B6F65CAB`; `GSE112057_sample_metadata.tsv` `FBE968758345CB964FD9AA263608D38B1E5C390FE14283E45D3210D8F8815A89`; `GSE112057_run_summary.tsv` `F00A66859D2E733BD7286909558F0B4EA5C91BAD649FBAE75C633685D2E852C`.

## Interpretation boundary

These are five locked within-accession E/U rows sharing one Control cohort. They add one independence group, do not change the frozen core_24 denominator, and are not cross-disease confirmatory conclusions. Because the Control covariate fields are incomplete and the JIA subgroups have different sex distributions, downstream interpretation must preserve the group-only model limitation.
