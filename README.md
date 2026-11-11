# Plasma metabolomic signatures of plant and animal protein and type 2 diabetes

Analysis code for elastic-net metabolomic signatures of plant protein, animal protein, and the plant-to-animal protein ratio in the Nurses' Health Study (NHS), NHSII, and Health Professionals Follow-up Study (HPFS), with replication in the Women's Health Initiative (WHI), and prospective associations with incident type 2 diabetes.

This repository contains scripts for the **primary analyses** reported in the manuscript (Cox models, metabolome-wide associations, biomarker partial correlations, dominance analysis, sequential R², figures, and supplementary tables). It does **not** include sensitivity scripts for elastic-net hyperparameters, selection-frequency reporting, or imputation diagnostics.

## Requirements

- R (≥ 4.x) with packages used in the scripts (e.g. `survival`, `dplyr`, `tidyverse`, `glmnet`, `domir`, `ggplot2`, `patchwork`, `tableone`, `VennDiagram`, as needed per script).
- Access to cohort analysis files on the Channing/Linux analysis environment (paths below are not public).

## Path configuration

Edit `R/00_paths.R` before running any script. Replace the placeholders with your local paths:

```r
repo <- "[path]"           # project root (Enet/, Endpoints/, Review/, Heatmap/, …)
biomarker_rdata <- "[path]" # full path to biomarker.RData (labs object `pool`)
met_names_csv <- "[path]"   # full path to metabolite name lookup CSV
```

Each analysis script begins with `source("00_paths.R")` and builds file locations with `file.path(repo, "…")`. Signature scores and imputed metabolomics data are produced upstream (elastic-net outputs and endpoint tables); this repository assumes those files exist under `repo`.

## Suggested run order (discovery)

1. `R/venn_figure1.R` — signature metabolites, Venn diagram, heatmap inputs (uses saved elastic-net results).
2. `R/table1.R` — Table 1.
3. `R/metabolite_t2d_cox_ppp.R` — per-metabolite Cox models for T2D (supplementary heatmap).
4. `R/heatmap_full_ppp_export.R` — full signature × protein × food × T2D table.
5. `R/script.R` — signature and self-reported protein Cox models (Table 2, supplementary self-report, controls-only sensitivity).
6. `R/food_contribution_signatures.R` — food-group dominance (`Review/food_contribution_signatures.csv`).
7. `R/signature_r2_diet_vs_covariates.R` — sequential R² (Supplementary Table 5).
8. `R/biomarker_phyto202509.R` — partial correlations with biomarkers.
9. `R/correlations.R` — signature–intake correlations (discovery + WHI).
10. `R/figure2_signature_biomarker.R` — Figure 2 (after biomarker CSV exists).
11. `R/figure3_intake_dominance.R` — Figure 3.
12. `R/figure4_metabolite_biomarker.R` — Figure 4.

WHI-specific steps are documented in the headers of `correlations.R` and related figure scripts.

## Data availability

Individual-level data from NHS, NHSII, HPFS, and WHI are available under controlled access; see the manuscript Data Availability statement.

## Citation

If you use this code, cite the corresponding peer-reviewed publication (add citation when available).

## License

Code is provided for reproducibility of published analyses. Contact the authors for questions about cohort data access.
