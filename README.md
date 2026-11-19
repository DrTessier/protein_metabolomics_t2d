# Plasma metabolomic signatures of plant and animal protein and type 2 diabetes

Analysis code for elastic-net metabolomic signatures of plant protein, animal protein, and the plant-to-animal protein ratio in the Nurses' Health Study (NHS), NHSII, and Health Professionals Follow-up Study (HPFS), with replication in the Women's Health Initiative (WHI), and prospective associations with incident type 2 diabetes.

This repository contains scripts for the analyses reported in the manuscript, including elastic-net signature derivation, Cox models, metabolome-wide associations, biomarker partial correlations, dominance analysis, sequential R², figures, and supplementary tables.

## Requirements

- R (≥ 4.x) with packages used in the scripts (e.g. `survival`, `dplyr`, `tidyverse`, `glmnet`, `domir`, `ggplot2`, `patchwork`, `tableone`, `VennDiagram`, as needed per script).
- Cohort data files, which are not included in this repository. Set their locations in `R/00_paths.R`.

## Path configuration

Edit `R/00_paths.R` before running any script. Replace the placeholders with your local paths:

```r
repo <- "[path]"           # project root (Enet/, Endpoints/, Review/, Heatmap/, …)
biomarker_rdata <- "[path]" # full path to biomarker.RData (labs object `pool`)
met_names_csv <- "[path]"   # full path to metabolite name lookup CSV
```

Each analysis script begins with `source("00_paths.R")` and builds file locations with `file.path(repo, "…")`. Signature scores and imputed metabolomics data are produced upstream (elastic-net outputs and endpoint tables); this repository assumes those files exist under `repo`.

## Elastic net (signature derivation)

From `R/` (after editing `00_paths.R`):

0. `impute_metabolites_qc.R` — optional if `imp_ztran_met_PR_multi_ep.RData` already exists under `repo`.
1. `enet_run_all.R` — or run steps separately:
   - `enet_step1_train_test_100reps.R` → `enet_run_all_step1()`
   - `enet_step2_leave_one_out_training.R` → `enet_run_all_step2()` (slow)
   - `enet_step3_merge_final_scores.R` → `enet_run_all_step3()`

The three signatures use the same steps. Step 1 saves the random-number state and step 2 continues from it. The final repetition is the first of the 100 runs with the most frequent metabolite count.

Outputs are the signature workspace and `signature_plant.csv`, `signature_animal.csv`, and `signature_ratio.csv`.

## Suggested run order (discovery)

1. `R/venn_figure1.R` — signature metabolites, Venn diagram, heatmap inputs (uses saved elastic-net results).
2. `R/table1.R` — Table 1.
3. `R/metabolite_t2d_cox_ppp.R` — per-metabolite Cox models for T2D (supplementary heatmap).
4. `R/heatmap_full_ppp_export.R` — full signature × protein × food × T2D table.
5. `R/script.R` — signature and self-reported protein Cox models (Table 2, supplementary self-report, controls-only sensitivity).
6. `R/food_contribution_signatures.R` — food-group dominance (`Review/food_contribution_signatures.csv`).
7. `R/signature_r2_diet_vs_covariates.R` — sequential R² (Supplementary Table 5).
8. `R/biomarker.R` — partial correlations with biomarkers.
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
