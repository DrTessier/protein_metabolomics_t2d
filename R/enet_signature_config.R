#### Configuration for plant, animal, and ratio elastic-net signatures (discovery cohort).
#### Sourced by enet_step1–3 scripts. Edit only if your saved file names differ.

enet_signatures <- list(
  plant = list(
    label         = "plant",
    enet_subdir   = "Plant",
    file_tag      = "vprot-PR-multi-ep",
    file_tag_alt  = "vprot-multi-ep", # LOO file name on some runs (getscore)
    intake        = "vprot_avg_p",
    test_prefix   = "vprot",
    results_obj   = "results_enet_vprot",
    coef_obj      = "met_coef_vprot",
    sig_var       = "vprot_p_sig",
    met_coef_obj  = c(met = "met_vprot", coef = "vprot_coef")
  ),
  animal = list(
    label         = "animal",
    enet_subdir   = "Animal",
    file_tag      = "aprot-PR-multi-ep",
    file_tag_alt  = "aprot-PR-multi-ep",
    intake        = "aprot_avg_p",
    test_prefix   = "aprot",
    results_obj   = "results_enet_aprot",
    coef_obj      = "met_coef_aprot",
    sig_var       = "aprot_p_sig",
    met_coef_obj  = c(met = "met_aprot", coef = "aprot_coef")
  ),
  ratio = list(
    label         = "ratio",
    enet_subdir   = "Protein_ratio",
    file_tag      = "PR_multi_ep",
    file_tag_alt  = "PR-multi-ep",
    intake        = "protein_ratio",
    test_prefix   = "pr",
    results_obj   = "results_enet_pr",
    coef_obj      = "met_coef_pr",
    sig_var       = "protein_ratio_sig",
    met_coef_obj  = c(met = "met_pr", coef = "pr_coef")
  )
)

N_METABOLITES <- 243L
ENET_ALPHA     <- 0.5
ENET_SEED      <- 1234L
ENET_NREP      <- 100L
TRAIN_FRAC     <- 0.7

step1_path <- function(cfg) {
  file.path(repo, "Enet", cfg$enet_subdir,
            paste0("Prediction_train_test_", cfg$file_tag, ".RData"))
}

loo_path <- function(cfg) {
  p1 <- file.path(repo, "Enet", cfg$enet_subdir,
                  paste0("Prediction_getscore_trainset_", cfg$file_tag, ".RData"))
  if (file.exists(p1)) return(p1)
  file.path(repo, "Enet", cfg$enet_subdir,
            paste0("Prediction_getscore_trainset_", cfg$file_tag_alt, ".RData"))
}

merge_path <- function(cfg) {
  file.path(repo, "Enet", cfg$enet_subdir,
            paste0("Prediction_merge_", cfg$file_tag, ".RData"))
}

results_path <- function(cfg) {
  file.path(repo, "Enet", cfg$enet_subdir,
            paste0("Prediction_results_", cfg$file_tag, ".RData"))
}

modal_final_repetition <- function(results, coef_df) {
  modal_n   <- as.numeric(names(which.max(table(results$NoMetabs))))
  final_rep <- which(results$NoMetabs == modal_n)[1L]
  rep_cols  <- grep("^inrep", names(coef_df), value = TRUE)
  list(modal_n = modal_n, final_rep = final_rep, final_col = rep_cols[final_rep])
}
