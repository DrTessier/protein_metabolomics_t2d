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
    cor_col       = "vprot_cor",
    results_obj   = "results_enet_vprot",
    coef_obj      = "met_coef_vprot",
    sig_var       = "vprot_p_sig",
    met_coef_obj  = c(met = "met_vprot", coef = "vprot_coef"),
    csv_name      = "signature_plant.csv"
  ),
  animal = list(
    label         = "animal",
    enet_subdir   = "Animal",
    file_tag      = "aprot-PR-multi-ep",
    file_tag_alt  = "aprot-PR-multi-ep",
    intake        = "aprot_avg_p",
    test_prefix   = "aprot",
    cor_col       = "aprot_cor",
    results_obj   = "results_enet_aprot",
    coef_obj      = "met_coef_aprot",
    sig_var       = "aprot_p_sig",
    met_coef_obj  = c(met = "met_aprot", coef = "aprot_coef"),
    csv_name      = "signature_animal.csv"
  ),
  ratio = list(
    label         = "ratio",
    enet_subdir   = "Protein_ratio",
    file_tag      = "PR_multi_ep",
    file_tag_alt  = "PR-multi-ep",
    intake        = "protein_ratio",
    test_prefix   = "pr",
    cor_col       = "pr_cor",
    results_obj   = "results_enet_pr",
    coef_obj      = "met_coef_pr",
    sig_var       = "protein_ratio_sig",
    met_coef_obj  = c(met = "met_pr", coef = "pr_coef"),
    csv_name      = "signature_ratio.csv"
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

modal_final_repetition <- function(results) {
  # First repetition whose metabolite count is the most frequent count.
  tab <- table(results$NoMetabs)
  modal_n <- as.numeric(names(tab)[which.max(tab)])
  final_rep <- which(results$NoMetabs == modal_n)[1]
  list(modal_n = modal_n, final_rep = final_rep,
       coef_col = paste0("inrep", final_rep),
       test_col = NULL)
}

save_with_seed <- function(envir, file) {
  assign(".Random.seed", get(".Random.seed", envir = .GlobalEnv), envir = envir)
  save(list = ls(all.names = TRUE, envir = envir), file = file, envir = envir)
}

restore_seed <- function(envir) {
  if (!exists(".Random.seed", envir = envir, inherits = FALSE)) {
    stop("Saved workspace has no .Random.seed; step 2 would not match step 1.", call. = FALSE)
  }
  assign(".Random.seed", get(".Random.seed", envir = envir), envir = .GlobalEnv)
}
