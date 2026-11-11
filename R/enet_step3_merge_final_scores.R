#### Step 3 — Final signature scores on full discovery cohort
#### Uses modal metabolite count rule; training scores = LOO; test scores = final repetition predict.
#### Outputs: Prediction_merge_<tag>.RData (object `signature`) and Prediction_results_<tag>.RData (coefs for figures)

source("00_paths.R")
source("enet_signature_config.R")

library(glmnet)

enet_run_step3 <- function(sig_name) {
  cfg <- enet_signatures[[sig_name]]

  imp_path <- file.path(repo, "imp_ztran_met_PR_multi_ep.RData")
  step1 <- step1_path(cfg)
  loo_f <- loo_path(cfg)
  if (!file.exists(imp_path)) stop("Missing ", imp_path, call. = FALSE)
  if (!file.exists(step1)) stop("Missing ", step1, call. = FALSE)
  if (!file.exists(loo_f)) stop("Missing LOO file: ", loo_f, " (run step 2).", call. = FALSE)

  load(imp_path)
  e1 <- new.env()
  load(step1, envir = e1)
  train    <- get("train", envir = e1)
  test     <- get("test", envir = e1)
  results  <- get(cfg$results_obj, envir = e1)
  met_coef <- get(cfg$coef_obj, envir = e1)
  load(loo_f)

  fin <- modal_final_repetition(results, met_coef)
  final_rep <- fin$final_rep
  final_col <- fin$final_col

  test_col <- paste0(cfg$test_prefix, final_rep)
  if (!test_col %in% names(test)) {
    stop("Test column ", test_col, " not found. Re-run step 1.", call. = FALSE)
  }

  if (!"newid" %in% names(AllCohorts_clean) || !"newid" %in% names(train)) {
    stop("Object newid required on AllCohorts_clean and train.", call. = FALSE)
  }

  loo_col <- grep("loo_signature", names(train), value = TRUE)[1]
  if (is.na(loo_col)) loo_col <- ncol(train)

  train_ids <- train$newid
  test_ids  <- test$newid
  train_scores <- setNames(train[[loo_col]], train$newid)
  test_scores  <- setNames(test[[test_col]], test$newid)

  full_score <- rep(NA_real_, nrow(AllCohorts_clean))
  names(full_score) <- AllCohorts_clean$newid
  full_score[names(train_scores)] <- train_scores
  full_score[names(test_scores)]  <- test_scores

  sig_df <- data.frame(
    newid = AllCohorts_clean$newid,
    stringsAsFactors = FALSE
  )
  sig_df[[cfg$sig_var]] <- as.numeric(full_score[sig_df$newid])

  if ("cohort" %in% names(AllCohorts_clean)) {
    sig_df$cohort <- AllCohorts_clean$cohort
  }

  signature <- sig_df

  # Coefficient vectors for Venn / heatmap (non-zero in final repetition)
  is_int <- met_coef$metabolite_names == "intercept"
  keep   <- !is_int & met_coef[[final_col]] != 0
  met_out <- met_coef$metabolite_names[keep]
  coef_out <- met_coef[[final_col]][keep]
  assign(cfg$met_coef_obj[["met"]], met_out, envir = environment())
  assign(cfg$met_coef_obj[["coef"]], coef_out, envir = environment())

  AllCohorts_clean_ztransformed <- AllCohorts_clean

  merge_file   <- merge_path(cfg)
  results_file <- results_path(cfg)

  save(signature, file = merge_file)

  met_nm <- cfg$met_coef_obj[["met"]]
  coef_nm <- cfg$met_coef_obj[["coef"]]
  assign(met_nm, met_out)
  assign(coef_nm, coef_out)
  assign(cfg$results_obj, results)
  assign(cfg$coef_obj, met_coef)

  save(
    AllCohorts_clean_ztransformed,
    list = c(met_nm, coef_nm, cfg$results_obj, cfg$coef_obj),
    file = results_file
  )

  cat(sig_name, ": modal n metabolites =", fin$modal_n,
      ", final repetition =", final_rep, "\n")
  cat("Wrote", merge_file, "\n")
  cat("Wrote", results_file, "\n")
  invisible(signature)
}

enet_run_all_step3 <- function() {
  for (nm in names(enet_signatures)) enet_run_step3(nm)
}
