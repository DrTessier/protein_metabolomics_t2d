#### Step 3 — choose the final repetition and merge training and test scores
#### Final repetition = first repetition with the most frequent metabolite count.
#### Test scores come from that repetition. Training scores come from the leave-one-out column.

source("00_paths.R")
source("enet_signature_config.R")

library(dplyr)

enet_run_step3 <- function(sig_name) {
  cfg <- enet_signatures[[sig_name]]
  loo_f <- file.path(repo, "Enet", cfg$enet_subdir,
                     paste0("Prediction_getscore_trainset_", cfg$file_tag, ".RData"))
  if (!file.exists(loo_f)) stop("Run step 2 first: ", loo_f, call. = FALSE)

  saved <- new.env()
  load(loo_f, envir = saved)
  train <- saved$train
  test <- saved$test
  results <- saved[[cfg$results_obj]]
  met_coef <- saved[[cfg$coef_obj]]

  fin <- modal_final_repetition(results)
  fin$test_col <- paste(cfg$test_prefix, fin$final_rep, sep = "")
  cat(cfg$label, ": metabolite count =", fin$modal_n, "\n")
  if (!fin$test_col %in% names(test)) {
    stop("Test column ", fin$test_col, " not found.", call. = FALSE)
  }

  predict_score <- test[, c("newid", cfg$intake, "cohort", fin$test_col)]
  colnames(predict_score)[4] <- cfg$sig_var

  met_keep <- met_coef[-1, ]
  met_keep <- met_keep[met_keep[[fin$coef_col]] != 0, ]
  met_vec <- as.character(met_keep$metabolite_names)
  coef_vec <- as.character(met_keep[[fin$coef_col]])

  results_env <- new.env()
  results_env$train <- train
  results_env$test <- test
  results_env[[cfg$results_obj]] <- results
  results_env[[cfg$coef_obj]] <- met_coef
  results_env[[paste0(cfg$test_prefix, "_predict_score")]] <- predict_score
  results_env[[cfg$met_coef_obj[["met"]]]] <- met_vec
  results_env[[cfg$met_coef_obj[["coef"]]]] <- coef_vec
  if (exists("AllCohorts_clean_ztransformed", envir = saved, inherits = FALSE)) {
    results_env$AllCohorts_clean_ztransformed <- saved$AllCohorts_clean_ztransformed
  }
  save(list = ls(results_env), file = results_path(cfg), envir = results_env)

  loo_col <- dim(train)[2]
  predict_score_training <- train[, c("newid", cfg$intake, "cohort")]
  predict_score_training[[cfg$sig_var]] <- train[, loo_col]
  predict_score_training <- predict_score_training[, c("newid", cfg$intake, "cohort", cfg$sig_var)]

  if (!identical(colnames(predict_score), colnames(predict_score_training))) {
    stop("Training and test score columns differ for ", sig_name, call. = FALSE)
  }

  predict_score$set <- "Testing"
  predict_score_training$set <- "Training"
  both <- rbind(predict_score, predict_score_training) %>% arrange(cohort, newid)
  signature <- cbind(both[, c("set", "cohort", "newid", cfg$intake, cfg$sig_var)])

  cat(cfg$label, "training correlation:\n")
  print(cor.test(signature[signature$set == "Training", cfg$sig_var],
                 signature[signature$set == "Training", cfg$intake]))
  cat(cfg$label, "test correlation:\n")
  print(cor.test(signature[signature$set == "Testing", cfg$sig_var],
                 signature[signature$set == "Testing", cfg$intake]))

  merge_env <- new.env()
  merge_env$signature <- signature
  save(list = "signature", file = merge_path(cfg), envir = merge_env)
  write.csv(signature, file.path(repo, "Enet", cfg$enet_subdir, cfg$csv_name),
            row.names = FALSE, na = "")
  cat(cfg$label, "step 3 complete\n")
  invisible(signature)
}

enet_run_all_step3 <- function() {
  for (nm in names(enet_signatures)) enet_run_step3(nm)
}
