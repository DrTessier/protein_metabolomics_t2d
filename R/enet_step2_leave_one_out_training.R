#### Step 2 — Leave-one-out signature scores on the training set (70% partition from step 1)
#### Input:  Prediction_train_test_<tag>.RData
#### Output: Prediction_getscore_trainset_<tag>.RData  (train with LOO score in the last column)
####
#### Warning: slow (~ nrow(train) × cost of cv.glmnet per row).

source("00_paths.R")
source("enet_signature_config.R")

library(glmnet)

enet_run_step2 <- function(sig_name) {
  cfg <- enet_signatures[[sig_name]]
  step1 <- step1_path(cfg)
  if (!file.exists(step1)) stop("Run step 1 first: ", step1, call. = FALSE)
  load(step1)

  met_cols <- seq_len(N_METABOLITES)
  n_train  <- nrow(train)
  loo_col  <- "loo_signature_score"

  train[[loo_col]] <- NA_real_

  for (i in seq_len(n_train)) {
    if (i %% 50L == 0L) cat(cfg$label, "LOO", i, "/", n_train, "\n")
    idx <- setdiff(seq_len(n_train), i)
    Training <- cv.glmnet(
      as.matrix(train[idx, met_cols]),
      train[idx, cfg$intake],
      nfolds = 10, alpha = ENET_ALPHA, family = "gaussian"
    )
    lambda_1se <- Training$lambda.1se
    train[i, loo_col] <- as.numeric(
      predict(Training, as.matrix(train[i, met_cols, drop = FALSE]),
              type = "response", s = lambda_1se)
    )
  }

  out_file <- file.path(repo, "Enet", cfg$enet_subdir,
                        paste0("Prediction_getscore_trainset_", cfg$file_tag, ".RData"))
  if (!file.exists(dirname(out_file))) dir.create(dirname(out_file), recursive = TRUE)
  save(train, file = out_file)
  cat("Wrote", out_file, "\n")
  invisible(train)
}

enet_run_all_step2 <- function() {
  for (nm in names(enet_signatures)) enet_run_step2(nm)
}
