#### Step 2 — leave-one-out scores for the training set
#### Continues the random-number stream saved at the end of step 1.
#### The score is added as the last column of the training data.

source("00_paths.R")
source("enet_signature_config.R")

library(glmnet)

enet_run_step2 <- function(sig_name) {
  cfg <- enet_signatures[[sig_name]]
  step1 <- step1_path(cfg)
  if (!file.exists(step1)) stop("Run step 1 first: ", step1, call. = FALSE)

  saved <- new.env()
  load(step1, envir = saved)
  restore_seed(saved)
  train <- saved$train

  for (i in 1:dim(train)[1]) {
    if (i %% 50 == 0) cat(cfg$label, "LOO", i, "/", nrow(train), "\n")

    Training <- cv.glmnet(as.matrix(train[-i, c(1:243)]), train[-i, cfg$intake],
                          nfolds = 10, alpha = ENET_ALPHA, family = "gaussian")
    lambda_min_10F <- Training$lambda.1se

    if (i == 1) {
      train[i, dim(train)[2] + 1] <- predict(
        Training, as.matrix(train[i, c(1:243)]), type = "response", s = lambda_min_10F
      )
    } else if (i > 1) {
      train[i, dim(train)[2]] <- predict(
        Training, as.matrix(train[i, c(1:243)]), type = "response", s = lambda_min_10F
      )
    }
  }

  saved$train <- train
  out_file <- file.path(repo, "Enet", cfg$enet_subdir,
                        paste0("Prediction_getscore_trainset_", cfg$file_tag, ".RData"))
  save_with_seed(saved, out_file)
  cat(cfg$label, "step 2 complete\n")
  invisible(out_file)
}

enet_run_all_step2 <- function() {
  for (nm in names(enet_signatures)) enet_run_step2(nm)
}
