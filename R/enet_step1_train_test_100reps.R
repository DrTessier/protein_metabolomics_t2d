#### Step 1 — 70/30 split, set.seed(1234), 100 repetitions, 10-fold CV, alpha 0.5, lambda.1se
#### The random-number state is saved with the workspace so step 2 continues from it.

source("00_paths.R")
source("enet_signature_config.R")

library(glmnet)

enet_run_step1 <- function(sig_name) {
  cfg <- enet_signatures[[sig_name]]
  if (is.null(cfg)) stop("Unknown signature: ", sig_name)

  imp_path <- file.path(repo, "imp_ztran_met_PR_multi_ep.RData")
  if (!file.exists(imp_path)) {
    stop("Missing imputed cohort: ", imp_path, call. = FALSE)
  }
  imp <- new.env()
  load(imp_path, envir = imp)
  if (!exists("AllCohorts_clean", envir = imp)) {
    stop("AllCohorts_clean not found in ", imp_path, call. = FALSE)
  }

  AllCohorts_clean_ztransformed <- imp$AllCohorts_clean

  set.seed(ENET_SEED)
  n <- floor(nrow(AllCohorts_clean_ztransformed) * TRAIN_FRAC)
  train_index <- sample(seq_len(nrow(AllCohorts_clean_ztransformed)), size = n)
  train <- AllCohorts_clean_ztransformed[train_index, ]
  test <- AllCohorts_clean_ztransformed[-train_index, ]

  results <- data.frame(repNo = NA, NoMetabs = NA, cor = NA, mse = NA)
  names(results)[3] <- cfg$cor_col
  met_coef <- data.frame(metabolite_names = c("intercept", colnames(train[, c(1:243)])))

  for (inrep in 1:100) {
    if (inrep %% 10 == 0) cat(cfg$label, "repetition", inrep, "\n")

    Training_CV <- cv.glmnet(as.matrix(train[, c(1:243)]), train[, cfg$intake],
                             nfolds = 10, alpha = ENET_ALPHA, family = "gaussian")
    lambda_min_10F <- Training_CV$lambda.1se

    Training_M <- glmnet(as.matrix(train[, c(1:243)]), train[, cfg$intake],
                         family = "gaussian", alpha = ENET_ALPHA)
    cmin <- coef(Training_M, s = lambda_min_10F)
    metabsmin <- data.frame(cmin[which(cmin[, 1] != 0), ])
    names(metabsmin) <- "Test.min"
    met_coef[, dim(met_coef)[2] + 1] <- cmin[, 1]
    names(met_coef)[dim(met_coef)[2]] <- paste("inrep", inrep, sep = "")

    test[, dim(test)[2] + 1] <- as.numeric(predict(
      Training_M, as.matrix(test[, c(1:243)]), type = "response", s = lambda_min_10F
    )[, 1])
    names(test)[dim(test)[2]] <- paste(cfg$test_prefix, inrep, sep = "")

    results[inrep, "repNo"] <- paste("inrep", inrep, sep = "")
    results[inrep, "NoMetabs"] <- dim(metabsmin)[1] - 1
    results[inrep, cfg$cor_col] <- cor(test[, cfg$intake], test[, dim(test)[2]])
    results[inrep, "mse"] <- min(Training_CV$cvm)
  }

  out_dir <- file.path(repo, "Enet", cfg$enet_subdir)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_file <- step1_path(cfg)

  saved <- new.env()
  saved$AllCohorts_clean_ztransformed <- AllCohorts_clean_ztransformed
  saved$train <- train
  saved$test <- test
  saved[[cfg$results_obj]] <- results
  saved[[cfg$coef_obj]] <- met_coef
  save_with_seed(saved, out_file)
  cat(cfg$label, "step 1 complete\n")
  invisible(out_file)
}

enet_run_all_step1 <- function() {
  for (nm in names(enet_signatures)) enet_run_step1(nm)
}
