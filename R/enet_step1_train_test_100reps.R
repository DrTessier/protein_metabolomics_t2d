#### Step 1 — Elastic net: 70/30 split (fixed seed), 100 CV repetitions, lambda.1se, alpha = 0.5
#### Input:  file.path(repo, "imp_ztran_met_PR_multi_ep.RData")  -> AllCohorts_clean
#### Output: Enet/<Plant|Animal|Protein_ratio>/Prediction_train_test_<tag>.RData
####
#### Run one signature:  enet_run_step1("plant")   # "animal" or "ratio"
#### Run all:            source("enet_step1_train_test_100reps.R"); enet_run_all_step1()

source("00_paths.R")
source("enet_signature_config.R")

library(glmnet)
library(dplyr)

enet_run_step1 <- function(sig_name) {
  cfg <- enet_signatures[[sig_name]]
  if (is.null(cfg)) stop("Unknown signature: ", sig_name)

  imp_path <- file.path(repo, "imp_ztran_met_PR_multi_ep.RData")
  if (!file.exists(imp_path)) {
    stop("Missing imputed cohort: ", imp_path,
         "\nRun impute_metabolites_qc.R first (or place the .RData at repo).", call. = FALSE)
  }
  load(imp_path)
  if (!exists("AllCohorts_clean")) stop("AllCohorts_clean not found in ", imp_path, call. = FALSE)

  dat <- AllCohorts_clean
  if (!cfg$intake %in% names(dat)) {
    stop("Intake column ", cfg$intake, " missing in AllCohorts_clean.", call. = FALSE)
  }

  met_cols <- seq_len(N_METABOLITES)
  if (ncol(dat) < max(met_cols)) {
    stop("Expected at least ", N_METABOLITES, " metabolite columns.", call. = FALSE)
  }

  set.seed(ENET_SEED)
  n <- floor(nrow(dat) * TRAIN_FRAC)
  train_index <- sample(seq_len(nrow(dat)), size = n)
  train <- dat[train_index, ]
  test  <- dat[-train_index, ]

  results <- data.frame(repNo = NA_character_, NoMetabs = NA_integer_,
                        cor_test = NA_real_, mse = NA_real_)
  results_obj <- cfg$results_obj
  coef_obj    <- cfg$coef_obj
  assign(results_obj, results, envir = environment())

  met_coef <- data.frame(
    metabolite_names = c("intercept", colnames(train[, met_cols])),
    stringsAsFactors = FALSE
  )

  for (inrep in seq_len(ENET_NREP)) {
    if (inrep %% 10L == 0L) cat(cfg$label, "repetition", inrep, "/", ENET_NREP, "\n")

    Training_CV <- cv.glmnet(
      as.matrix(train[, met_cols]),
      train[[cfg$intake]],
      nfolds = 10, alpha = ENET_ALPHA, family = "gaussian"
    )
    lambda_1se <- Training_CV$lambda.1se

    Training_M <- glmnet(
      as.matrix(train[, met_cols]),
      train[[cfg$intake]],
      family = "gaussian", alpha = ENET_ALPHA
    )
    cmin <- coef(Training_M, s = lambda_1se)
    met_coef[[paste0("inrep", inrep)]] <- as.numeric(cmin)

    pred_test <- as.numeric(
      predict(Training_M, as.matrix(test[, met_cols]), type = "response", s = lambda_1se)[, 1]
    )
    test[[paste0(cfg$test_prefix, inrep)]] <- pred_test

    n_met <- sum(cmin[-1, 1] != 0)
    results[inrep, "repNo"]   <- paste0("inrep", inrep)
    results[inrep, "NoMetabs"] <- n_met
    results[inrep, "cor_test"] <- cor(test[[cfg$intake]], pred_test, use = "complete.obs")
    results[inrep, "mse"]      <- min(Training_CV$cvm)
  }

  assign(results_obj, results)
  assign(coef_obj, met_coef)

  out_dir <- file.path(repo, "Enet", cfg$enet_subdir)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  out_file <- step1_path(cfg)
  save(train, test, list = c(results_obj, coef_obj), file = out_file)
  cat("Wrote", out_file, "\n")
  invisible(list(train = train, test = test, results = get(results_obj), coef = get(coef_obj)))
}

enet_run_all_step1 <- function() {
  for (nm in names(enet_signatures)) enet_run_step1(nm)
}

# When sourced interactively, uncomment:
# enet_run_all_step1()
