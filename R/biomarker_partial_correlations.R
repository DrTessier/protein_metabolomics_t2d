#### Biomarkers: merge discovery signatures + T2D with external lab files, then partial correlations
#### Labs: biomarker.RData object `pool` (NOT biomarker_con.RData).
#### Merge: newid only.
#### Run: source(file.path(repo, "biomarker_phyto202509.R"))

source("00_paths.R")
library(dplyr)
review <- file.path(repo, "Review")
rdata  <- biomarker_rdata

FORCE_REMERGE <- TRUE
PARTIAL_COV <- c("agemo", "cohort", "calor_avg")
sig_vars <- c(plant = "vprot_p_sig", animal = "aprot_p_sig", ratio = "protein_ratio_sig")
merged_rdata <- file.path(review, "AllCohorts_mortality_with_biomarker_phyto202509.RData")

# Available labs from the biomarker `bio` object (do not add others). BMI is already on the cohort table.
bio_labels <- c(
  `C-peptide`  = "adj_lncpep",
  HbA1c        = "adj_lna1c",
  `Total chol` = "adj_lnchol",
  LDL          = "adj_lnldl",
  `TAG/HDL-C`  = "adj_lnratio",
  CRP          = "adj_lncrp",
  Adiponectin  = "adj_lnadipo",
  Leptin       = "adj_lnleptin"
)
biomarker_labels <- c(BMI = "bmi_avg", bio_labels)

# ---- Discovery cohort: 3 signatures + PPP T2D (same logic as script.R) ----
build_discovery_cohort_ppp <- function(repo) {
  load(file.path(repo, "Enet/Plant/Prediction_results_vprot-PR-multi-ep.RData"))
  load(file.path(repo, "Enet/Plant/Prediction_merge_vprot-PR-multi-ep.RData"))
  signature_vprot <- signature
  load(file.path(repo, "Enet/Animal/Prediction_merge_aprot-PR-multi-ep.RData"))
  signature_aprot <- signature
  load(file.path(repo, "Enet/Protein_ratio/Prediction_merge_PR_multi_ep.RData"))
  signature_PR <- signature
  AllCohorts3 <- merge(AllCohorts_clean_ztransformed, signature_vprot[, c("newid", "vprot_p_sig")], by = "newid", all.x = TRUE)
  AllCohorts3 <- merge(AllCohorts3, signature_aprot[, c("newid", "aprot_p_sig")], by = "newid", all.x = TRUE)
  AllCohorts3 <- merge(AllCohorts3, signature_PR[, c("newid", "protein_ratio_sig")], by = "newid", all.x = TRUE)
  t2d_cols <- c("newid", "db2case", "tdb2", "dead", "tdead", "tlong", "deadage_yr",
                "birthday", "age_end", "age_start", "lastq", "followup", "cutoffdead", "exrec", "ses")
  ppp_paths <- file.path(repo, c("Endpoints/nhs1coxPPP.csv", "Endpoints/nhs2coxPPP.csv", "Endpoints/hpcoxPPP.csv"))
  if (any(!file.exists(ppp_paths))) stop("Missing PPP CSV under Endpoints/.", call. = FALSE)
  nhs1 <- read.csv(ppp_paths[1], header = TRUE)
  nhs2 <- read.csv(ppp_paths[2], header = TRUE)
  hpfs <- read.csv(ppp_paths[3], header = TRUE)
  ep <- rbind(
    nhs1[!duplicated(nhs1$newid), t2d_cols, drop = FALSE],
    nhs2[!duplicated(nhs2$newid), t2d_cols, drop = FALSE],
    hpfs[!duplicated(hpfs$newid), t2d_cols, drop = FALSE]
  )
  endpoints_ppp <- ep[ep$tdb2 > 0, , drop = FALSE]
  AllCohorts3_clean <- AllCohorts3[, !(names(AllCohorts3) %in% setdiff(t2d_cols, "newid")), drop = FALSE]
  merge(AllCohorts3_clean, endpoints_ppp, by = "newid", all.x = TRUE)
}

load_cohort <- function(repo) {
  endpoints_rdata <- file.path(repo, "Endpoints/Endpoints_dataset-multi-ep.RData")
  if (file.exists(endpoints_rdata)) {
    env <- new.env()
    load(endpoints_rdata, envir = env)
    cat("Cohort:", endpoints_rdata, " rows:", nrow(env$AllCohorts_mortality), "\n")
    return(list(dat = env$AllCohorts_mortality, source = endpoints_rdata))
  }
  cat("Cohort: Enet + PPP CSVs\n")
  list(dat = build_discovery_cohort_ppp(repo), source = "Enet + PPP")
}

partial_cor_test <- function(x, y, cov_df) {
  ok <- complete.cases(data.frame(x = x, y = y, cov_df))
  n <- sum(ok)
  if (n < 10) return(list(n = n, r = NA_real_, p = NA_real_))
  cov_ok <- cov_df[ok, , drop = FALSE]
  rx <- residuals(lm(x[ok] ~ ., data = cov_ok))
  ry <- residuals(lm(y[ok] ~ ., data = cov_ok))
  ct <- stats::cor.test(rx, ry)
  list(n = n, r = unname(ct$estimate), p = ct$p.value)
}

if (!dir.exists(review)) dir.create(review, recursive = TRUE)

# =============================================================================
# 1. Load pool
# =============================================================================
if (!file.exists(rdata)) stop("Missing ", rdata, call. = FALSE)
bio_env <- new.env()
load(rdata, envir = bio_env)
pool <- bio_env$pool
bio_lab_vars <- if (exists("bio", envir = bio_env)) unname(bio_env$bio) else unname(bio_labels)
cat("\n=== 1. pool ===\n")
cat("rows:", nrow(pool), " unique newid:", length(unique(pool$newid)), "\n")
cat("bio:", paste(bio_lab_vars, collapse = ", "), "\n")

# =============================================================================
# 2. Merge on newid only + spot-check 2 shared variables
# =============================================================================
if (!FORCE_REMERGE && file.exists(merged_rdata)) {
  menv <- new.env()
  load(merged_rdata, envir = menv)
  dat_bio <- menv$AllCohorts_mortality_with_biomarker
  cat("\n=== 2. Using saved merge ===\n")
} else {
  cat("\n=== 2. Merge on newid ===\n")
  cohort_info <- load_cohort(repo)
  dat <- as.data.frame(cohort_info$dat)
  pool_m <- as.data.frame(pool)
  if (anyDuplicated(pool_m$newid)) {
    cat("pool duplicate newid — keeping first row per person\n")
    pool_m <- pool_m[!duplicated(pool_m$newid), , drop = FALSE]
  }

  # Spot-check copies (not merge keys): cohort and endpoint if present on both sides
  if ("cohort" %in% names(pool_m)) pool_m$cohort_pool <- as.character(pool_m$cohort)
  if ("endpoint" %in% names(pool_m)) pool_m$endpoint_pool <- as.character(pool_m$endpoint)

  overlap <- setdiff(intersect(names(pool_m), names(dat)), c("newid", "cohort_pool", "endpoint_pool"))
  if (length(overlap)) {
    cat("Renaming overlapping pool columns:", paste(overlap, collapse = ", "), "\n")
    names(pool_m)[names(pool_m) %in% overlap] <- paste0(names(pool_m)[names(pool_m) %in% overlap], "_pool")
  }

  n_id <- length(intersect(dat$newid, pool_m$newid))
  cat("newid overlap:", n_id, "\n")
  n_before <- nrow(dat)
  dat_bio <- merge(dat, pool_m, by = "newid", all.x = TRUE, sort = FALSE)
  cat("Rows:", n_before, "->", nrow(dat_bio), "\n")

  # Spot check: cohort / endpoint among people with a lab
  lab1 <- bio_lab_vars[bio_lab_vars %in% names(dat_bio)][1]
  has_lab <- if (!is.na(lab1)) !is.na(dat_bio[[lab1]]) else !is.na(dat_bio$cohort_pool)
  cat("N with lab (", lab1, "): ", sum(has_lab), "\n", sep = "")
  if ("cohort" %in% names(dat_bio) && "cohort_pool" %in% names(dat_bio)) {
    sub <- dat_bio[has_lab, ]
    cat("Spot-check cohort (among those with lab):\n")
    print(table(cohort = as.character(sub$cohort), pool = as.character(sub$cohort_pool), useNA = "ifany"))
  }
  if ("endpoint" %in% names(dat_bio) && "endpoint_pool" %in% names(dat_bio)) {
    sub <- dat_bio[has_lab, ]
    n_end <- sum(as.character(sub$endpoint) == as.character(sub$endpoint_pool), na.rm = TRUE)
    cat("Spot-check endpoint equal (among those with lab):", n_end, "/", nrow(sub), "\n")
  }
  set.seed(1)
  pick <- intersect(c("newid", "cohort", "cohort_pool", "endpoint", "endpoint_pool", lab1), names(dat_bio))
  idx <- c(which(has_lab)[seq_len(min(20, sum(has_lab)))],
           which(!has_lab)[seq_len(min(10, sum(!has_lab)))])
  write.csv(dat_bio[idx, pick, drop = FALSE],
            file.path(review, "biomarker_merge_spotcheck_sample_phyto202509.csv"),
            row.names = FALSE)

  merge_save_env <- new.env(parent = emptyenv())
  list2env(list(
    AllCohorts_mortality_with_biomarker = dat_bio,
    bio_lab_vars = bio_lab_vars,
    merge_keys = "newid",
    cohort_source = cohort_info$source
  ), envir = merge_save_env)
  save(list = ls(envir = merge_save_env), envir = merge_save_env, file = merged_rdata)
  cat("Saved:", merged_rdata, "\n")
}

# =============================================================================
# 3. N with signatures + available biomarkers
# =============================================================================
miss_sig <- setdiff(unname(sig_vars), names(dat_bio))
if (length(miss_sig)) stop("Missing signatures: ", paste(miss_sig, collapse = ", "))
sig_ok <- complete.cases(dat_bio[, unname(sig_vars), drop = FALSE])
cat("\n=== 3. Available biomarkers (bio vector + BMI) ===\n")
cat("All 3 signatures non-missing:", sum(sig_ok), "\n")
inv <- data.frame(
  label = names(biomarker_labels),
  column = unname(biomarker_labels),
  n_with_signature = NA_integer_,
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(inv))) {
  col <- inv$column[i]
  inv$n_with_signature[i] <- if (col %in% names(dat_bio)) sum(sig_ok & !is.na(dat_bio[[col]])) else 0L
}
print(inv, row.names = FALSE)
write.csv(inv, file.path(review, "biomarker_n_phyto202509.csv"), row.names = FALSE)

# =============================================================================
# 4. Partial correlations: 3 signatures vs available biomarkers
# =============================================================================
cat("\n=== 4. Partial correlations ===\n")
cat("Adjustment:", paste(PARTIAL_COV, collapse = ", "), "\n")
cov_cols <- intersect(PARTIAL_COV, names(dat_bio))
if (length(cov_cols) < length(PARTIAL_COV)) {
  stop("Missing covariates: ", paste(setdiff(PARTIAL_COV, cov_cols), collapse = ", "), call. = FALSE)
}
markers_use <- biomarker_labels[biomarker_labels %in% names(dat_bio)]
rows_pc <- list()
for (s in names(sig_vars)) {
  for (m_label in names(markers_use)) {
    cov_df <- dat_bio[, cov_cols, drop = FALSE]
    if ("cohort" %in% cov_cols) cov_df$cohort <- as.factor(cov_df$cohort)
    pc <- partial_cor_test(dat_bio[[sig_vars[[s]]]], dat_bio[[markers_use[[m_label]]]], cov_df)
    rows_pc[[length(rows_pc) + 1]] <- data.frame(
      signature = s,
      biomarker_label = m_label,
      biomarker_column = unname(markers_use[m_label]),
      n = pc$n,
      partial_r = pc$r,
      p_value = pc$p,
      stringsAsFactors = FALSE
    )
  }
}
partial_tab <- do.call(rbind, rows_pc)
partial_tab$q_fdr <- p.adjust(partial_tab$p_value, method = "BH")
partial_tab$star <- ifelse(is.na(partial_tab$p_value), "", ifelse(partial_tab$p_value <= 0.05, "*", ""))
fmt_p <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "<0.001", format(round(p, 3), nsmall = 3)))
}
partial_tab$p_fmt <- fmt_p(partial_tab$p_value)
partial_tab$q_fmt <- fmt_p(partial_tab$q_fdr)
print(partial_tab[, c("signature", "biomarker_label", "n", "partial_r", "p_fmt", "q_fmt", "star")],
      row.names = FALSE)
write.csv(partial_tab, file.path(review, "biomarker_partial_correlations_phyto202509.csv"), row.names = FALSE)

sig_lab <- c(plant = "Plant protein signature", animal = "Animal protein signature",
             ratio = "Plant-to-animal ratio signature")
supp <- data.frame(
  Signature = unname(sig_lab[partial_tab$signature]),
  Biomarker = partial_tab$biomarker_label,
  n = partial_tab$n,
  Partial_r = partial_tab$partial_r,
  P = partial_tab$p_value,
  Q_FDR = partial_tab$q_fdr,
  stringsAsFactors = FALSE
)
out_supp <- file.path(review, "SupplementaryTable_biomarker_partial_correlations.csv")
op <- options(digits = 15, scipen = 0)
on.exit(options(op), add = TRUE)
write.csv(supp, out_supp, row.names = FALSE)
cat("Saved supplement:", out_supp, "\n")
cat("FDR: Benjamini-Hochberg across", sum(!is.na(partial_tab$p_value)), "tests.\n")

if (requireNamespace("ggplot2", quietly = TRUE) && any(!is.na(partial_tab$partial_r))) {
  options(device = function(...) grDevices::pdf(file = NULL))
  library(ggplot2)
  sig_order <- c("plant", "animal", "ratio")
  sig_plot_labels <- c(
    plant = "Plant protein\nsignature",
    animal = "Animal protein\nsignature",
    ratio = "Plant-to-animal\nratio signature"
  )
  heat_colours <- c("#0046FF", "#5876FF", "#7F97FF", "#A9B9FF", "#D4DBFF", "#FFFFFF",
                    "#FFE9F3", "#FFB0D6", "#FF7FBC", "#FF5AA8", "#FF2F92")
  plot_df <- partial_tab %>%
    filter(!is.na(partial_r)) %>%
    mutate(
      signature = factor(signature, levels = sig_order),
      biomarker_label = factor(biomarker_label, levels = rev(names(markers_use)))
    )
  plot_stars <- plot_df[nzchar(plot_df$star), , drop = FALSE]
  lim <- min(max(abs(plot_df$partial_r), na.rm = TRUE), 1)
  if (!is.finite(lim) || lim == 0) lim <- 0.5
  p_bio <- ggplot(plot_df, aes(x = signature, y = biomarker_label, fill = partial_r)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(data = plot_stars, aes(label = star), size = 8 / .pt, vjust = 0.75) +
    scale_fill_gradientn(colours = heat_colours, limits = c(-lim, lim), oob = scales::squish, name = "Partial r") +
    scale_x_discrete(labels = sig_plot_labels[sig_order], expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 7) +
    theme(axis.text.x = element_text(size = 7, lineheight = 0.95),
          axis.text.y = element_text(size = 7),
          legend.position = "bottom",
          legend.key.height = grid::unit(2, "mm"))
  out_png <- file.path(review, "Figure_biomarker_partial_cor_phyto202509.png")
  out_pdf <- file.path(review, "Figure_biomarker_partial_cor_phyto202509.pdf")
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(out_png, p_bio, width = 120, height = 90, units = "mm", dpi = 300, device = ragg::agg_png)
  } else {
    ggsave(out_png, p_bio, width = 120, height = 90, units = "mm", dpi = 300, type = "cairo")
  }
  if (capabilities("cairo")) {
    ggsave(out_pdf, p_bio, width = 120, height = 90, units = "mm", device = grDevices::cairo_pdf)
  } else {
    ggsave(out_pdf, p_bio, width = 120, height = 90, units = "mm")
  }
  cat("Saved:", out_png, "\n", out_pdf, "\n")
}
