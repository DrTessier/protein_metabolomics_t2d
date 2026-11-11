#### Individual metabolite Cox models for incident T2D using PPP endpoints
#### Recreates the columns in Endpoints/Results_indmetab_db2.csv
####   Metabolite, HR.mod2, p.value.mod2, p.value.mod2.bon
#### Does NOT overwrite the original CSV (writes *_PPP.csv).
####
#### Model 2 = MV Cox, same lifestyle block as the food–metabolite models (includes alco_avg):
####   Surv(tdb2, db2case) ~ metabolite + fast + white + multivit + hchol + cholmed +
####   htn + hbpmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg + alco_avg +
####   calor_avg + strata(cacoall, cohort, endpoint, agemo)
#### Metabolites are already z-scored (per 1 SD).
####
#### Run: source(file.path(repo, "metabolite_t2d_cox_ppp.R"))

source("00_paths.R")
library(survival)
library(dplyr)

options(device = function(...) grDevices::pdf(file = NULL))
review <- file.path(repo, "Review")
if (!dir.exists(review)) dir.create(review, recursive = TRUE)

old_csv <- file.path(repo, "Endpoints/Results_indmetab_db2.csv")
out_csv <- file.path(repo, "Endpoints/Results_indmetab_db2_PPP.csv")
out_rev <- file.path(review, "Results_indmetab_db2_PPP.csv")
cmp_csv <- file.path(review, "Results_indmetab_db2_PPP_vs_original.csv")

# =============================================================================
# Discovery cohort + PPP T2D (same merge as script.R)
# =============================================================================
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
dat <- merge(AllCohorts3_clean, endpoints_ppp, by = "newid", all.x = TRUE)

factor_vars <- intersect(c("white", "fast", "multivit", "htn", "hchol", "hbpmed", "smoke_avg",
                           "cholmed", "cohort", "cacoall", "endpoint"), names(dat))
for (v in factor_vars) dat[[v]] <- as.factor(dat[[v]])

covar <- c("fast", "white", "multivit", "hchol", "cholmed", "htn", "hbpmed", "ses",
           "smoke_avg", "AHEI_mod_avg", "th_avg", "bmi_avg", "alco_avg", "calor_avg")
miss_cov <- setdiff(c(covar, "tdb2", "db2case", "cacoall", "cohort", "endpoint", "agemo"), names(dat))
if (length(miss_cov)) stop("Missing covariates: ", paste(miss_cov, collapse = ", "), call. = FALSE)

# Metabolite columns: HMDB* on the imputed z-scored table
met_cols <- grep("^HMDB", names(dat), value = TRUE, ignore.case = TRUE)
if (length(met_cols) == 0) {
  stop("No HMDB* metabolite columns found on the discovery table.", call. = FALSE)
}

ok <- complete.cases(dat[, c("tdb2", "db2case", covar, "cacoall", "cohort", "endpoint", "agemo")])
ok <- ok & dat$tdb2 > 0 & !is.na(dat$db2case)
dat_cox <- dat[ok, , drop = FALSE]
n_cox <- nrow(dat_cox)
n_evt <- sum(dat_cox$db2case == 1, na.rm = TRUE)
cat("PPP T2D Cox sample: n =", n_cox, " events =", n_evt, "\n")
cat("Metabolites:", length(met_cols), "\n")

# =============================================================================
# Per-metabolite Cox (model 2)
# =============================================================================
fit_one <- function(met) {
  fml <- as.formula(paste0(
    "Surv(tdb2, db2case) ~ `", met, "` + fast + white + multivit + hchol + cholmed + ",
    "htn + hbpmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg + alco_avg + calor_avg + ",
    "strata(cacoall, cohort, endpoint, agemo)"
  ))
  out <- tryCatch({
    fit <- coxph(fml, data = dat_cox, ties = "efron")
    s <- summary(fit)
    rn <- rownames(s$coefficients)
    i <- which(rn == met | rn == paste0("`", met, "`"))[1]
    if (is.na(i)) i <- 1L
    list(
      HR.mod2 = unname(exp(s$coefficients[i, "coef"])),
      p.value.mod2 = unname(s$coefficients[i, "Pr(>|z|)"]),
      n = n_cox,
      n_events = n_evt,
      ok = TRUE
    )
  }, error = function(e) {
    list(HR.mod2 = NA_real_, p.value.mod2 = NA_real_, n = n_cox, n_events = n_evt, ok = FALSE)
  })
  data.frame(Metabolite = met, HR.mod2 = out$HR.mod2, p.value.mod2 = out$p.value.mod2,
             n = out$n, n_events = out$n_events, converged = out$ok, stringsAsFactors = FALSE)
}

rows <- vector("list", length(met_cols))
for (i in seq_along(met_cols)) {
  if (i %% 50 == 0) cat("  ", i, "/", length(met_cols), "\n", sep = "")
  rows[[i]] <- fit_one(met_cols[i])
}
res <- do.call(rbind, rows)
res$p.value.mod2.bon <- p.adjust(res$p.value.mod2, method = "bonferroni")
res <- res[order(res$p.value.mod2), ]

csv_out <- res[, c("Metabolite", "HR.mod2", "p.value.mod2", "p.value.mod2.bon", "n", "n_events", "converged")]
write.csv(csv_out, out_csv, row.names = FALSE)
write.csv(csv_out, out_rev, row.names = FALSE)
cat("Wrote:\n ", out_csv, "\n ", out_rev, "\n", sep = "")
cat("Converged:", sum(res$converged), "/", nrow(res),
    "  Bonferroni p<0.05:", sum(res$p.value.mod2.bon <= 0.05, na.rm = TRUE),
    "  p<0.05:", sum(res$p.value.mod2 <= 0.05, na.rm = TRUE), "\n")

# =============================================================================
# Optional: compare to original CSV if present (heatmap still uses the original until rebuilt)
# =============================================================================
if (file.exists(old_csv)) {
  old <- read.csv(old_csv, stringsAsFactors = FALSE)
  keep_old <- intersect(c("Metabolite", "HR.mod2", "p.value.mod2", "p.value.mod2.bon"), names(old))
  old <- old[, keep_old, drop = FALSE]
  names(old)[names(old) != "Metabolite"] <- paste0(names(old)[names(old) != "Metabolite"], "_old")
  both <- merge(csv_out, old, by = "Metabolite", all = TRUE)
  if ("HR.mod2_old" %in% names(both) && "HR.mod2" %in% names(both)) {
    both$logHR_diff <- log(both$HR.mod2) - log(both$HR.mod2_old)
    both$star_old <- ifelse(!is.na(both$p.value.mod2.bon_old) & both$p.value.mod2.bon_old <= 0.05, "**",
                            ifelse(!is.na(both$p.value.mod2_old) & both$p.value.mod2_old <= 0.05, "*", ""))
    both$star_ppp <- ifelse(!is.na(both$p.value.mod2.bon) & both$p.value.mod2.bon <= 0.05, "**",
                            ifelse(!is.na(both$p.value.mod2) & both$p.value.mod2 <= 0.05, "*", ""))
    both$star_changed <- both$star_old != both$star_ppp
    write.csv(both, cmp_csv, row.names = FALSE)
    cat("Comparison vs original:\n ", cmp_csv, "\n", sep = "")
    cat("  Metabolites in both:", sum(!is.na(both$HR.mod2) & !is.na(both$HR.mod2_old)), "\n")
    cat("  Star changes (*, **):", sum(both$star_changed, na.rm = TRUE), "\n")
    cat("  Median |log HR difference|:",
        round(median(abs(both$logHR_diff), na.rm = TRUE), 4), "\n")
  }
} else {
  cat("Original", old_csv, "not found — skipped comparison.\n")
}

cat("\nNext: merge Results_indmetab_db2_PPP.csv into the heatmap RData (do not overwrite the original CSV until you have checked the comparison).\n")
