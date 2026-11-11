library(domir)

#### Relative contribution of food groups to the metabolomic signatures (dominance analysis)
#### For each signature: share of the variance explained by food groups attributable to each food group,
#### with covariates included in every model so that their contribution is not attributed to foods

# Discovery cohort: imputed data + signature scores (same merges as script.R)
endpoints_rdata <- "/udd/nhtes/git/protein-metab-mortality/Endpoints/Endpoints_dataset-multi-ep.RData"
if (file.exists(endpoints_rdata)) {
  load(endpoints_rdata)
} else {
  load("/udd/nhtes/git/protein-metab-mortality/Enet/Plant/Prediction_results_vprot-PR-multi-ep.RData")
  load("/udd/nhtes/git/protein-metab-mortality/Enet/Plant/Prediction_merge_vprot-PR-multi-ep.RData")
  signature_vprot <- signature
  load("/udd/nhtes/git/protein-metab-mortality/Enet/Animal/Prediction_merge_aprot-PR-multi-ep.RData")
  signature_aprot <- signature
  load("/udd/nhtes/git/protein-metab-mortality/Enet/Protein_ratio/Prediction_merge_PR_multi_ep.RData")
  signature_PR <- signature
  AllCohorts_mortality <- AllCohorts_clean_ztransformed
  AllCohorts_mortality <- merge(AllCohorts_mortality, signature_vprot[, c("newid", "vprot_p_sig")], by = "newid", all.x = TRUE)
  AllCohorts_mortality <- merge(AllCohorts_mortality, signature_aprot[, c("newid", "aprot_p_sig")], by = "newid", all.x = TRUE)
  AllCohorts_mortality <- merge(AllCohorts_mortality, signature_PR[, c("newid", "protein_ratio_sig")], by = "newid", all.x = TRUE)
}

# Same food groups (servings/day) as the food–metabolite heatmap
foods <- c(processed_meat = "procm_avg",
           red_meat       = "rmeat_avg",
           poultry        = "poul_avg",
           fish           = "fish_avg",
           egg            = "eggs_avg",
           dairy          = "dairy_avg",
           legumes        = "legu_avg",
           nuts           = "nuts_avg",
           refined_grains = "rgrain_avg",
           whole_grains   = "wgrain_avg")

foods_missing <- setdiff(foods, colnames(AllCohorts_mortality))
if (length(foods_missing) > 0) {
  stop("Food columns missing from AllCohorts_mortality: ",
       paste(foods_missing, collapse = ", "), call. = FALSE)
}

signatures <- c(plant = "vprot_p_sig", animal = "aprot_p_sig", ratio = "protein_ratio_sig")

# "minimal" (main): age, cohort (captures sex) and total energy; signatures were derived unadjusted
# "energy": total energy only
# "full" (sensitivity): same covariates as the food–metabolite models
adjustment <- "minimal"

covariates <- switch(adjustment,
  minimal = c("agemo", "cohort", "calor_avg"),
  energy  = "calor_avg",
  full   = c("fast", "agemo", "white", "multivit", "htn", "hchol", "hbpmed", "cholmed",
             "calor_avg", "cacoall", "cohort", "endpoint", "bmi_avg", "alco_avg",
             "AHEI_mod_avg", "smoke_avg", "th_avg")
)

model_vars <- c(unname(signatures), unname(foods), covariates)
if (any(duplicated(colnames(AllCohorts_mortality)))) {
  stop("Duplicate column names in AllCohorts_mortality — resolve before dominance analysis.",
       call. = FALSE)
}

dat <- AllCohorts_mortality[complete.cases(AllCohorts_mortality[, model_vars]), model_vars, drop = FALSE]
dat <- as.data.frame(dat)
rownames(dat) <- NULL
cat("N =", nrow(dat), "\n")

factor <- intersect(c("fast", "white", "multivit", "htn", "hchol", "hbpmed", "cholmed", "smoke_avg",
                      "cacoall", "cohort", "endpoint"), covariates)
for (i in seq_along(factor)) {
  dat[[factor[i]]] <- as.factor(dat[[factor[i]]])
}

# domir passes only the formula to .fct — always fit with analysis data `dat` from closure
r2_lm <- function(fml) summary(lm(fml, data = dat))$r.squared

food_contribution <- function(signature) {

  # Covariates must appear in .obj and in .all (domir requirement)
  rhs <- paste(c(unname(foods), covariates), collapse = " + ")
  dom <- domir(
    as.formula(paste(signature, "~", rhs)),
    r2_lm,
    .all = as.formula(paste("~", paste(covariates, collapse = " + ")))
  )

  # Direction of association for each food, from the model including all foods and covariates
  full_fit <- lm(as.formula(paste(signature, "~", paste(c(unname(foods), covariates), collapse = " + "))),
                 data = dat)
  direction <- ifelse(coef(full_fit)[foods] > 0, "+", "-")

  data.frame(signature        = names(signatures)[signatures == signature],
             food_group       = names(foods),
             R2_contribution  = round(dom$General_Dominance[foods], 4),
             pct_contribution = round(100 * dom$General_Dominance[foods] / sum(dom$General_Dominance[foods]), 1),
             direction        = direction,
             R2_foods_total   = round(sum(dom$General_Dominance[foods]), 4),
             R2_covariates    = round(dom$Value_All, 4))
}

contribution_table <- do.call(rbind, lapply(signatures, food_contribution))
rownames(contribution_table) <- NULL
print(contribution_table)

out_dir <- "/udd/nhtes/git/protein-metab-mortality/Review"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_file <- switch(adjustment,
  minimal = "food_contribution_signatures.csv",
  energy  = "food_contribution_signatures_energyadj.csv",
  full    = "food_contribution_signatures_fulladj.csv"
)
write.csv(contribution_table, file.path(out_dir, out_file), row.names = FALSE)

# ---- Benchmark: food groups vs the protein intake each signature was built to predict ----
intakes <- c(plant = "vprot_avg_p", animal = "aprot_avg_p", ratio = "protein_ratio")

bench_vars <- c(model_vars, unname(intakes))
bench_dat <- as.data.frame(
  AllCohorts_mortality[complete.cases(AllCohorts_mortality[, bench_vars]), bench_vars, drop = FALSE]
)
for (v in factor) bench_dat[[v]] <- as.factor(bench_dat[[v]])
cat("Benchmark N =", nrow(bench_dat), "\n")

r2_beyond <- function(y, x) {
  r2 <- function(terms) summary(lm(reformulate(terms, y), data = bench_dat))$r.squared
  r2(c(x, covariates)) - r2(covariates)
}

benchmark_table <- do.call(rbind, lapply(names(signatures), function(s) {
  r2_intake <- r2_beyond(signatures[[s]], intakes[[s]])
  r2_foods  <- r2_beyond(signatures[[s]], unname(foods))
  r2_both   <- r2_beyond(signatures[[s]], c(unname(foods), intakes[[s]]))
  data.frame(signature = s,
             R2_intake = round(r2_intake, 4),
             R2_foods  = round(r2_foods, 4),
             R2_both   = round(r2_both, 4),
             # Share of the intake-explained variance that the food groups also explain
             pct_intake_captured_by_foods = round(100 * (r2_intake + r2_foods - r2_both) / r2_intake, 1))
}))
print(benchmark_table)

bench_file <- sub("food_contribution_signatures", "food_vs_intake_benchmark", out_file)
write.csv(benchmark_table, file.path(out_dir, bench_file), row.names = FALSE)
