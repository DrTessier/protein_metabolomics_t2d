#### Sequential R2: diet, then other factors, for each signature and its self-reported intake.
#### Same order as the comparison table: age/energy/cohort, food groups, BMI,
#### lifestyle, then clinical diagnoses and medications.
#### Food groups enter before the other factors, so shared variance is attributed to diet.
#### A second table uses the protein intake itself as the diet step for the signature only.
####
#### Food variables are servings/day. Replace `foods` with protein grams when those
#### variables exist; the self-reported column is then the fair reconstruction of intake.
####
#### Run: source(file.path(repo, "signature_r2_diet_vs_covariates.R"))
source("00_paths.R")
endpoints_rdata <- file.path(repo, "Endpoints/Endpoints_dataset-multi-ep.RData")
if (file.exists(endpoints_rdata)) {
  load(endpoints_rdata)
} else {
  stop("Missing ", endpoints_rdata, call. = FALSE)
}

signatures <- c(plant = "vprot_p_sig", animal = "aprot_p_sig", ratio = "protein_ratio_sig")
intakes   <- c(plant = "vprot_avg_p", animal = "aprot_avg_p", ratio = "protein_ratio")
foods <- c("procm_avg", "rmeat_avg", "poul_avg", "fish_avg", "eggs_avg", "dairy_avg",
           "legu_avg", "nuts_avg", "rgrain_avg", "wgrain_avg")

base       <- c("agemo", "cohort", "calor_avg")
bmi        <- "bmi_avg"
lifestyle  <- c("smoke_avg", "alco_avg", "th_avg", "fast", "white", "multivit", "ses")
clinical   <- c("htn", "hbpmed", "hchol", "cholmed")

need <- unique(c(unname(signatures), unname(intakes), foods, base, bmi, lifestyle, clinical))
miss <- setdiff(need, colnames(AllCohorts_mortality))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)

dat <- as.data.frame(AllCohorts_mortality[complete.cases(AllCohorts_mortality[, need]), need])
dat$cohort <- as.factor(dat$cohort)
cat("N =", nrow(dat), "\n")

r2_of <- function(y, terms) {
  summary(lm(reformulate(terms, response = y), data = dat))$r.squared
}

step_terms <- list(
  "Age, energy, cohort" = base,
  "+ food groups" = c(base, foods),
  "+ BMI" = c(base, foods, bmi),
  "+ smoking, alcohol, activity, fasting, race, multivitamins, SES" = c(base, foods, bmi, lifestyle),
  "+ hypertension, hypercholesterolemia, and medications" = c(base, foods, bmi, lifestyle, clinical)
)

sequential_table <- function(y, label) {
  r2 <- vapply(step_terms, function(terms) r2_of(y, terms), numeric(1))
  increase <- c(NA_real_, diff(r2))
  data.frame(
    outcome = label,
    step = names(step_terms),
    R2 = round(r2, 3),
    increase = round(increase, 3),
    stringsAsFactors = FALSE
  )
}

food_rows <- lapply(names(signatures), function(s) {
  rbind(
    sequential_table(signatures[[s]], paste0(s, " signature")),
    sequential_table(intakes[[s]], paste0(s, " self-reported intake"))
  )
})
food_tab <- do.call(rbind, food_rows)
rownames(food_tab) <- NULL

# Signature only: diet step is the protein intake the score was derived to track,
# then the same covariate steps. Self-reported intake is not an outcome here.
intake_steps <- list(
  "Age, energy, cohort" = base,
  "+ protein intake" = c(base, "INTAKE"),
  "+ BMI" = c(base, "INTAKE", bmi),
  "+ smoking, alcohol, activity, fasting, race, multivitamins, SES" = c(base, "INTAKE", bmi, lifestyle),
  "+ hypertension, hypercholesterolemia, and medications" = c(base, "INTAKE", bmi, lifestyle, clinical)
)
intake_rows <- lapply(names(signatures), function(s) {
  terms_list <- lapply(intake_steps, function(terms) gsub("^INTAKE$", intakes[[s]], terms))
  r2 <- vapply(terms_list, function(terms) r2_of(signatures[[s]], terms), numeric(1))
  increase <- c(NA_real_, diff(r2))
  data.frame(
    outcome = paste0(s, " signature"),
    step = names(intake_steps),
    R2 = round(r2, 3),
    increase = round(increase, 3),
    stringsAsFactors = FALSE
  )
})
intake_tab <- do.call(rbind, intake_rows)
rownames(intake_tab) <- NULL

cat("\n=== Food groups, then other factors (signature and self-reported intake) ===\n")
print(food_tab, row.names = FALSE)
cat("\n=== Protein intake, then other factors (signature only) ===\n")
print(intake_tab, row.names = FALSE)
out_dir <- file.path(repo, "Review")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(food_tab, file.path(out_dir, "signature_r2_foods_then_covariates.csv"), row.names = FALSE)
write.csv(intake_tab, file.path(out_dir, "signature_r2_intake_then_covariates.csv"), row.names = FALSE)
cat("\nWrote Review/signature_r2_foods_then_covariates.csv\n")
cat("Wrote Review/signature_r2_intake_then_covariates.csv\n")
