source("00_paths.R")
library(tidyverse)
library(dplyr)
library(survival)


# Mediators: https://www.sciencedirect.com/science/article/abs/pii/S0002870321000739
# for direct effect estimation: adjust
# for total effect estimation: do not adjust for mediators


# Get covariates
load(file.path(repo, "Enet/Plant/Prediction_results_vprot-PR-multi-ep.RData"))

# Get metabolomic signature and
load(file.path(repo, "Enet/Plant/Prediction_merge_vprot-PR-multi-ep.RData"))
signature_vprot <- signature
load(file.path(repo, "Enet/Animal/Prediction_merge_aprot-PR-multi-ep.RData"))
signature_aprot <- signature
load(file.path(repo, "Enet/Protein_ratio/Prediction_merge_PR_multi_ep.RData"))
signature_PR <- signature
gdata::keep(AllCohorts_clean_ztransformed, signature_vprot, signature_aprot, signature_PR, sure=T)


# AllCohorts_clean_ztransformed2 <- AllCohorts_clean_ztransformed[,c(1, 246:307)]

AllCohorts <- merge(AllCohorts_clean_ztransformed, signature_vprot[,c("newid", "vprot_p_sig")], by.x="newid", by.y="newid", all.x=T )

AllCohorts2 <- merge(AllCohorts, signature_aprot[,c("newid", "aprot_p_sig")], by.x="newid", by.y="newid", all.x=T )

AllCohorts3 <- merge(AllCohorts2, signature_PR[,c("newid","protein_ratio_sig")], by.x="newid", by.y="newid", all.x=T )

##########
##### using files from healthy lifestyle metab project as I have NA values for time to event... not sure why
##########

# Columns needed for T2D analysis
t2d_cols <- c("newid",
              "db2case", "tdb2",          # T2D outcome
              "dead", "tdead",             # all-cause mortality (censoring/longevity)
              "tlong",                     # time for longevity analysis
              "deadage_yr",                # age at death
              "birthday", "age_end", "age_start",  # age variables
              "lastq", "followup", "cutoffdead",    # follow-up timing
              "exrec",                     # exclusion flag
              "ses")                       # socioeconomic status

hpfs <- read.csv(file.path(repo, "Endpoints/hpcoxPPP.csv"), header = T)
hpfs_unique <- hpfs[!duplicated(hpfs$newid), t2d_cols]

nhs1 <- read.csv(file.path(repo, "Endpoints/nhs1coxPPP.csv"), header = T)
nhs1_unique <- nhs1[!duplicated(nhs1$newid), t2d_cols]

nhs2 <- read.csv(file.path(repo, "Endpoints/nhs2coxPPP.csv"), header = T)
nhs2_unique <- nhs2[!duplicated(nhs2$newid), t2d_cols]

AllCohorts_endpoints <- bind_rows(nhs1_unique, nhs2_unique, hpfs_unique)

# Exclude participants with T2D at baseline (tdb2 <= 0)
count_tdb2_baseline <- sum(AllCohorts_endpoints$tdb2 <= 0, na.rm = TRUE)
cat("Participants excluded (T2D at baseline or missing tdb2):", count_tdb2_baseline, "\n")

endpoints_ppp <- AllCohorts_endpoints[AllCohorts_endpoints$tdb2 > 0, ]

# Drop from AllCohorts3 any columns that will come from PPP files to avoid .x/.y conflicts
AllCohorts3_clean <- AllCohorts3[, !(colnames(AllCohorts3) %in% setdiff(t2d_cols, "newid"))]

# Merge metabolomics + signatures (left) with PPP endpoints (right) — no duplicate columns
AllCohorts_mortality <- merge(AllCohorts3_clean, endpoints_ppp, by = "newid", all.x = TRUE)




##########
##### Convert lifestyle score and signature to Z scores for comparable results (same unit increment)
##########

AllCohorts_mortality$aprot_avg_p_std=(AllCohorts_mortality$aprot_avg_p-mean(AllCohorts_mortality$aprot_avg_p, na.rm=T))/sd(AllCohorts_mortality$aprot_avg_p, na.rm=T)
AllCohorts_mortality$vprot_avg_p_std=(AllCohorts_mortality$vprot_avg_p-mean(AllCohorts_mortality$vprot_avg_p, na.rm=T))/sd(AllCohorts_mortality$vprot_avg_p, na.rm=T)
AllCohorts_mortality$protein_ratio_std=(AllCohorts_mortality$protein_ratio-mean(AllCohorts_mortality$protein_ratio, na.rm=T))/sd(AllCohorts_mortality$protein_ratio, na.rm=T)

AllCohorts_mortality$aprot_p_sig_std=(AllCohorts_mortality$aprot_p_sig-mean(AllCohorts_mortality$aprot_p_sig, na.rm=T))/sd(AllCohorts_mortality$aprot_p_sig, na.rm=T)
AllCohorts_mortality$vprot_p_sig_std=(AllCohorts_mortality$vprot_p_sig-mean(AllCohorts_mortality$vprot_p_sig, na.rm=T))/sd(AllCohorts_mortality$vprot_p_sig, na.rm=T)
AllCohorts_mortality$protein_ratio_sig_std=(AllCohorts_mortality$protein_ratio_sig-mean(AllCohorts_mortality$protein_ratio_sig, na.rm=T))/sd(AllCohorts_mortality$protein_ratio_sig, na.rm=T)

AllCohorts_mortality$age_2cat <- ifelse(AllCohorts_mortality$agemo >= 660, 1, 0)

### Right censoring for longevity analysis
AllCohorts_mortality$dead_85y <- ifelse(AllCohorts_mortality$dead == 1 & AllCohorts_mortality$deadage_yr < 85, 1, 0) 
# 1890 died after 85 years -> I find 90 a bit more interesting than 85y. I feel like 85y is just not premature mortality but not necessarily longevity
AllCohorts_mortality$age_end <- AllCohorts_mortality$age_end/12
AllCohorts_mortality$long85y <- ifelse((AllCohorts_mortality$dead == 0 & AllCohorts_mortality$age_end >= 85) | (AllCohorts_mortality$dead == 1 & AllCohorts_mortality$deadage_yr >= 85), 1, 0)
# SAME AS AllCohorts_mortality$age_end >= 85

AllCohorts_mortality$age_yr <- AllCohorts_mortality$agemo/12
AllCohorts_mortality$age_5yr <- AllCohorts_mortality$age_yr/5

#format factor variables
factor=c("white", "fast", "age_2cat", "multivit", "htn", "hchol", "hbpmed", "smoke_avg", "cholmed", "diab", "cohort",
         "cacoall", "endpoint")
for (i in 1:length(factor)){
  AllCohorts_mortality[, factor[i]]=as.factor(AllCohorts_mortality[, factor[i]])
}

# new_cvd <- subset(AllCohorts_mortality, tcvd >= 0)
# new_db <- subset(AllCohorts_mortality, tdb2 >= 0)
# new_chd <- subset(AllCohorts_mortality, tchd >= 0)
# new_stroke <- subset(AllCohorts_mortality, tstr >= 0)
# new_cancer <- subset(AllCohorts_mortality, tcancer >= 0)


save.image(file.path(repo, "Endpoints/Endpoints_dataset-multi-ep.RData"))


#### MEDIATION FOR DIABETES — product-of-coefficients method (Baron & Kenny)
# Exposure (X): aprot_avg_p_std (animal protein % energy, self-reported)
# Mediator (M): aprot_p_sig_std (animal protein metabolomic signature)
# Outcome (Y): incident T2D (tdb2, db2case)

library(survival)

# a path: effect of exposure on mediator
model.M <- lm(aprot_p_sig_std ~ aprot_avg_p_std + fast + white + multivit +
                htn + hchol + hbpmed + cholmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                alco_avg + calor_avg + cacoall + cohort + endpoint + agemo, data = AllCohorts_mortality)

# c path: total effect of exposure on outcome (NO mediator)
model.Y.total <- coxph(Surv(tdb2, db2case) ~ aprot_avg_p_std + fast + white + multivit +
                         htn + hchol + hbpmed + cholmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                         alco_avg + calor_avg + strata(cacoall, cohort, endpoint, agemo), data = AllCohorts_mortality)

# c' path: direct effect of exposure on outcome (WITH mediator)
model.Y.direct <- coxph(Surv(tdb2, db2case) ~ aprot_avg_p_std + fast + white + multivit +
                          htn + hchol + hbpmed + cholmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                          alco_avg + calor_avg + strata(cacoall, cohort, endpoint, agemo) + aprot_p_sig_std, data = AllCohorts_mortality)

# Coefficients (log-hazard scale)
a_path        <- coef(model.M)["aprot_avg_p_std"]           # exposure -> mediator
b_path        <- coef(model.Y.direct)["aprot_p_sig_std"]    # mediator -> outcome (adjusted)
total_effect  <- coef(model.Y.total)["aprot_avg_p_std"]     # c  (total)
direct_effect <- coef(model.Y.direct)["aprot_avg_p_std"]    # c' (direct)
indirect_effect <- a_path * b_path                          # a*b (indirect)

percent_mediated <- (indirect_effect / total_effect) * 100

cat("Total effect (c):    ", round(total_effect, 4),   " | HR =", round(exp(total_effect), 3), "\n")
cat("Direct effect (c'):  ", round(direct_effect, 4),  " | HR =", round(exp(direct_effect), 3), "\n")
cat("Indirect effect (a*b):", round(indirect_effect, 4), "\n")
cat("Percentage mediated: ", round(percent_mediated, 1), "%\n")


#################################################################################
#### TABLE 2 — Metabolite signature scores × T2D (NHS, NHSII, HPFS)         ####
#### Plant protein, Animal protein, Plant-to-animal protein ratio signatures ####
#### Model 2/3: includes alco_avg (modified AHEI excludes alcohol).           ####
#################################################################################

### --- PLANT PROTEIN SIGNATURE (vprot_p_sig_std) ---

t2d_vprot_m1 <- coxph(Surv(tdb2, db2case) ~ vprot_p_sig_std + ses + fast +
                         strata(cacoall, cohort, endpoint, agemo), data = AllCohorts_mortality)

t2d_vprot_m2 <- coxph(Surv(tdb2, db2case) ~ vprot_p_sig_std + fast + white + multivit +
                         hchol + cholmed + htn + hbpmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                         alco_avg + calor_avg + strata(cacoall, cohort, endpoint, agemo), data = AllCohorts_mortality)

t2d_vprot_m3 <- coxph(Surv(tdb2, db2case) ~ vprot_p_sig_std + fast + white + multivit +
                         hchol + cholmed + htn + hbpmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                         alco_avg + calor_avg + strata(cacoall, cohort, endpoint, agemo) + vprot_avg_p_std, data = AllCohorts_mortality)

cox.zph(t2d_vprot_m1); cox.zph(t2d_vprot_m2); cox.zph(t2d_vprot_m3)

extract_hr <- function(model) {
  coef  <- summary(model)$coefficients[1, 1]
  se    <- summary(model)$coefficients[1, 3]
  pval  <- summary(model)$coefficients[1, 5]
  hr    <- round(exp(coef), 2)
  lower <- round(exp(coef - 1.96 * se), 2)
  upper <- round(exp(coef + 1.96 * se), 2)
  result <- paste0(format(hr, nsmall=2), " (", format(lower, nsmall=2), ", ", format(upper, nsmall=2), ")")
  pval_r <- round(pval, 3)
  return(list(result=result, pval=pval_r))
}

vprot_t2d_1 <- extract_hr(t2d_vprot_m1)
vprot_t2d_2 <- extract_hr(t2d_vprot_m2)
vprot_t2d_3 <- extract_hr(t2d_vprot_m3)


### --- ANIMAL PROTEIN SIGNATURE (aprot_p_sig_std) ---

t2d_aprot_m1 <- coxph(Surv(tdb2, db2case) ~ aprot_p_sig_std + ses + fast +
                         strata(cacoall, cohort, endpoint, agemo), data = AllCohorts_mortality)

t2d_aprot_m2 <- coxph(Surv(tdb2, db2case) ~ aprot_p_sig_std + fast + white + multivit +
                         hchol + cholmed + htn + hbpmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                         alco_avg + calor_avg + strata(cacoall, cohort, endpoint, agemo), data = AllCohorts_mortality)

t2d_aprot_m3 <- coxph(Surv(tdb2, db2case) ~ aprot_p_sig_std + fast + white + multivit +
                         hchol + cholmed + htn + hbpmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                         alco_avg + calor_avg + strata(cacoall, cohort, endpoint, agemo) + aprot_avg_p_std, data = AllCohorts_mortality)

cox.zph(t2d_aprot_m1); cox.zph(t2d_aprot_m2); cox.zph(t2d_aprot_m3)

aprot_t2d_1 <- extract_hr(t2d_aprot_m1)
aprot_t2d_2 <- extract_hr(t2d_aprot_m2)
aprot_t2d_3 <- extract_hr(t2d_aprot_m3)


### --- PROTEIN RATIO SIGNATURE (protein_ratio_sig_std) ---

t2d_PR_m1 <- coxph(Surv(tdb2, db2case) ~ protein_ratio_sig_std + ses + fast +
                      strata(cacoall, cohort, endpoint, agemo), data = AllCohorts_mortality)

t2d_PR_m2 <- coxph(Surv(tdb2, db2case) ~ protein_ratio_sig_std + fast + white + multivit +
                      hchol + cholmed + htn + hbpmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                      alco_avg + calor_avg + strata(cacoall, cohort, endpoint, agemo), data = AllCohorts_mortality)

t2d_PR_m3 <- coxph(Surv(tdb2, db2case) ~ protein_ratio_sig_std + fast + white + multivit +
                      hchol + cholmed + htn + hbpmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                      alco_avg + calor_avg + strata(cacoall, cohort, endpoint, agemo) + protein_ratio_std, data = AllCohorts_mortality)

cox.zph(t2d_PR_m1); cox.zph(t2d_PR_m2); cox.zph(t2d_PR_m3)

PR_t2d_1 <- extract_hr(t2d_PR_m1)
PR_t2d_2 <- extract_hr(t2d_PR_m2)
PR_t2d_3 <- extract_hr(t2d_PR_m3)


### --- BUILD TABLE 2 ---

Table2 <- rbind(
  c("",                      "HR (95%CI)",          "P value",            "HR (95%CI)",          "P value"),
  c("Plant protein",         "",                    "",                   "",                    ""),
  c("Age-adjusted",          vprot_t2d_1$result,    vprot_t2d_1$pval,    "",                    ""),
  c("MV",                    vprot_t2d_2$result,    vprot_t2d_2$pval,    "",                    ""),
  c("MV+mutual adjustment",  vprot_t2d_3$result,    vprot_t2d_3$pval,    "",                    ""),
  c("Animal protein",        "",                    "",                   "",                    ""),
  c("Age-adjusted",          aprot_t2d_1$result,    aprot_t2d_1$pval,    "",                    ""),
  c("MV",                    aprot_t2d_2$result,    aprot_t2d_2$pval,    "",                    ""),
  c("MV+mutual adjustment",  aprot_t2d_3$result,    aprot_t2d_3$pval,    "",                    ""),
  c("Plant-to-animal ratio", "",                    "",                   "",                    ""),
  c("Age-adjusted",          PR_t2d_1$result,       PR_t2d_1$pval,       "",                    ""),
  c("MV",                    PR_t2d_2$result,       PR_t2d_2$pval,       "",                    ""),
  c("MV+mutual adjustment",  PR_t2d_3$result,       PR_t2d_3$pval,       "",                    "")
)

colnames(Table2) <- c("Model", "NHS/NHSII/HPFS HR (95%CI)", "P value", "WHI HR (95%CI)", "P value")

write.table(Table2, 
            file = file.path(repo, "Endpoints/Table2_signature_T2D.txt"),
            sep = "*", quote = FALSE, row.names = FALSE)

print(Table2)


#################################################################################
#### SUPPLEMENT — Self-reported protein intake × T2D (NHS, NHSII, HPFS)      ####
#### Same models as Table 2. Model 2/3 include alco_avg.                     ####
#### Model 3 adds the corresponding metabolomic signature.                   ####
#################################################################################

sr_vprot_m1 <- coxph(Surv(tdb2, db2case) ~ vprot_avg_p_std + ses + fast +
                       strata(cacoall, cohort, endpoint, agemo), data = AllCohorts_mortality)
sr_vprot_m2 <- coxph(Surv(tdb2, db2case) ~ vprot_avg_p_std + fast + white + multivit +
                       hchol + cholmed + htn + hbpmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                       alco_avg + calor_avg + strata(cacoall, cohort, endpoint, agemo), data = AllCohorts_mortality)
sr_vprot_m3 <- coxph(Surv(tdb2, db2case) ~ vprot_avg_p_std + fast + white + multivit +
                       hchol + cholmed + htn + hbpmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                       alco_avg + calor_avg + strata(cacoall, cohort, endpoint, agemo) + vprot_p_sig_std, data = AllCohorts_mortality)

sr_aprot_m1 <- coxph(Surv(tdb2, db2case) ~ aprot_avg_p_std + ses + fast +
                       strata(cacoall, cohort, endpoint, agemo), data = AllCohorts_mortality)
sr_aprot_m2 <- coxph(Surv(tdb2, db2case) ~ aprot_avg_p_std + fast + white + multivit +
                       hchol + cholmed + htn + hbpmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                       alco_avg + calor_avg + strata(cacoall, cohort, endpoint, agemo), data = AllCohorts_mortality)
sr_aprot_m3 <- coxph(Surv(tdb2, db2case) ~ aprot_avg_p_std + fast + white + multivit +
                       hchol + cholmed + htn + hbpmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                       alco_avg + calor_avg + strata(cacoall, cohort, endpoint, agemo) + aprot_p_sig_std, data = AllCohorts_mortality)

sr_PR_m1 <- coxph(Surv(tdb2, db2case) ~ protein_ratio_std + ses + fast +
                    strata(cacoall, cohort, endpoint, agemo), data = AllCohorts_mortality)
sr_PR_m2 <- coxph(Surv(tdb2, db2case) ~ protein_ratio_std + fast + white + multivit +
                    hchol + cholmed + htn + hbpmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                    alco_avg + calor_avg + strata(cacoall, cohort, endpoint, agemo), data = AllCohorts_mortality)
sr_PR_m3 <- coxph(Surv(tdb2, db2case) ~ protein_ratio_std + fast + white + multivit +
                    hchol + cholmed + htn + hbpmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                    alco_avg + calor_avg + strata(cacoall, cohort, endpoint, agemo) + protein_ratio_sig_std, data = AllCohorts_mortality)

sr_v1 <- extract_hr(sr_vprot_m1); sr_v2 <- extract_hr(sr_vprot_m2); sr_v3 <- extract_hr(sr_vprot_m3)
sr_a1 <- extract_hr(sr_aprot_m1); sr_a2 <- extract_hr(sr_aprot_m2); sr_a3 <- extract_hr(sr_aprot_m3)
sr_r1 <- extract_hr(sr_PR_m1);    sr_r2 <- extract_hr(sr_PR_m2);    sr_r3 <- extract_hr(sr_PR_m3)

TableS_selfreport <- rbind(
  c("",                      "HR (95%CI)",     "P value"),
  c("Plant protein",         "",               ""),
  c("Age-adjusted",          sr_v1$result,     sr_v1$pval),
  c("MV",                    sr_v2$result,     sr_v2$pval),
  c("MV+signature",          sr_v3$result,     sr_v3$pval),
  c("Animal protein",        "",               ""),
  c("Age-adjusted",          sr_a1$result,     sr_a1$pval),
  c("MV",                    sr_a2$result,     sr_a2$pval),
  c("MV+signature",          sr_a3$result,     sr_a3$pval),
  c("Plant-to-animal ratio", "",               ""),
  c("Age-adjusted",          sr_r1$result,     sr_r1$pval),
  c("MV",                    sr_r2$result,     sr_r2$pval),
  c("MV+signature",          sr_r3$result,     sr_r3$pval)
)
colnames(TableS_selfreport) <- c("Model", "HR per 1 SD (95%CI)", "P value")

write.table(TableS_selfreport,
            file = file.path(repo, "Endpoints/TableS_selfreported_T2D.txt"),
            sep = "*", quote = FALSE, row.names = FALSE)
print(TableS_selfreport)


#################################################################################
#### SENSITIVITY — signature × T2D among original controls only (cacoall == 2) ####
#### Model 1 and Model 2 as in Table 2. Model 2 includes alco_avg.            ####
#### cacoall is constant, so it is omitted from the strata.                   ####
#################################################################################

controls <- AllCohorts_mortality[AllCohorts_mortality$cacoall == 2 &
                                   !is.na(AllCohorts_mortality$cacoall), ]
cat("Controls: n =", nrow(controls),
    " events =", sum(controls$db2case == 1, na.rm = TRUE), "\n")

ctrl_v1 <- coxph(Surv(tdb2, db2case) ~ vprot_p_sig_std + ses + fast +
                   strata(cohort, endpoint, agemo), data = controls)
ctrl_a1 <- coxph(Surv(tdb2, db2case) ~ aprot_p_sig_std + ses + fast +
                   strata(cohort, endpoint, agemo), data = controls)
ctrl_r1 <- coxph(Surv(tdb2, db2case) ~ protein_ratio_sig_std + ses + fast +
                   strata(cohort, endpoint, agemo), data = controls)

ctrl_v2 <- coxph(Surv(tdb2, db2case) ~ vprot_p_sig_std + fast + white + multivit +
                   hchol + cholmed + htn + hbpmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                   alco_avg + calor_avg + strata(cohort, endpoint, agemo), data = controls)
ctrl_a2 <- coxph(Surv(tdb2, db2case) ~ aprot_p_sig_std + fast + white + multivit +
                   hchol + cholmed + htn + hbpmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                   alco_avg + calor_avg + strata(cohort, endpoint, agemo), data = controls)
ctrl_r2 <- coxph(Surv(tdb2, db2case) ~ protein_ratio_sig_std + fast + white + multivit +
                   hchol + cholmed + htn + hbpmed + ses + smoke_avg + AHEI_mod_avg + th_avg + bmi_avg +
                   alco_avg + calor_avg + strata(cohort, endpoint, agemo), data = controls)

cv1 <- extract_hr(ctrl_v1); ca1 <- extract_hr(ctrl_a1); cr1 <- extract_hr(ctrl_r1)
cv2 <- extract_hr(ctrl_v2); ca2 <- extract_hr(ctrl_a2); cr2 <- extract_hr(ctrl_r2)
TableS_controls <- rbind(
  c("Plant protein signature",         "Age-adjusted", cv1$result, cv1$pval),
  c("Plant protein signature",         "MV",           cv2$result, cv2$pval),
  c("Animal protein signature",        "Age-adjusted", ca1$result, ca1$pval),
  c("Animal protein signature",        "MV",           ca2$result, ca2$pval),
  c("Plant-to-animal ratio signature", "Age-adjusted", cr1$result, cr1$pval),
  c("Plant-to-animal ratio signature", "MV",           cr2$result, cr2$pval)
)
colnames(TableS_controls) <- c("Signature", "Model", "HR per 1 SD (95%CI)", "P value")
write.table(TableS_controls,
            file = file.path(repo, "Endpoints/TableS_T2D_controls_only.txt"),
            sep = "*", quote = FALSE, row.names = FALSE)
print(TableS_controls)

