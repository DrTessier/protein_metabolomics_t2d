#### Imputation and within-batch z-scoring of metabolites (discovery cohort input for elastic net)
#### Requires: merged cohort CSV and metabolite QC file at paths under `repo` (see below).
#### Output: file.path(repo, "imp_ztran_met_PR_multi_ep.RData") with object AllCohorts_clean

source("00_paths.R")

library(tidyverse)
library(dplyr)
library(glmnet)
library(gtools)
library(missRanger)
library(moments)

cohort_csv <- file.path(repo, "AllCohorts_protein_PR_excl.csv")
qc_rdata   <- file.path(repo, "readin20240606.RData")

if (!file.exists(cohort_csv)) {
  stop("Place merged discovery CSV at: ", cohort_csv, call. = FALSE)
}
if (!file.exists(qc_rdata)) {
  stop("Place metabolite QC RData (object f) at: ", qc_rdata, call. = FALSE)
}

load(qc_rdata)
gdata::keep(f, sure = TRUE)

AllCohorts <- read.csv(cohort_csv, header = TRUE, sep = ",")

hpfs <- read.csv(file.path(repo, "Endpoints/hpcox.csv"), header = TRUE)
hpfs_unique <- hpfs[!duplicated(hpfs$newid), ]

nhs1 <- read.csv(file.path(repo, "Endpoints/nhs1cox.csv"), header = TRUE)
nhs1_unique <- nhs1[!duplicated(nhs1$newid), ]

nhs2 <- read.csv(file.path(repo, "Endpoints/nhs2cox.csv"), header = TRUE)
nhs2_unique <- nhs2[!duplicated(nhs2$newid), ]

AllCohorts_endpoints <- bind_rows(nhs1_unique, nhs2_unique, hpfs_unique)
AllCohorts <- merge(AllCohorts, AllCohorts_endpoints, by = "newid", all.x = TRUE)

AllCohorts <- AllCohorts[
  (AllCohorts$tcvd > 0) & (AllCohorts$tdb2 > 0) & (AllCohorts$tstr > 0) &
    (AllCohorts$tcancer > 0) & (AllCohorts$tchd > 0),
]

all_metabolites <- data.frame(
  met_var = colnames(AllCohorts[, 2:401]),
  missing = round(colMeans(is.na(AllCohorts[, 2:401])), 4)
)
metabolites_below25_HMDB <- as.character(all_metabolites[all_metabolites$missing < 0.25, "met_var"])
metabolites_below25 <- AllCohorts[, colnames(AllCohorts) %in% metabolites_below25_HMDB]

metskewness <- data.frame(hmdb = colnames(metabolites_below25), skew = NA_real_)
for (i in seq_len(ncol(metabolites_below25))) {
  metskewness[i, 2] <- skewness(metabolites_below25[, i], na.rm = TRUE)
}
metnormal <- as.character(metskewness %>% filter(abs(skew) < 2) %>% pull(hmdb))
metskew   <- as.character(metskewness %>% filter(abs(skew) > 2) %>% pull(hmdb))
temp_normal <- metabolites_below25 %>% select(all_of(metnormal))
temp_skew   <- metabolites_below25 %>% select(all_of(metskew))
temp_skew_ln <- log(temp_skew)
temp_new <- cbind(temp_normal, temp_skew_ln)
temp_new$run_id <- AllCohorts$run_id

z_score <- function(metab) {
  if (length(metab) < 1) return(metab)
  avg <- mean(metab, na.rm = TRUE)
  std <- sd(metab, na.rm = TRUE)
  for (i in seq_along(metab)) {
    if (!is.na(metab[i])) metab[i] <- (metab[i] - avg) / std
  }
  metab
}

z_tranb <- function(dataset, metab_col) {
  z.data <- dataset
  for (j in seq_len(metab_col)) {
    z.data[, j] <- ave(dataset[, j], dataset$run_id, FUN = z_score)
  }
  z.data
}

AllCohorts_ztransformed <- z_tranb(temp_new, 243)
imputed_metabolites <- missRanger(AllCohorts_ztransformed[, 1:243], num.trees = 100, seed = 1234)
AllCohorts_clean <- cbind(imputed_metabolites, AllCohorts[, c(1, 402:508)])

AllCohorts_clean <- AllCohorts_clean[
  (AllCohorts_clean$tcvd > 0) & (AllCohorts_clean$tdb2 > 0) &
    (AllCohorts_clean$tstr > 0) & (AllCohorts_clean$tcancer > 0) & (AllCohorts_clean$tchd > 0),
]

save.image(file.path(repo, "imp_ztran_met_PR_multi_ep.RData"))
cat("Wrote", file.path(repo, "imp_ztran_met_PR_multi_ep.RData"), "\n")
