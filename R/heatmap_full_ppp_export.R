#### Full signature-metabolite heatmap table with PPP type 2 diabetes results
#### Food, signature, mortality, longevity, and CVD columns stay as in
#### Heatmap/heatmap_multi-ep.RData. Only db_coeff2 / db_pval2 / db_pval2bon
#### are replaced from Endpoints/Results_indmetab_db2_PPP.csv.
#### Does not overwrite heatmap_multi-ep.RData or Results_indmetab_db2.csv.
####
#### Run after metabolite_t2d_cox_ppp.R:
####   source(file.path(repo, "heatmap_full_ppp_export.R"))
source("00_paths.R")
review <- file.path(repo, "Review")
if (!dir.exists(review)) dir.create(review, recursive = TRUE)

heatmap_rdata <- file.path(repo, "Heatmap/heatmap_multi-ep.RData")
ppp_csv <- file.path(repo, "Endpoints/Results_indmetab_db2_PPP.csv")
if (!file.exists(heatmap_rdata)) stop("Missing ", heatmap_rdata, call. = FALSE)
if (!file.exists(ppp_csv)) stop("Missing ", ppp_csv, " — run metabolite_t2d_cox_ppp.R first.", call. = FALSE)

load(heatmap_rdata)
if (!exists("data_heatmap_final")) stop("heatmap RData has no data_heatmap_final.", call. = FALSE)

ppp <- read.csv(ppp_csv, stringsAsFactors = FALSE)
need <- c("Metabolite", "HR.mod2", "p.value.mod2", "p.value.mod2.bon")
if (!all(need %in% names(ppp))) {
  stop("PPP CSV missing: ", paste(setdiff(need, names(ppp)), collapse = ", "), call. = FALSE)
}

hit <- match(data_heatmap_final$Metabolite, ppp$Metabolite)
data_heatmap_final$db_coeff2   <- ifelse(!is.na(hit), ppp$HR.mod2[hit],           data_heatmap_final$db_coeff2)
data_heatmap_final$db_pval2    <- ifelse(!is.na(hit), ppp$p.value.mod2[hit],      data_heatmap_final$db_pval2)
data_heatmap_final$db_pval2bon <- ifelse(!is.na(hit), ppp$p.value.mod2.bon[hit],  data_heatmap_final$db_pval2bon)

cat("Rows in full heatmap:", nrow(data_heatmap_final), "\n")
cat("T2D columns replaced from PPP:", sum(!is.na(hit)), "\n")
cat("Rows still without a PPP match:", sum(is.na(hit)), "\n")

out_full <- file.path(review, "heatmap_full_PPP.csv")
write.csv(data_heatmap_final, out_full, row.names = FALSE, na = "")

## Labeled supplemental table: estimate, P, Bonferroni P, star
## Star: * P <= 0.05; ** Bonferroni P <= 0.05. Signature weights have no P.
star_from_p <- function(p, pbon) {
  ifelse(!is.na(pbon) & pbon <= 0.05, "**", ifelse(!is.na(p) & p <= 0.05, "*", ""))
}
blocks <- c(
  PP_signature = "vprot_coef_numeric",
  AP_signature = "aprot_coef_numeric",
  PR_signature = "pr_coef_numeric",
  PP = "vprot", AP = "aprot", PR = "prot_ratio",
  Legumes = "legu", Nuts = "nuts", Refined_grains = "rgrain", Whole_grains = "wgrain",
  Processed_meat = "procm", Red_meat = "rmeat", Poultry = "poul", Fish = "fish",
  Egg = "egg", Dairy = "dairy",
  Total_mortality = "mort", Longevity = "longev",
  Type_2_diabetes = "db", CVD = "cvd"
)
miss <- character(0)
supp <- data.frame(
  HMDB = data_heatmap_final$Metabolite,
  Metabolite = if ("name" %in% names(data_heatmap_final)) data_heatmap_final$name else data_heatmap_final$Metabolite,
  stringsAsFactors = FALSE
)
if ("cate" %in% names(data_heatmap_final)) supp$Class <- data_heatmap_final$cate

for (lab in names(blocks)) {
  key <- blocks[[lab]]
  if (lab %in% c("PP_signature", "AP_signature", "PR_signature")) {
    if (!key %in% names(data_heatmap_final)) { miss <- c(miss, key); next }
    supp[[lab]] <- data_heatmap_final[[key]]
  } else {
    est <- paste0(key, "_coeff2")
    pv  <- paste0(key, "_pval2")
    bon <- paste0(key, "_pval2bon")
    if (!all(c(est, pv, bon) %in% names(data_heatmap_final))) {
      miss <- c(miss, key); next
    }
    supp[[lab]] <- data_heatmap_final[[est]]
    supp[[paste0(lab, "_P")]] <- data_heatmap_final[[pv]]
    supp[[paste0(lab, "_P_Bonferroni")]] <- data_heatmap_final[[bon]]
    supp[[paste0(lab, "_star")]] <- star_from_p(data_heatmap_final[[pv]], data_heatmap_final[[bon]])
  }
}
if (length(miss)) cat("Skipped missing columns:", paste(miss, collapse = ", "), "\n")

supp <- supp[order(supp$PP_signature, decreasing = TRUE, na.last = TRUE), ]
out_supp <- file.path(review, "SupplementaryTable_heatmap_full_PPP.csv")
write.csv(supp, out_supp, row.names = FALSE, na = "")

cat("Wrote:\n ", out_full, "\n ", out_supp, "\n", sep = "")
cat("Type 2 diabetes column is the PPP hazard ratio per 1 SD (db_coeff2 / Type_2_diabetes).\n")
cat("heatmap_multi-ep.RData was not overwritten.\n")
