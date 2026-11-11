#### Figure 4 — Panel A: 17 shared-metabolite heatmap (same plot as figure_heatmap_shared17_vertical.R)
####                 Panel B: signature vs biomarker partial-r heatmap
#### Does not edit figure_heatmap_shared17_vertical.R or figure4_combined_vertical.R
####
#### Needs: Heatmap/heatmap_multi-ep.RData,
####        Endpoints/Results_indmetab_db2_PPP.csv (T2D column; from metabolite_t2d_cox_ppp.R),
####        Review/biomarker_partial_correlations_phyto202509.csv
#### Heatmap plot is inlined (same as figure_heatmap_shared17_vertical.R) so this file runs alone.
#### Does not overwrite heatmap_multi-ep.RData.
#### Run: source(file.path(repo, "figure4_metabolite_biomarker.R"))
####
#### Caption:
#### Figure 4. Association of shared signature metabolites with protein-rich foods,
#### type 2 diabetes, and clinical biomarkers.
#### PP, plant protein; AP, animal protein; PR, plant-to-animal protein ratio.
#### (A) Heatmap of the 17 metabolites shared by the PP, AP and PR signatures.
#### Colour is the column-standardized estimate (z-score): elastic-net signature
#### weights; multivariable linear-regression β for plant and animal foods
#### (servings/d); and log hazard ratio for incident type 2 diabetes (T2D) per
#### 1-SD metabolite. Food and T2D models were stratified by age, study cohort,
#### original sub-study, and case–control status within the original sub-study,
#### and adjusted for fasting status, ethnicity, multivitamin use, hypertension,
#### antihypertensive medication use, hypercholesterolemia, lipid-lowering
#### medication use, socioeconomic status, total energy intake, BMI, modified AHEI,
#### smoking, alcohol, and physical activity (PPP incident T2D Cox; metabolite_t2d_cox_ppp.R).
#### (B) Partial Pearson correlations of the three signatures with circulating
#### biomarkers, residualized on age, cohort, and total energy intake.
#### *P < 0.05; **Bonferroni-adjusted P < 0.05.

source("00_paths.R")
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

options(device = function(...) grDevices::pdf(file = NULL))
review <- file.path(repo, "Review")
fig_font <- 7
fig_family <- "Arial"

# =============================================================================
# Panel A — 17-metabolite heatmap (same code as figure_heatmap_shared17_vertical.R)
# =============================================================================
heatmap_rdata <- file.path(repo, "Heatmap/heatmap_multi-ep.RData")
if (!file.exists(heatmap_rdata)) stop("Missing ", heatmap_rdata, call. = FALSE)
load(heatmap_rdata)

ppp_csv <- file.path(repo, "Endpoints/Results_indmetab_db2_PPP.csv")
if (!file.exists(ppp_csv)) {
  stop("Missing ", ppp_csv, " — run metabolite_t2d_cox_ppp.R first.", call. = FALSE)
}
ppp <- read.csv(ppp_csv, stringsAsFactors = FALSE)
need_ppp <- c("Metabolite", "HR.mod2", "p.value.mod2", "p.value.mod2.bon")
if (!all(need_ppp %in% names(ppp))) {
  stop("PPP CSV missing columns: ", paste(setdiff(need_ppp, names(ppp)), collapse = ", "),
       call. = FALSE)
}
hit <- match(data_heatmap_final$Metabolite, ppp$Metabolite)
data_heatmap_final$db_coeff2   <- ifelse(!is.na(hit), ppp$HR.mod2[hit],           data_heatmap_final$db_coeff2)
data_heatmap_final$db_pval2    <- ifelse(!is.na(hit), ppp$p.value.mod2[hit],      data_heatmap_final$db_pval2)
data_heatmap_final$db_pval2bon <- ifelse(!is.na(hit), ppp$p.value.mod2.bon[hit],  data_heatmap_final$db_pval2bon)
cat("PPP T2D HRs overlaid on", sum(!is.na(hit)), "metabolites (heatmap RData not overwritten).\n")

in_all_three <- with(
  data_heatmap_final,
  !is.na(vprot_coef_numeric) & !is.na(aprot_coef_numeric) & !is.na(pr_coef_numeric)
)
shared <- data_heatmap_final[in_all_three, ]
cat("Shared metabolites:", nrow(shared),
    "  T2D HRs from PPP:", sum(!is.na(shared$db_coeff2)), "\n")

label_col <- if ("name" %in% names(shared) && any(!is.na(shared$name))) "name" else "Metabolite"
met_label <- shared[[label_col]]
met_label[is.na(met_label) | met_label == ""] <- shared$Metabolite
met_label <- make.unique(as.character(met_label))

cols_w <- c("PP signature", "AP signature", "PR signature")
cols_plant <- c("Legumes", "Nuts", "Refined grains", "Whole grains")
cols_animal <- c("Processed meat", "Red meat", "Poultry", "Fish", "Eggs", "Dairy")
col_vars <- c(cols_w, cols_plant, cols_animal, "Type 2 diabetes")
food_p <- c("legu", "nuts", "rgrain", "wgrain", "procm", "rmeat", "poul", "fish", "egg", "dairy")

scale_cols <- function(m) {
  m <- scale(m, center = FALSE)
  m[is.na(m)] <- 0
  m
}
mat_w <- as.matrix(shared[, c("vprot_coef_numeric", "aprot_coef_numeric", "pr_coef_numeric")])
mat_f <- as.matrix(shared[, paste0(food_p, "_coeff2")])
mat_t2d <- log(as.matrix(shared[, "db_coeff2", drop = FALSE]))
mat <- cbind(scale_cols(mat_w), scale_cols(mat_f), scale_cols(mat_t2d))
colnames(mat) <- col_vars

star_from_p <- function(p, pbon) {
  ifelse(!is.na(pbon) & pbon <= 0.05, "**", ifelse(!is.na(p) & p <= 0.05, "*", ""))
}
stars <- cbind(
  matrix("", nrow(shared), length(cols_w)),
  sapply(c(food_p, "db"), function(f) star_from_p(shared[[paste0(f, "_pval2")]],
                                                  shared[[paste0(f, "_pval2bon")]]))
)
colnames(stars) <- col_vars
met_ord <- met_label[order(shared$vprot_coef_numeric, decreasing = TRUE)]
met_levels <- rev(met_ord)
n_met <- length(met_levels)

# Same layout as panel B: unit spacing within a group, 1.5-unit gap between groups
x_pos <- setNames(
  c(1:3, 4.5 + 0:3, 9.0 + 0:5, 15.5),
  col_vars
)
grp_df_a <- data.frame(
  group = c("Signature weights", "Plant foods", "Animal foods"),
  xmin  = c(1, 4.5, 9.0) - 0.5,
  xmax  = c(3, 7.5, 14.0) + 0.5,
  stringsAsFactors = FALSE
)
grp_df_a$xmid <- (grp_df_a$xmin + grp_df_a$xmax) / 2
x_lim_shared <- c(min(x_pos) - 0.5, max(x_pos) + 0.5)
x_exp_shared <- expansion(mult = 0)
xlab_h <- 2.6
bar_h  <- 0.9
bar_gap <- 0.05

heat_dat_v <- data.frame(
  metabolite = factor(rep(met_label, times = length(col_vars)), levels = met_levels),
  variable   = factor(rep(col_vars, each = n_met), levels = col_vars),
  value      = as.vector(mat),
  star       = as.vector(stars),
  stringsAsFactors = FALSE
)
heat_dat_v$x <- unname(x_pos[as.character(heat_dat_v$variable)])
heat_dat_v$y <- as.integer(heat_dat_v$metabolite)

lim <- min(max(abs(heat_dat_v$value), na.rm = TRUE), 2.5)
heat_colours <- c("#0046FF", "#5876FF", "#7F97FF", "#A9B9FF", "#D4DBFF", "#FFFFFF",
                  "#FFE9F3", "#FFB0D6", "#FF7FBC", "#FF5AA8", "#FF2F92")

heat_axis_theme <- theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_text(size = fig_font - 1, family = fig_family, colour = "black"),
  legend.position = "bottom",
  legend.title = element_text(size = fig_font, family = fig_family),
  legend.text = element_text(size = fig_font - 1, family = fig_family),
  plot.margin = margin(3, 4, 2, 7, "mm")
)

header_ys <- function(n_row, lab_h = xlab_h) {
  y_xlab <- n_row + 0.55
  y_bar_min <- y_xlab + lab_h + bar_gap
  y_bar_max <- y_bar_min + bar_h
  list(xlab = y_xlab, bar_min = y_bar_min, bar_max = y_bar_max,
       top = y_bar_max + 0.12, extra = y_bar_max + 0.12 - 0.5)
}

ys_a <- header_ys(n_met)
grp_df_a$ymin <- ys_a$bar_min
grp_df_a$ymax <- ys_a$bar_max
grp_df_a$ytxt <- (ys_a$bar_min + ys_a$bar_max) / 2
col_lab_a <- data.frame(
  x = unname(x_pos), lab = names(x_pos), y = ys_a$xlab,
  stringsAsFactors = FALSE
)

p_heat_a <- ggplot(heat_dat_v, aes(x = x, y = y, fill = value)) +
  geom_tile(width = 1, height = 1, colour = NA) +
  geom_text(aes(label = star), size = 6 / .pt, vjust = 0.75, family = fig_family) +
  geom_text(data = col_lab_a, aes(x = x, y = y, label = lab),
            inherit.aes = FALSE, angle = 90, hjust = 0, vjust = 0.5,
            size = (fig_font - 1) / .pt, family = fig_family) +
  geom_rect(data = grp_df_a,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "grey82", colour = "white", linewidth = 0.3) +
  geom_text(data = grp_df_a,
            aes(x = xmid, y = ytxt, label = group),
            inherit.aes = FALSE, size = 6 / .pt, family = fig_family) +
  scale_x_continuous(limits = x_lim_shared, expand = x_exp_shared) +
  scale_y_continuous(breaks = seq_len(n_met), labels = met_levels,
                     limits = c(0.5, ys_a$top), expand = c(0, 0)) +
  scale_fill_gradientn(colours = heat_colours, limits = c(-lim, lim), breaks = -2:2,
                       oob = scales::squish, name = "Z-score") +
  coord_fixed(ratio = 1, clip = "off") +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                barwidth = unit(40, "mm"), barheight = unit(2.5, "mm"),
                                frame.colour = NA, ticks.colour = "white")) +
  theme_classic(base_size = fig_font, base_family = fig_family) +
  heat_axis_theme

# =============================================================================
# Panel B — partial r heatmap (available biomarkers only)
# =============================================================================
pc_file <- file.path(review, "biomarker_partial_correlations_phyto202509.csv")
if (!file.exists(pc_file)) {
  stop("Run biomarker_phyto202509.R first (missing ", pc_file, ").", call. = FALSE)
}
pc <- read.csv(pc_file, stringsAsFactors = FALSE)

sig_row <- c(
  ratio  = "Plant-to-animal\nratio signature",
  animal = "Animal\nsignature",
  plant  = "Plant\nsignature"
)
sig_y <- c(ratio = 1, animal = 2, plant = 3)
n_sig <- length(sig_y)

bio_meta <- data.frame(
  biomarker_label = c("HbA1c", "C-peptide", "TAG/HDL-C", "CRP",
                      "Adiponectin", "Leptin", "Total chol", "LDL"),
  bio_short = c("HbA1c", "C-peptide", "TG/HDL", "CRP",
                "Adiponectin", "Leptin", "Tot chol", "LDL"),
  group = c("Glycaemia & IR", "Glycaemia & IR", "Glycaemia & IR",
            "Inflammation", "Adipokines", "Adipokines", "Lipids", "Lipids"),
  x = c(1, 2, 3, 5.6, 8.0, 9.0, 10.6, 11.6),
  stringsAsFactors = FALSE
)
heat_dat <- pc %>%
  inner_join(bio_meta, by = "biomarker_label") %>%
  mutate(
    y = unname(sig_y[signature]),
    lab = star_from_p(p_value, p.adjust(p_value, method = "bonferroni"))
  )

bio_lim <- ceiling(max(abs(heat_dat$partial_r), na.rm = TRUE) * 20) / 20
if (!is.finite(bio_lim) || bio_lim < 0.2) bio_lim <- 0.25

grp_df <- bio_meta %>%
  group_by(group) %>%
  summarise(xmin = min(x) - 0.5, xmax = max(x) + 0.5, .groups = "drop") %>%
  mutate(
    xmid = (xmin + xmax) / 2,
    half = ifelse((xmax - xmin) < 1.2, 1.25, (xmax - xmin) / 2),
    xmin = xmid - half,
    xmax = xmid + half,
    group = factor(group, levels = c("Glycaemia & IR", "Inflammation",
                                     "Adipokines", "Lipids"))
  )

ys_b <- header_ys(n_sig, lab_h = 1.95)
grp_df$ymin <- ys_b$bar_min
grp_df$ymax <- ys_b$bar_max
grp_df$ytxt <- (ys_b$bar_min + ys_b$bar_max) / 2
col_lab_b <- data.frame(
  x = bio_meta$x, lab = bio_meta$bio_short, y = ys_b$xlab,
  stringsAsFactors = FALSE
)

p_heat_b <- ggplot(heat_dat, aes(x = x, y = y, fill = partial_r)) +
  geom_tile(width = 1, height = 1, colour = NA) +
  geom_text(aes(label = lab), size = 6 / .pt, vjust = 0.75, family = fig_family) +
  geom_text(data = col_lab_b, aes(x = x, y = y, label = lab),
            inherit.aes = FALSE, angle = 90, hjust = 0, vjust = 0.5,
            size = (fig_font - 1) / .pt, family = fig_family) +
  geom_rect(data = grp_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "grey82", colour = "white", linewidth = 0.3) +
  geom_text(data = grp_df,
            aes(x = xmid, y = ytxt, label = group),
            inherit.aes = FALSE, size = 6 / .pt, family = fig_family) +
  scale_x_continuous(limits = x_lim_shared, expand = x_exp_shared) +
  scale_y_continuous(breaks = unname(sig_y), labels = unname(sig_row),
                     limits = c(0.5, ys_b$top), expand = c(0, 0)) +
  scale_fill_gradientn(
    colours = heat_colours, limits = c(-bio_lim, bio_lim), oob = squish,
    name = "Partial r"
  ) +
  coord_fixed(ratio = 1, clip = "off") +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                barwidth = unit(40, "mm"), barheight = unit(2.5, "mm"),
                                frame.colour = NA, ticks.colour = "white")) +
  theme_classic(base_size = fig_font, base_family = fig_family) +
  heat_axis_theme +
  theme(axis.text.y = element_text(size = fig_font, family = fig_family,
                                   colour = "black", lineheight = 0.95))

# =============================================================================
# Combine: A = 17-metabolite heatmap, B = biomarkers
# Same x limits; grey titles drawn in the same plot as the tiles so they share x
# =============================================================================
fig4 <- (p_heat_a / p_heat_b) +
  plot_layout(heights = c(ys_a$extra, ys_b$extra)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 9, face = "bold", family = fig_family,
                            hjust = 1, vjust = 1),
    plot.tag.location = "panel",
    plot.tag.position = c(-0.06, 1),
    plot.margin = margin(3, 4, 2, 7, "mm")
  )

# npj Science of Food: RGB, >=300 dpi; prefer vector PDF
fig_w <- 180
fig_h <- 210
stem <- file.path(review, "Figure4_metabolite_biomarker")

save_npj_figure <- function(plot, stem, width_mm, height_mm) {
  pdf_file <- paste0(stem, ".pdf")
  if (capabilities("cairo")) {
    ggsave(pdf_file, plot, width = width_mm, height = height_mm, units = "mm",
           device = grDevices::cairo_pdf, bg = "white", family = fig_family)
  } else {
    ggsave(pdf_file, plot, width = width_mm, height = height_mm, units = "mm", bg = "white")
  }

  tiff_fun <- function(file, dpi) {
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_tiff(file, width = width_mm, height = height_mm, units = "mm",
                     res = dpi, compression = "lzw", background = "white")
      print(plot)
      grDevices::dev.off()
    } else {
      ggsave(file, plot, width = width_mm, height = height_mm, units = "mm",
             dpi = dpi, device = "tiff", compression = "lzw", bg = "white")
    }
  }
  tiff_fun(paste0(stem, "_300dpi.tiff"), 300)
  tiff_fun(paste0(stem, "_600dpi.tiff"), 600)

  png_file <- paste0(stem, "_300dpi.png")
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(png_file, plot, width = width_mm, height = height_mm, units = "mm",
           dpi = 300, device = ragg::agg_png, bg = "white")
  } else {
    ggsave(png_file, plot, width = width_mm, height = height_mm, units = "mm",
           dpi = 300, type = "cairo", bg = "white")
  }
  file.copy(png_file, paste0(stem, ".png"), overwrite = TRUE)
  cat("Saved (npj Science of Food):\n ", pdf_file, "  [vector — use this for submission]\n ",
      paste0(stem, "_300dpi.tiff"), "\n ", paste0(stem, "_600dpi.tiff"), "\n ",
      png_file, "\n", sep = "")
}

save_npj_figure(fig4, stem, fig_w, fig_h)
