#### Figure 3 — Panel A: signature vs self-reported intake (forests)
####                 Panel B: food-group dominance (same plot as figure4_combined_vertical.R)
#### Does not edit figure4_combined_vertical.R or food_contribution_signatures.R
####
#### Needs: Enet + WHI signature files (as in correlations.R),
####        Review/food_contribution_signatures.csv
#### Run: source(file.path(repo, "figure3_intake_dominance.R"))

source("00_paths.R")
library(ggplot2)
library(dplyr)
library(patchwork)

options(device = function(...) grDevices::pdf(file = NULL))
review <- file.path(repo, "Review")
fig_font <- 7
# npj Science of Food: Arial or Helvetica (not Calibri)
fig_family <- "Arial"

if (!dir.exists(review)) dir.create(review, recursive = TRUE)

# =============================================================================
# Panel A — Pearson r (95% CI) vs self-reported intake, by cohort
# Same files, columns, and cohort splits as correlations.R (NHS, NHSII, HPFS, WHI only)
# =============================================================================
sig_facet <- c(
  plant  = "Plant protein\nsignature",
  animal = "Animal protein\nsignature",
  ratio  = "Plant-to-animal\nratio signature"
)

# Paper Figure 2 scatter (correlations.R labels / manuscript)
paper_r <- data.frame(
  signature = c(rep("plant", 4), rep("animal", 4), rep("ratio", 4)),
  cohort    = rep(c("NHS", "NHSII", "HPFS", "WHI"), 3),
  r_paper   = c(0.35, 0.39, 0.48, 0.31,
                0.26, 0.28, 0.34, 0.23,
                0.28, 0.36, 0.41, 0.25),
  stringsAsFactors = FALSE
)

cor_one <- function(x, y, cohort, signature) {
  ct <- cor.test(x, y)
  data.frame(
    signature = signature,
    cohort = cohort,
    n = sum(complete.cases(x, y)),
    r = unname(ct$estimate),
    lo = ct$conf.int[1],
    hi = ct$conf.int[2],
    stringsAsFactors = FALSE
  )
}

format_cohort <- function(x) {
  dplyr::case_when(
    as.character(x) %in% c("nhs1", "NHS")   ~ "NHS",
    as.character(x) %in% c("nhs2", "NHSII") ~ "NHSII",
    as.character(x) %in% c("hpfs", "HPFS")  ~ "HPFS",
    as.character(x) == "WHI"                ~ "WHI",
    TRUE ~ as.character(x)
  )
}

# Plant — identical objects/columns to correlations.R
load(file.path(repo, "WHI-replication/vprot_predicted_score_WHI.RData"))
signature_vprot$vprot_p_sig <- signature_vprot$sig
signature_vprot$newid <- signature_vprot$ID
signature_vprot$cohort <- "WHI"
v_whi <- signature_vprot[, c("newid", "vprot_avg_p", "cohort", "vprot_p_sig")]
v_har <- read.csv(file.path(repo, "Enet/Plant/signature_vprot_0606-PR.csv"))
v_har <- v_har[, c("newid", "vprot_avg_p", "cohort", "vprot_p_sig")]
v_all <- rbind(v_har, v_whi)
v_all$cohort <- format_cohort(v_all$cohort)

# Animal
load(file.path(repo, "WHI-replication/aprot_predicted_score_WHI.RData"))
signature$aprot_p_sig <- signature$sig
signature$newid <- signature$ID
signature$cohort <- "WHI"
a_whi <- signature[, c("newid", "aprot_avg_p", "cohort", "aprot_p_sig")]
a_har <- read.csv(file.path(repo, "Enet/Animal/signature_aprot_0606-PR.csv"))
a_har <- a_har[, c("newid", "aprot_avg_p", "cohort", "aprot_p_sig")]
a_all <- rbind(a_har, a_whi)
a_all$cohort <- format_cohort(a_all$cohort)

# Ratio
load(file.path(repo, "WHI-replication/PR_predicted_score_WHI.RData"))
signature_pr$protein_ratio_sig <- signature_pr$sig
signature_pr$newid <- signature_pr$ID
signature_pr$cohort <- "WHI"
r_whi <- signature_pr[, c("newid", "protein_ratio", "cohort", "protein_ratio_sig")]
r_har <- read.csv(file.path(repo, "Enet/Protein_ratio/signature_pr_0606.csv"))
r_har <- r_har[, c("newid", "protein_ratio", "cohort", "protein_ratio_sig")]
r_all <- rbind(r_har, r_whi)
r_all$cohort <- format_cohort(r_all$cohort)

coh_ord <- c("Pooled", "NHS", "NHSII", "HPFS", "WHI")
rows <- list()
for (coh in c("NHS", "NHSII", "HPFS", "WHI")) {
  vv <- v_all[v_all$cohort == coh, ]
  aa <- a_all[a_all$cohort == coh, ]
  rr <- r_all[r_all$cohort == coh, ]
  rows[[length(rows) + 1]] <- cor_one(vv$vprot_avg_p, vv$vprot_p_sig, coh, "plant")
  rows[[length(rows) + 1]] <- cor_one(aa$aprot_avg_p, aa$aprot_p_sig, coh, "animal")
  rows[[length(rows) + 1]] <- cor_one(rr$protein_ratio, rr$protein_ratio_sig, coh, "ratio")
}
# Pooled = NHS + NHSII + HPFS (not WHI); same Pearson r + Fisher z 95% CI
v_disc <- v_all[v_all$cohort != "WHI", ]
a_disc <- a_all[a_all$cohort != "WHI", ]
r_disc <- r_all[r_all$cohort != "WHI", ]
rows[[length(rows) + 1]] <- cor_one(v_disc$vprot_avg_p, v_disc$vprot_p_sig, "Pooled", "plant")
rows[[length(rows) + 1]] <- cor_one(a_disc$aprot_avg_p, a_disc$aprot_p_sig, "Pooled", "animal")
rows[[length(rows) + 1]] <- cor_one(r_disc$protein_ratio, r_disc$protein_ratio_sig, "Pooled", "ratio")

cor_dat <- do.call(rbind, rows)
cor_dat$cohort <- factor(cor_dat$cohort, levels = coh_ord)
cor_dat$signature <- factor(cor_dat$signature, levels = c("plant", "animal", "ratio"))
cor_dat$is_whi <- cor_dat$cohort == "WHI"
cor_dat$r_2dec <- round(cor_dat$r, 2)

chk <- merge(cor_dat[cor_dat$cohort != "Pooled", ], paper_r, by = c("signature", "cohort"), sort = FALSE)
cat("\n=== Intake correlations vs paper Figure 2 ===\n")
print(chk[, c("signature", "cohort", "n", "r", "r_2dec", "r_paper", "lo", "hi")], row.names = FALSE)
if (any(chk$r_2dec != chk$r_paper)) {
  warning("Rounded r does not match the Figure 2 scatter / manuscript values.", call. = FALSE)
  print(chk[chk$r_2dec != chk$r_paper, c("signature", "cohort", "r", "r_2dec", "r_paper")],
        row.names = FALSE)
}
cat("\n=== Pooled (NHS+NHSII+HPFS) ===\n")
print(cor_dat[cor_dat$cohort == "Pooled", c("signature", "n", "r", "r_2dec", "lo", "hi")],
      row.names = FALSE)

write.csv(cor_dat, file.path(review, "figure3_intake_correlations.csv"), row.names = FALSE)

make_forest <- function(sig, show_y = FALSE) {
  d <- cor_dat[cor_dat$signature == sig, ]
  p <- ggplot(d, aes(x = cohort, y = r, colour = is_whi)) +
    geom_hline(yintercept = seq(0, 0.6, 0.1), colour = "grey85", linewidth = 0.3) +
    geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.5) +
    geom_point(size = 1.9) +
    scale_colour_manual(
      values = c("FALSE" = "#0000CD", "TRUE" = "#E69F00"),
      labels = c("FALSE" = "Discovery", "TRUE" = "External"),
      name = NULL
    ) +
    scale_x_discrete(expand = expansion(mult = 0.14)) +
    scale_y_continuous(breaks = seq(0, 0.6, 0.1)) +
    coord_cartesian(ylim = c(0, 0.6)) +
    labs(title = unname(sig_facet[sig]), x = NULL, y = NULL) +
    theme_classic(base_size = fig_font, base_family = fig_family) +
    theme(
      plot.title = element_text(size = fig_font, face = "bold", hjust = 0.5, lineheight = 0.95),
      panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.4),
      axis.line = element_blank(),
      axis.text = element_text(size = fig_font, colour = "black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = fig_font),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(t = 0, b = 1, unit = "mm"),
      plot.margin = margin(1, 2, 1, 1, "mm")
    )
  if (show_y) {
    p <- p +
      labs(y = "Pearson r") +
      theme(
        axis.text.y = element_text(size = fig_font, colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        axis.title.y = element_text(size = fig_font, angle = 90, vjust = 0.5, hjust = 0.5,
                                    margin = margin(r = 2, unit = "mm"))
      )
  }
  p
}

p_forests <- (make_forest("plant", TRUE) | make_forest("animal") | make_forest("ratio")) +
  plot_layout(guides = "collect", tag_level = "new") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.margin = margin(t = 0, b = 1, unit = "mm")
  )

# =============================================================================
# Panel B — food-group dominance (copied from figure4_combined_vertical.R, plot unchanged)
# =============================================================================
dom_csv <- file.path(review, "food_contribution_signatures.csv")
if (!file.exists(dom_csv)) {
  stop("Run food_contribution_signatures.R first (missing ", dom_csv, ").", call. = FALSE)
}
tab <- read.csv(dom_csv)

food_order <- c("legumes", "nuts", "refined_grains", "whole_grains",
                "processed_meat", "red_meat", "poultry", "fish", "egg", "dairy")
food_labels <- c(legumes = "Legumes", nuts = "Nuts", whole_grains = "Whole grains",
                 refined_grains = "Refined grains", processed_meat = "Processed meat",
                 red_meat = "Red meat", poultry = "Poultry", fish = "Fish",
                 egg = "Eggs", dairy = "Dairy")
food_colors <- c(legumes = "#B5005F", nuts = "#FF1493", whole_grains = "#FF7FBC",
                 refined_grains = "#FFC4E1", processed_meat = "#00008B", red_meat = "#0000CD",
                 poultry = "#5876FF", fish = "#7F97FF", egg = "#A9B9FF", dairy = "#D4DBFF")
dark_text_foods <- c("whole_grains", "refined_grains", "fish", "egg", "dairy")
signature_order <- c("plant", "animal", "ratio")
signature_labels <- c(plant  = "Plant protein\nsignature",
                      animal = "Animal protein\nsignature",
                      ratio  = "Plant-to-animal\nratio signature")

plot_dat <- tab %>%
  filter(signature %in% signature_order) %>%
  mutate(food_group = factor(food_group, levels = rev(food_order)),
         signature  = factor(signature, levels = signature_order)) %>%
  arrange(signature, food_group) %>%
  group_by(signature) %>%
  mutate(pct  = 100 * R2_contribution / sum(R2_contribution),
         ymax = cumsum(pct), ymin = ymax - pct, ymid = (ymin + ymax) / 2) %>%
  ungroup() %>%
  mutate(x        = as.numeric(signature),
         sign     = ifelse(direction == "+", "+", "\u2212"),
         label    = ifelse(pct >= 3, paste0(round(pct), "% ", sign), ""),
         text_col = ifelse(food_group %in% dark_text_foods, "black", "white"))

r2_dat <- plot_dat %>%
  distinct(signature, x, R2_foods_total) %>%
  mutate(label = sprintf("R\u00b2 = %.3f", R2_foods_total))

make_dom_bar <- function(sig, show_y = FALSE) {
  d <- plot_dat[as.character(plot_dat$signature) == sig, ]
  r2_lab <- r2_dat$label[as.character(r2_dat$signature) == sig]
  p <- ggplot(d) +
    geom_rect(aes(xmin = 0.62, xmax = 1.38, ymin = ymin, ymax = ymax, fill = food_group),
              colour = "white", linewidth = 0.3) +
    geom_text(aes(x = 1, y = ymid, label = label, colour = text_col),
              size = 6 / .pt, family = fig_family) +
    scale_colour_identity() +
    scale_fill_manual(values = food_colors, labels = food_labels, breaks = food_order, name = NULL) +
    scale_x_continuous(breaks = 1, labels = r2_lab, limits = c(0.15, 1.85), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"),
                       expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_legend(ncol = 2)) +
    theme_classic(base_size = fig_font, base_family = fig_family) +
    theme(
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text = element_text(size = fig_font, colour = "black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = fig_font),
      legend.key.size = unit(3, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(2, 2, 1, 1, "mm")
    )
  if (show_y) {
    p <- p +
      labs(y = "Share of variance explained by food groups (%)") +
      theme(
        axis.text.y = element_text(size = fig_font, colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        axis.title.y = element_text(size = fig_font, angle = 90, vjust = 0.5, hjust = 0.5,
                                    margin = margin(r = 1.5, unit = "mm"))
      )
  } else {
    p <- p + theme(axis.line.y = element_blank())
  }
  p
}

p_bars <- (make_dom_bar("plant", TRUE) | make_dom_bar("animal") | make_dom_bar("ratio")) +
  plot_layout(guides = "collect", tag_level = "new") &
  theme(legend.position = "bottom", legend.direction = "horizontal")

# =============================================================================
# Combine and save — A = forests, B = dominance; same 3 columns so bars sit under A
# =============================================================================
panel_a <- wrap_elements(p_forests) +
  labs(tag = "A") +
  theme(plot.tag = element_text(size = 9, face = "bold", family = fig_family),
        plot.tag.position = c(0, 1),
        plot.margin = margin(2, 2, 1, 6, "mm"))
panel_b <- wrap_elements(p_bars) +
  labs(tag = "B") +
  theme(plot.tag = element_text(size = 9, face = "bold", family = fig_family),
        plot.tag.position = c(0, 1),
        plot.margin = margin(1, 2, 2, 6, "mm"))

fig2 <- (panel_a / panel_b) + plot_layout(heights = c(1, 1.15))

# npj Science of Food (Nature Portfolio): RGB, >=300 dpi; prefer vector PDF/EPS.
# https://www.nature.com/npjscifood/for-authors-and-referees/submission-guidelines
# Export PDF (vector, preferred) + TIFF 300/600 dpi LZW + PNG 300 dpi so we do not re-run later.
fig_w <- 180
fig_h <- 168
stem <- file.path(review, "Figure3_intake_dominance")

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
  # Keep a short name for checking in Review/
  file.copy(png_file, paste0(stem, ".png"), overwrite = TRUE)
  cat("Saved (npj Science of Food):\n ", pdf_file, "  [vector — use this for submission]\n ",
      paste0(stem, "_300dpi.tiff"), "\n ", paste0(stem, "_600dpi.tiff"), "\n ",
      png_file, "\n", sep = "")
}

save_npj_figure(fig2, stem, fig_w, fig_h)
