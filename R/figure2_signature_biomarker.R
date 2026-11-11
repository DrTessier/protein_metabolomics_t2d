#### Figure 2 — signature vs intake (Panel A) and signature vs biomarkers (Panel B)
#### Stacked: 3 forest plots on top, grouped partial-r heatmap below (full 180 mm width)
#### Biomarkers are those actually available (no WC, glucose, HOMA-IR, TNF-α, HDL, TAG)
####
#### Run on Channing after biomarker_phyto202509.R has written the partial-r CSV:
####   source(file.path(repo, "figure2_signature_biomarker.R"))

source("00_paths.R")
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

options(device = function(...) grDevices::pdf(file = NULL))
review <- file.path(repo, "Review")
fig_font <- 7

sig_facet <- c(
  plant  = "Plant signature",
  animal = "Animal signature",
  ratio  = "Plant-to-animal ratio"
)
sig_row <- c(
  plant  = "Plant\nsignature",
  animal = "Animal\nsignature",
  ratio  = "Plant-to-animal\nratio signature"
)
y_lab_width <- "Plant-to-animal\nratio signature"
heat_colours <- c("#0046FF", "#5876FF", "#7F97FF", "#A9B9FF", "#D4DBFF", "#FFFFFF",
                  "#FFE9F3", "#FFB0D6", "#FF7FBC", "#FF5AA8", "#FF2F92")

# =============================================================================
# Panel A — Pearson r (95% CI) of signature vs self-reported intake, by cohort
# Same files as correlations.R, plus a pooled discovery estimate (NHS+NHSII+HPFS)
# =============================================================================
cor_one <- function(x, y, cohort, signature) {
  ok <- complete.cases(x, y)
  ct <- cor.test(x[ok], y[ok])
  data.frame(
    signature = signature,
    cohort = cohort,
    n = sum(ok),
    r = unname(ct$estimate),
    lo = ct$conf.int[1],
    hi = ct$conf.int[2],
    stringsAsFactors = FALSE
  )
}

recode_coh <- function(x) {
  dplyr::recode(as.character(x),
                nhs1 = "NHS", nhs2 = "NHSII", hpfs = "HPFS",
                NHS = "NHS", NHSII = "NHSII", HPFS = "HPFS", WHI = "WHI")
}

load_intake_cors <- function() {
  load(file.path(repo, "WHI-replication/vprot_predicted_score_WHI.RData"))
  signature_vprot$vprot_p_sig <- signature_vprot$sig
  signature_vprot$cohort <- "WHI"
  v_whi <- signature_vprot[, c("vprot_avg_p", "cohort", "vprot_p_sig")]
  v_har <- read.csv(file.path(repo, "Enet/Plant/signature_vprot_0606-PR.csv"))
  v_har <- v_har[, c("vprot_avg_p", "cohort", "vprot_p_sig")]
  v_all <- rbind(v_har, v_whi)
  v_all$cohort <- recode_coh(v_all$cohort)

  load(file.path(repo, "WHI-replication/aprot_predicted_score_WHI.RData"))
  signature$aprot_p_sig <- signature$sig
  signature$cohort <- "WHI"
  a_whi <- signature[, c("aprot_avg_p", "cohort", "aprot_p_sig")]
  a_har <- read.csv(file.path(repo, "Enet/Animal/signature_aprot_0606-PR.csv"))
  a_har <- a_har[, c("aprot_avg_p", "cohort", "aprot_p_sig")]
  a_all <- rbind(a_har, a_whi)
  a_all$cohort <- recode_coh(a_all$cohort)

  load(file.path(repo, "WHI-replication/PR_predicted_score_WHI.RData"))
  signature_pr$protein_ratio_sig <- signature_pr$sig
  signature_pr$cohort <- "WHI"
  r_whi <- signature_pr[, c("protein_ratio", "cohort", "protein_ratio_sig")]
  r_har <- read.csv(file.path(repo, "Enet/Protein_ratio/signature_pr_0606.csv"))
  r_har <- r_har[, c("protein_ratio", "cohort", "protein_ratio_sig")]
  r_all <- rbind(r_har, r_whi)
  r_all$cohort <- recode_coh(r_all$cohort)

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
  v_disc <- v_all[v_all$cohort != "WHI", ]
  a_disc <- a_all[a_all$cohort != "WHI", ]
  r_disc <- r_all[r_all$cohort != "WHI", ]
  rows[[length(rows) + 1]] <- cor_one(v_disc$vprot_avg_p, v_disc$vprot_p_sig, "Pooled", "plant")
  rows[[length(rows) + 1]] <- cor_one(a_disc$aprot_avg_p, a_disc$aprot_p_sig, "Pooled", "animal")
  rows[[length(rows) + 1]] <- cor_one(r_disc$protein_ratio, r_disc$protein_ratio_sig, "Pooled", "ratio")

  out <- do.call(rbind, rows)
  out$cohort <- factor(out$cohort, levels = coh_ord)
  out$signature <- factor(out$signature, levels = c("plant", "animal", "ratio"))
  out$is_whi <- out$cohort == "WHI"
  out
}

if (!dir.exists(review)) dir.create(review, recursive = TRUE)
cor_dat <- load_intake_cors()
write.csv(cor_dat, file.path(review, "figure2_intake_correlations.csv"), row.names = FALSE)
print(cor_dat, row.names = FALSE)

pooled_r <- cor_dat[cor_dat$cohort == "Pooled", c("signature", "r")]

make_forest <- function(sig, show_y = FALSE) {
  d <- cor_dat[cor_dat$signature == sig, ]
  href <- pooled_r$r[as.character(pooled_r$signature) == sig]
  p <- ggplot(d, aes(x = cohort, y = r, colour = is_whi)) +
    geom_hline(yintercept = href, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.5) +
    geom_point(size = 1.9) +
    scale_colour_manual(
      values = c("FALSE" = "#0000CD", "TRUE" = "#E69F00"),
      labels = c("FALSE" = "Discovery", "TRUE" = "WHI (external)"),
      name = NULL
    ) +
    scale_x_discrete(expand = expansion(mult = 0.14)) +
    scale_y_continuous(breaks = seq(0.2, 0.5, 0.1)) +
    coord_cartesian(ylim = c(0.15, 0.55)) +
    labs(title = unname(sig_facet[sig]), x = NULL, y = NULL) +
    theme_classic(base_size = fig_font) +
    theme(
      plot.title = element_text(size = fig_font, face = "bold", hjust = 0.5),
      panel.border = element_rect(colour = "grey35", fill = NA, linewidth = 0.4),
      axis.line = element_blank(),
      axis.text = element_text(size = fig_font, colour = "black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = fig_font),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(1, 2, 1, 1, "mm")
    )
  if (show_y) {
    p <- p +
      labs(y = "Pearson r") +
      theme(
        axis.text.y = element_text(size = fig_font, colour = "black"),
        axis.ticks.y = element_line(colour = "black"),
        axis.title.y = element_text(size = fig_font)
      )
  }
  p
}

p_a <- (make_forest("plant", TRUE) | make_forest("animal") | make_forest("ratio")) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Agreement with self-reported intake, by cohort (Pearson r with 95% CIs)",
    theme = theme(plot.title = element_text(size = fig_font, face = "bold", hjust = 0))
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

# =============================================================================
# Panel B — partial r heatmap (available biomarkers only)
# =============================================================================
pc_file <- file.path(review, "biomarker_partial_correlations_phyto202509.csv")
if (!file.exists(pc_file)) {
  stop("Run biomarker_phyto202509.R first (missing ", pc_file, ").", call. = FALSE)
}
pc <- read.csv(pc_file, stringsAsFactors = FALSE)

# Numeric x so groups can be spaced; facet_grid 1-column strips were clipping
bio_meta <- data.frame(
  biomarker_label = c("BMI", "HbA1c", "C-peptide", "TAG/HDL-C", "CRP",
                      "Adiponectin", "Leptin", "Total chol", "LDL"),
  bio_short = c("BMI", "HbA1c", "C-peptide", "TG/HDL", "CRP",
                "Adiponectin", "Leptin", "Tot chol", "LDL"),
  group = c("Adiposity", "Glycaemia & IR", "Glycaemia & IR", "Glycaemia & IR",
            "Inflammation", "Adipokines", "Adipokines", "Lipids", "Lipids"),
  x = c(1, 2.5, 3.5, 4.5, 6.0, 7.5, 8.5, 10.0, 11.0),
  stringsAsFactors = FALSE
)
x_lim <- c(min(bio_meta$x) - 0.55, max(bio_meta$x) + 0.55)
x_exp <- expansion(mult = 0)

heat_dat <- pc %>%
  inner_join(bio_meta, by = "biomarker_label") %>%
  mutate(
    sig_f = factor(signature, levels = c("ratio", "animal", "plant")),
    lab = sprintf("%.2f", partial_r)
  )

lim <- ceiling(max(abs(heat_dat$partial_r), na.rm = TRUE) * 20) / 20
if (!is.finite(lim) || lim < 0.2) lim <- 0.25

grp_df <- bio_meta %>%
  group_by(group) %>%
  summarise(xmin = min(x) - 0.48, xmax = max(x) + 0.48, .groups = "drop") %>%
  mutate(
    xmid = (xmin + xmax) / 2,
    group = factor(group, levels = c("Adiposity", "Glycaemia & IR", "Inflammation",
                                     "Adipokines", "Lipids"))
  )

heat_theme_y <- theme(
  axis.line = element_blank(),
  axis.ticks.y = element_blank(),
  axis.title = element_blank(),
  axis.text.y = element_text(size = fig_font, colour = "black", lineheight = 0.95)
)

p_head <- ggplot(grp_df) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0.12, ymax = 0.88),
            fill = "grey82", colour = "white", linewidth = 0.3) +
  geom_text(aes(x = xmid, y = 0.5, label = group), size = 6 / .pt) +
  scale_x_continuous(limits = x_lim, expand = x_exp) +
  scale_y_continuous(breaks = 0.5, labels = y_lab_width, limits = c(0, 1)) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = fig_font) +
  heat_theme_y +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = fig_font, colour = NA, lineheight = 0.95),
    plot.margin = margin(2, 4, 0, 2, "mm")
  )

p_heat <- ggplot(heat_dat, aes(x = x, y = sig_f, fill = partial_r)) +
  geom_tile(width = 0.92, height = 0.9, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = lab), size = 6.5 / .pt, colour = "black") +
  scale_x_continuous(breaks = bio_meta$x, labels = bio_meta$bio_short,
                     limits = x_lim, expand = x_exp) +
  scale_y_discrete(labels = sig_row, expand = expansion(mult = 0.04)) +
  scale_fill_gradientn(
    colours = heat_colours, limits = c(-lim, lim), oob = squish,
    name = "Partial r"
  ) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(title.position = "left", title.vjust = 1,
                                barwidth = unit(45, "mm"), barheight = unit(2.4, "mm"),
                                frame.colour = NA, ticks.colour = "white")) +
  theme_classic(base_size = fig_font) +
  heat_theme_y +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = fig_font, colour = "black",
                               angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = fig_font),
    legend.text = element_text(size = fig_font - 1),
    plot.margin = margin(1, 4, 2, 2, "mm")
  )

p_b <- (p_head / p_heat) +
  plot_layout(heights = c(0.16, 1)) +
  plot_annotation(
    title = "Partial correlations with biomarkers",
    theme = theme(plot.title = element_text(size = fig_font, face = "bold", hjust = 0))
  )

# =============================================================================
# Combine and save (stacked so both panels get the full 180 mm width)
# =============================================================================
fig2 <- (p_a / p_b) +
  plot_layout(heights = c(1.05, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 9, face = "bold"),
        plot.margin = margin(3, 3, 3, 3, "mm"))

out_png <- file.path(review, "Figure2_signature_intake_biomarker.png")
out_pdf <- file.path(review, "Figure2_signature_intake_biomarker.pdf")
fig_w <- 180
fig_h <- 155
if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(out_png, fig2, width = fig_w, height = fig_h, units = "mm", dpi = 300, device = ragg::agg_png)
} else {
  ggsave(out_png, fig2, width = fig_w, height = fig_h, units = "mm", dpi = 300, type = "cairo")
}
if (capabilities("cairo")) {
  ggsave(out_pdf, fig2, width = fig_w, height = fig_h, units = "mm", device = grDevices::cairo_pdf)
} else {
  ggsave(out_pdf, fig2, width = fig_w, height = fig_h, units = "mm")
}
cat("Saved:\n ", out_png, "\n ", out_pdf, "\n", sep = "")
