library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)

##############
### Animal protein
#############

load("/udd/nhtes/git/protein-metab-mortality/WHI-replication/aprot_predicted_score_WHI.RData")

#### Format for main plot by cohorts including WHI to be ran in /udd/nhtes/git/lifestyle-metab-mortality/Enet/Correlation_plot_signature_vs_lifestylescore_OPPOSITE.R
signature$aprot_p_sig <- signature$sig
signature$newid <- signature$ID
signature$cohort <- "WHI"
merge_WHI <- signature[,c("newid", "aprot_avg_p", "cohort", "aprot_p_sig")]

merge_harvard <- read.csv("/udd/nhtes/git/protein-metab-mortality/Enet/Animal/signature_aprot_0606-PR.csv")
merge_harvard <- merge_harvard[,c("newid", "aprot_avg_p", "cohort", "aprot_p_sig")]

### Merge both datasets

merge_allcohorts <- rbind(merge_harvard, merge_WHI)

# format for by cohort plot
merge_allcohorts=merge_allcohorts %>% mutate(cohort=case_when(cohort=="nhs1" ~ "NHS",
                                                cohort=="nhs2" ~ "NHSII",
                                                cohort=="hpfs" ~ "HPFS",
                                                cohort=="WHI" ~ "WHI" ))

#format for lifestyle factors vs signature
merge_allcohorts <- merge_allcohorts %>% gather(type, actual_score_aprot, c("aprot_avg_p")) %>% #actual score
  gather(type2, sig_aprot, c("aprot_p_sig")) %>% #training and testing signature
  mutate(type2=case_when(type2=="aprot_p_sig" ~ "aprot_avg_p")) %>%
  filter(type==type2) %>%
  dplyr::select(-type2)

#### By cohort

NHS_group_aprot <- merge_allcohorts[which(merge_allcohorts$cohort=="NHS"),]
NHSII_group_aprot <- merge_allcohorts[which(merge_allcohorts$cohort=="NHSII"),]
HPFS_group_aprot <- merge_allcohorts[which(merge_allcohorts$cohort=="HPFS"),]
WHI_group_aprot <- merge_allcohorts[which(merge_allcohorts$cohort=="WHI"),]

cor.test(NHS_group_aprot$sig_aprot,   NHS_group_aprot$actual_score_aprot) 
cor.test(NHSII_group_aprot$sig_aprot, NHSII_group_aprot$actual_score_aprot)
cor.test(HPFS_group_aprot$sig_aprot,  HPFS_group_aprot$actual_score_aprot) 
cor.test(WHI_group_aprot$sig_aprot,   WHI_group_aprot$actual_score_aprot) 


#### test 
plot_NHS_aprot <- ggplot(NHS_group_aprot, aes(x = actual_score_aprot, y = sig_aprot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  # ggtitle("NHS Cohort") +
  xlab("Animal protein intake") +
  ylab("Animal protein signature") +
    geom_text(x = Inf, y = -Inf, label = "r=0.26, P < 0.0001", hjust = 1, vjust = 0, size = 5, color = "black") +
        theme(
    text = element_text(size = 16),  # Adjust the overall text size
    axis.title = element_text(size = 16),  # Adjust axis title size
    axis.text = element_text(size = 16),  # Adjust axis text size
    plot.title = element_text(size = 16),  # Adjust plot title size
    plot.subtitle = element_text(size = 14),  # Adjust plot subtitle size
    plot.caption = element_text(size = 16)  # Adjust plot caption size
  ) 


plot_NHSII_aprot <- ggplot(NHSII_group_aprot, aes(x = actual_score_aprot, y = sig_aprot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  # ggtitle("NHSII Cohort") +
  xlab("Animal protein intake") +
  ylab("Animal protein signature") +
    geom_text(x = Inf, y = -Inf, label = "r=0.28, P < 0.0001", hjust = 1, vjust = 0, size = 5, color = "black") +
        theme(
    text = element_text(size = 16),  # Adjust the overall text size
    axis.title = element_text(size = 16),  # Adjust axis title size
    axis.text = element_text(size = 16),  # Adjust axis text size
    plot.title = element_text(size = 16),  # Adjust plot title size
    plot.subtitle = element_text(size = 14),  # Adjust plot subtitle size
    plot.caption = element_text(size = 16)  # Adjust plot caption size
  ) 


plot_HPFS_aprot <- ggplot(HPFS_group_aprot, aes(x = actual_score_aprot, y = sig_aprot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  # ggtitle("HPFS Cohort") +
  xlab("Animal protein intake") +
  ylab("Animal protein signature") +
    geom_text(x = Inf, y = -Inf, label = "r=0.34, P < 0.0001", hjust = 1, vjust = 0, size = 5, color = "black") +
        theme(
    text = element_text(size = 16),  # Adjust the overall text size
    axis.title = element_text(size = 16),  # Adjust axis title size
    axis.text = element_text(size = 16),  # Adjust axis text size
    plot.title = element_text(size = 16),  # Adjust plot title size
    plot.subtitle = element_text(size = 14),  # Adjust plot subtitle size
    plot.caption = element_text(size = 16)  # Adjust plot caption size
  ) 


plot_WHI_aprot <- ggplot(WHI_group_aprot, aes(x = actual_score_aprot, y = sig_aprot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  # ggtitle("WHI Cohort") +
  xlab("Animal protein intake") +
  ylab("Animal protein signature") +
    geom_text(x = Inf, y = -Inf, label = "r=0.23, P < 0.0001", hjust = 1, vjust = 0, size = 5, color = "black") +
        theme(
    text = element_text(size = 16),  # Adjust the overall text size
    axis.title = element_text(size = 16),  # Adjust axis title size
    axis.text = element_text(size = 16),  # Adjust axis text size
    plot.title = element_text(size = 16),  # Adjust plot title size
    plot.subtitle = element_text(size = 14),  # Adjust plot subtitle size
    plot.caption = element_text(size = 16)  # Adjust plot caption size
  )  

grid1 <- grid.arrange(plot_NHS_aprot, plot_NHSII_aprot, plot_HPFS_aprot, plot_WHI_aprot, ncol = 1)

ggsave("/udd/nhtes/git/protein-metab-mortality/Enet/Correlation/cohort_scatter_plots_aprot.png", 
       plot = grid1, 
       width = 7, height = 25, dpi = 300)



###############
### Plant protein
###############

load("/udd/nhtes/git/protein-metab-mortality/WHI-replication/vprot_predicted_score_WHI.RData")

#### Format for main plot by cohorts including WHI to be ran in /udd/nhtes/git/lifestyle-metab-mortality/Enet/Correlation_plot_signature_vs_lifestylescore_OPPOSITE.R
signature_vprot$vprot_p_sig <- signature_vprot$sig
signature_vprot$newid <- signature_vprot$ID
signature_vprot$cohort <- "WHI"
merge_WHI <- signature_vprot[,c("newid", "vprot_avg_p", "cohort", "vprot_p_sig")]

merge_harvard <- read.csv("/udd/nhtes/git/protein-metab-mortality/Enet/Plant/signature_vprot_0606-PR.csv")
merge_harvard <- merge_harvard[,c("newid", "vprot_avg_p", "cohort", "vprot_p_sig")]

### Merge both datasets

merge_allcohorts <- rbind(merge_harvard, merge_WHI)

# format for by cohort plot
merge_allcohorts=merge_allcohorts %>% mutate(cohort=case_when(cohort=="nhs1" ~ "NHS",
                                                cohort=="nhs2" ~ "NHSII",
                                                cohort=="hpfs" ~ "HPFS",
                                                cohort=="WHI" ~ "WHI" ))

#format for lifestyle factors vs signature
merge_allcohorts <- merge_allcohorts %>% gather(type, actual_score_vprot, c("vprot_avg_p")) %>% #actual score
  gather(type2, sig_vprot, c("vprot_p_sig")) %>% #training and testing signature
  mutate(type2=case_when(type2=="vprot_p_sig" ~ "vprot_avg_p")) %>%
  filter(type==type2) %>%
  dplyr::select(-type2)

#### By cohort

NHS_group_vprot <- merge_allcohorts[which(merge_allcohorts$cohort=="NHS"),]
NHSII_group_vprot <- merge_allcohorts[which(merge_allcohorts$cohort=="NHSII"),]
HPFS_group_vprot <- merge_allcohorts[which(merge_allcohorts$cohort=="HPFS"),]
WHI_group_vprot <- merge_allcohorts[which(merge_allcohorts$cohort=="WHI"),]

cor.test(NHS_group_vprot$sig_vprot,NHS_group_vprot$actual_score_vprot) 
cor.test(NHSII_group_vprot$sig_vprot,NHSII_group_vprot$actual_score_vprot)
cor.test(HPFS_group_vprot$sig_vprot,HPFS_group_vprot$actual_score_vprot) 
cor.test(WHI_group_vprot$sig_vprot,WHI_group_vprot$actual_score_vprot) 


#### test 
plot_NHS_vprot <- ggplot(NHS_group_vprot, aes(x = actual_score_vprot, y = sig_vprot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  # ggtitle("NHS Cohort") +
  xlab("Plant protein intake") +
  ylab("Plant protein signature") +
    geom_text(x = Inf, y = -Inf, label = "r=0.35, P < 0.0001", hjust = 1, vjust = 0, size = 5, color = "black") +
        theme(
    text = element_text(size = 16),  # Adjust the overall text size
    axis.title = element_text(size = 16),  # Adjust axis title size
    axis.text = element_text(size = 16),  # Adjust axis text size
    plot.title = element_text(size = 16),  # Adjust plot title size
    plot.subtitle = element_text(size = 14),  # Adjust plot subtitle size
    plot.caption = element_text(size = 16)  # Adjust plot caption size
  ) 


plot_NHSII_vprot <- ggplot(NHSII_group_vprot, aes(x = actual_score_vprot, y = sig_vprot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  # ggtitle("NHSII Cohort") +
  xlab("Plant protein intake") +
  ylab("Plant protein signature") +
    geom_text(x = Inf, y = -Inf, label = "r=0.39, P < 0.0001", hjust = 1, vjust = 0, size = 5, color = "black") +
        theme(
    text = element_text(size = 16),  # Adjust the overall text size
    axis.title = element_text(size = 16),  # Adjust axis title size
    axis.text = element_text(size = 16),  # Adjust axis text size
    plot.title = element_text(size = 16),  # Adjust plot title size
    plot.subtitle = element_text(size = 14),  # Adjust plot subtitle size
    plot.caption = element_text(size = 16)  # Adjust plot caption size
  ) 


plot_HPFS_vprot <- ggplot(HPFS_group_vprot, aes(x = actual_score_vprot, y = sig_vprot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  # ggtitle("HPFS Cohort") +
  xlab("Plant protein intake") +
  ylab("Plant protein signature") +
    geom_text(x = Inf, y = -Inf, label = "r=0.48, P < 0.0001", hjust = 1, vjust = 0, size = 5, color = "black") +
        theme(
    text = element_text(size = 16),  # Adjust the overall text size
    axis.title = element_text(size = 16),  # Adjust axis title size
    axis.text = element_text(size = 16),  # Adjust axis text size
    plot.title = element_text(size = 16),  # Adjust plot title size
    plot.subtitle = element_text(size = 14),  # Adjust plot subtitle size
    plot.caption = element_text(size = 16)  # Adjust plot caption size
  ) 


plot_WHI_vprot <- ggplot(WHI_group_vprot, aes(x = actual_score_vprot, y = sig_vprot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  # ggtitle("WHI Cohort") +
  xlab("Plant protein intake") +
  ylab("Plant protein signature") +
    geom_text(x = Inf, y = -Inf, label = "r=0.31, P < 0.0001", hjust = 1, vjust = 0, size = 5, color = "black") +
        theme(
    text = element_text(size = 16),  # Adjust the overall text size
    axis.title = element_text(size = 16),  # Adjust axis title size
    axis.text = element_text(size = 16),  # Adjust axis text size
    plot.title = element_text(size = 16),  # Adjust plot title size
    plot.subtitle = element_text(size = 14),  # Adjust plot subtitle size
    plot.caption = element_text(size = 16)  # Adjust plot caption size
  ) 


grid2 <- grid.arrange(plot_NHS_vprot, plot_NHSII_vprot, plot_HPFS_vprot, plot_WHI_vprot, ncol = 1)

ggsave("/udd/nhtes/git/protein-metab-mortality/Enet/Correlation/cohort_scatter_plots_vprot.png", 
       plot = grid2, 
       width = 7, height = 25, dpi = 300)




###############
### Protein Ratio
###############

load("/udd/nhtes/git/protein-metab-mortality/WHI-replication/PR_predicted_score_WHI.RData")

#### Format for main plot by cohorts including WHI to be ran in /udd/nhtes/git/lifestyle-metab-mortality/Enet/Correlation_plot_signature_vs_lifestylescore_OPPOSITE.R
signature_pr$protein_ratio_sig <- signature_pr$sig
signature_pr$newid <- signature_pr$ID
signature_pr$cohort <- "WHI"
merge_WHI <- signature_pr[,c("newid", "protein_ratio", "cohort", "protein_ratio_sig")]

merge_harvard <- read.csv("/udd/nhtes/git/protein-metab-mortality/Enet/Protein_ratio/signature_pr_0606.csv")
merge_harvard <- merge_harvard[,c("newid", "protein_ratio", "cohort", "protein_ratio_sig")]

### Merge both datasets

merge_allcohorts <- rbind(merge_harvard, merge_WHI)

# format for by cohort plot
merge_allcohorts=merge_allcohorts %>% mutate(cohort=case_when(cohort=="nhs1" ~ "NHS",
                                                cohort=="nhs2" ~ "NHSII",
                                                cohort=="hpfs" ~ "HPFS",
                                                cohort=="WHI" ~ "WHI" ))

#format for lifestyle factors vs signature
merge_allcohorts <- merge_allcohorts %>% gather(type, actual_score_pr, c("protein_ratio")) %>% #actual score
  gather(type2, sig_pr, c("protein_ratio_sig")) %>% #training and testing signature
  mutate(type2=case_when(type2=="protein_ratio_sig" ~ "protein_ratio")) %>%
  filter(type==type2) %>%
  dplyr::select(-type2)

#### By cohort

NHS_group_pr <- merge_allcohorts[which(merge_allcohorts$cohort=="NHS"),]
NHSII_group_pr <- merge_allcohorts[which(merge_allcohorts$cohort=="NHSII"),]
HPFS_group_pr <- merge_allcohorts[which(merge_allcohorts$cohort=="HPFS"),]
WHI_group_pr <- merge_allcohorts[which(merge_allcohorts$cohort=="WHI"),]

cor.test(NHS_group_pr$sig_pr,NHS_group_pr$actual_score_pr) 
cor.test(NHSII_group_pr$sig_pr,NHSII_group_pr$actual_score_pr)
cor.test(HPFS_group_pr$sig_pr,HPFS_group_pr$actual_score_pr) 
cor.test(WHI_group_pr$sig_pr,WHI_group_pr$actual_score_pr) 


#### test 
plot_NHS_pr <- ggplot(NHS_group_pr, aes(x = actual_score_pr, y = sig_pr)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  # ggtitle("NHS Cohort") +
  xlab("Protein ratio (P:A)") +
  ylab("Protein ratio signature") +
    geom_text(x = Inf, y = -Inf, label = "r=0.28, P < 0.0001", hjust = 1, vjust = 0, size = 5, color = "black") + 
theme(
    text = element_text(size = 16),  # Adjust the overall text size
    axis.title = element_text(size = 16),  # Adjust axis title size
    axis.text = element_text(size = 16),  # Adjust axis text size
    plot.title = element_text(size = 16),  # Adjust plot title size
    plot.subtitle = element_text(size = 14),  # Adjust plot subtitle size
    plot.caption = element_text(size = 16)  # Adjust plot caption size
  )

plot_NHSII_pr <- ggplot(NHSII_group_pr, aes(x = actual_score_pr, y = sig_pr)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  # ggtitle("NHSII Cohort") +
  xlab("Protein ratio (P:A)") +
  ylab("Protein ratio signature") +
    geom_text(x = Inf, y = -Inf, label = "r=0.36, P < 0.0001", hjust = 1, vjust = 0, size = 5, color = "black") +
theme(
    text = element_text(size = 16),  # Adjust the overall text size
    axis.title = element_text(size = 16),  # Adjust axis title size
    axis.text = element_text(size = 16),  # Adjust axis text size
    plot.title = element_text(size = 16),  # Adjust plot title size
    plot.subtitle = element_text(size = 14),  # Adjust plot subtitle size
    plot.caption = element_text(size = 16)  # Adjust plot caption size
  )

plot_HPFS_pr <- ggplot(HPFS_group_pr, aes(x = actual_score_pr, y = sig_pr)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  # ggtitle("HPFS Cohort") +
  xlab("Protein ratio (P:A)") +
  ylab("Protein ratio signature") +
    geom_text(x = Inf, y = -Inf, label = "r=0.41, P < 0.0001", hjust = 1, vjust = 0, size = 5, color = "black") +
    theme(
    text = element_text(size = 16),  # Adjust the overall text size
    axis.title = element_text(size = 16),  # Adjust axis title size
    axis.text = element_text(size = 16),  # Adjust axis text size
    plot.title = element_text(size = 16),  # Adjust plot title size
    plot.subtitle = element_text(size = 14),  # Adjust plot subtitle size
    plot.caption = element_text(size = 16)  # Adjust plot caption size
  )


plot_WHI_pr <- ggplot(WHI_group_pr, aes(x = actual_score_pr, y = sig_pr)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  # ggtitle("WHI Cohort") +
  xlab("Protein ratio (P:A)") +
  ylab("Protein ratio signature") +
    geom_text(x = Inf, y = -Inf, label = "r=0.25, P < 0.0001", hjust = 1, vjust = 0, size = 5, color = "black") +
    theme(
    text = element_text(size = 16),  # Adjust the overall text size
    axis.title = element_text(size = 16),  # Adjust axis title size
    axis.text = element_text(size = 16),  # Adjust axis text size
    plot.title = element_text(size = 16),  # Adjust plot title size
    plot.subtitle = element_text(size = 14),  # Adjust plot subtitle size
    plot.caption = element_text(size = 16)  # Adjust plot caption size
  ) 


grid3 <- grid.arrange(plot_NHS_pr, plot_NHSII_pr, plot_HPFS_pr, plot_WHI_pr, ncol = 1)

ggsave("/udd/nhtes/git/protein-metab-mortality/Enet/Correlation/cohort_scatter_plots_pr.png", 
       plot = grid3, 
       width = 7, height = 25, dpi = 300)


grid_final <- grid.arrange(grid2, grid1, grid3, ncol = 3, heights = c(4, 1))
ggsave("/udd/nhtes/git/protein-metab-mortality/Enet/Correlation/cohort_scatter_plots_all.png", 
       plot = grid_final, 
       width = 19, height = 25, dpi = 300)


ggsave("/udd/nhtes/git/protein-metab-mortality/Enet/Correlation/cohort_scatter_plots_all.svg", 
       plot = grid_final, 
       width = 19, height = 25, dpi = 300)






#### Extract all correlations into a clean table ####

get_r  <- function(x, y) round(cor.test(x, y, use = "complete.obs")$estimate, 2)
get_p  <- function(x, y) {
  p <- cor.test(x, y, use = "complete.obs")$p.value
  ifelse(p < 0.0001, "< 0.0001", round(p, 4))
}

corr_results <- data.frame(
  cohort    = c("NHS", "NHSII", "HPFS", "WHI"),
  plant_r   = c(get_r(NHS_group_vprot$actual_score_vprot,   NHS_group_vprot$sig_vprot),
                get_r(NHSII_group_vprot$actual_score_vprot, NHSII_group_vprot$sig_vprot),
                get_r(HPFS_group_vprot$actual_score_vprot,  HPFS_group_vprot$sig_vprot),
                get_r(WHI_group_vprot$actual_score_vprot,   WHI_group_vprot$sig_vprot)),
  plant_p   = c(get_p(NHS_group_vprot$actual_score_vprot,   NHS_group_vprot$sig_vprot),
                get_p(NHSII_group_vprot$actual_score_vprot, NHSII_group_vprot$sig_vprot),
                get_p(HPFS_group_vprot$actual_score_vprot,  HPFS_group_vprot$sig_vprot),
                get_p(WHI_group_vprot$actual_score_vprot,   WHI_group_vprot$sig_vprot)),
  animal_r  = c(get_r(NHS_group_aprot$actual_score_aprot,   NHS_group_aprot$sig_aprot),
                get_r(NHSII_group_aprot$actual_score_aprot, NHSII_group_aprot$sig_aprot),
                get_r(HPFS_group_aprot$actual_score_aprot,  HPFS_group_aprot$sig_aprot),
                get_r(WHI_group_aprot$actual_score_aprot,   WHI_group_aprot$sig_aprot)),
  animal_p  = c(get_p(NHS_group_aprot$actual_score_aprot,   NHS_group_aprot$sig_aprot),
                get_p(NHSII_group_aprot$actual_score_aprot, NHSII_group_aprot$sig_aprot),
                get_p(HPFS_group_aprot$actual_score_aprot,  HPFS_group_aprot$sig_aprot),
                get_p(WHI_group_aprot$actual_score_aprot,   WHI_group_aprot$sig_aprot)),
  ratio_r   = c(get_r(NHS_group_pr$actual_score_pr,   NHS_group_pr$sig_pr),
                get_r(NHSII_group_pr$actual_score_pr, NHSII_group_pr$sig_pr),
                get_r(HPFS_group_pr$actual_score_pr,  HPFS_group_pr$sig_pr),
                get_r(WHI_group_pr$actual_score_pr,   WHI_group_pr$sig_pr)),
  ratio_p   = c(get_p(NHS_group_pr$actual_score_pr,   NHS_group_pr$sig_pr),
                get_p(NHSII_group_pr$actual_score_pr, NHSII_group_pr$sig_pr),
                get_p(HPFS_group_pr$actual_score_pr,  HPFS_group_pr$sig_pr),
                get_p(WHI_group_pr$actual_score_pr,   WHI_group_pr$sig_pr))
)

print(corr_results)

write.csv(corr_results,
          "/udd/nhtes/git/protein-metab-mortality/Results_correlations_Fig2.csv",
          row.names = FALSE)
