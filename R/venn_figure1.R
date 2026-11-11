source("00_paths.R")
library(dplyr)
library(readxl)
library(ggplot2)
library(grid)
library(gtable)
library(gridExtra)
library(lattice)
library(tidyr)
library(scales)
library(extrafont)
library(Cairo)
library(svglite)

### venn diagram https://www.rdocumentation.org/packages/VennDiagram/versions/1.7.3/topics/calculate.overlap

# readin metabolites names
met_names <- read.csv(met_names_csv, header=TRUE, sep=',')

#readin metabolites dataset
load(file.path(repo, "Enet/Animal/Prediction_results_aprot-PR-multi-ep.RData"))

load(file.path(repo, "Enet/Plant/Prediction_results_vprot-PR-multi-ep.RData"))

load(file.path(repo, "Enet/Protein_ratio/Prediction_results_PR_multi_ep.RData"))

# readin selected metabolites
aprot_enet_heatmap_data <- data.frame(Metabolite=met_aprot, aprot_coef_numeric=as.numeric(aprot_coef))

vprot_enet_heatmap_data <- data.frame(Metabolite=met_vprot, vprot_coef_numeric=as.numeric(vprot_coef))

pr_enet_heatmap_data <- data.frame(Metabolite=met_pr, pr_coef_numeric=as.numeric(pr_coef))


#readin beta and ORs of association between all metabolites and lifestyle score anf individual factors

vprot_p <- read.csv(file.path(repo, "Results_vprot_p_ztran-multi-ep.csv"), header=T) 
aprot_p <- read.csv(file.path(repo, "Results_aprot_p_ztran-multi-ep.csv"), header=T) 
prot_ratio <- read.csv(file.path(repo, "Results_prot_ratio_ztran-multi-ep.csv"), header=T) 
legu <- read.csv(file.path(repo, "Results_legu_ztran-multi-ep.csv"), header=T) 
nuts <- read.csv(file.path(repo, "Results_nuts_ztran-multi-ep.csv"), header=T) 
rgrain <- read.csv(file.path(repo, "Results_rgrain_ztran-multi-ep.csv"), header=T) 
wgrain <- read.csv(file.path(repo, "Results_wgrain_ztran-multi-ep.csv"), header=T) 
procm <- read.csv(file.path(repo, "Results_procm_ztran-multi-ep.csv"), header=T) 
rmeat <- read.csv(file.path(repo, "Results_rmeat_ztran-multi-ep.csv"), header=T) 
poul <- read.csv(file.path(repo, "Results_poul_ztran-multi-ep.csv"), header=T) 
fish <- read.csv(file.path(repo, "Results_fish_ztran-multi-ep.csv"), header=T) 
egg <- read.csv(file.path(repo, "Results_egg_ztran-multi-ep.csv"), header=T) 
dairy <- read.csv(file.path(repo, "Results_dairy_ztran-multi-ep.csv"), header=T) 

mort <- read.csv(file.path(repo, "Endpoints/Results_indmetab_dead.csv"), header=T) 
db <- read.csv(file.path(repo, "Endpoints/Results_indmetab_db2.csv"), header=T) 
cvd <- read.csv(file.path(repo, "Endpoints/Results_indmetab_cvd.csv"), header=T) 
longev <- read.csv(file.path(repo, "Endpoints/Results_indmetab_longev.csv"), header=T) 


# vprot_p <- read.csv(file.path(repo, "Results_vprot_p_ztran.csv"), header=T) 
# aprot_p <- read.csv(file.path(repo, "Results_aprot_p_ztran.csv"), header=T) 
# prot_ratio <- read.csv(file.path(repo, "Results_prot_ratio_ztran.csv"), header=T) 
# legu <- read.csv(file.path(repo, "Results_legu_ztran.csv"), header=T) 
# nuts <- read.csv(file.path(repo, "Results_nuts_ztran.csv"), header=T) 
# rgrain <- read.csv(file.path(repo, "Results_rgrain_ztran.csv"), header=T) 
# wgrain <- read.csv(file.path(repo, "Results_wgrain_ztran.csv"), header=T) 
# procm <- read.csv(file.path(repo, "Results_procm_ztran.csv"), header=T) 
# rmeat <- read.csv(file.path(repo, "Results_rmeat_ztran.csv"), header=T) 
# poul <- read.csv(file.path(repo, "Results_poul_ztran.csv"), header=T) 
# fish <- read.csv(file.path(repo, "Results_fish_ztran.csv"), header=T) 
# egg <- read.csv(file.path(repo, "Results_egg_ztran.csv"), header=T) 
# dairy <- read.csv(file.path(repo, "Results_dairy_ztran.csv"), header=T) 

## merge the selected metabolites and their coefficients from enet with beta for sources of protein intake and food groups

data_heatmap1 <- merge(vprot_enet_heatmap_data, aprot_enet_heatmap_data, by.x="Metabolite", by.y ="Metabolite", all=T) # keep all metabolites from both dataframe. If it puts missing and its problematic, replace them with zeros

data_heatmap2 <- merge(data_heatmap1, pr_enet_heatmap_data, by.x="Metabolite", by.y ="Metabolite", all=T) # keep all metabolites from both dataframe. If it puts missing and its problematic, replace them with zeros

data_heatmap_final <- merge(data_heatmap2, met_names[,c("HMDB", "rep", "name", "cate")], by.x="Metabolite", by.y ="HMDB", all.x=T)

#write excel file with overlaps between Enet-selected metabolites for aprot, vprot and PR
write.csv(data_heatmap_final, file.path(repo, "Heatmap/protein_intake_overlap-multi-ep.csv"), na="")


data_heatmap3 <- merge(data_heatmap2, vprot_p[,c("Metabolite","Beta.mod2", "p.value.mod2", "p.value.mod2.bon")], by.x="Metabolite", by.y ="Metabolite", all.x=T)
colnames(data_heatmap3) 
names(data_heatmap3)[names(data_heatmap3) == "Beta.mod2"] <- "vprot_coeff2"
names(data_heatmap3)[names(data_heatmap3) == "p.value.mod2"] <- "vprot_pval2"
names(data_heatmap3)[names(data_heatmap3) == "p.value.mod2.bon"] <- "vprot_pval2bon"

data_heatmap4 <- merge(data_heatmap3, aprot_p[,c("Metabolite","Beta.mod2", "p.value.mod2", "p.value.mod2.bon")], by.x="Metabolite", by.y ="Metabolite", all.x=T)
colnames(data_heatmap4) 
names(data_heatmap4)[names(data_heatmap4) == "Beta.mod2"] <- "aprot_coeff2"
names(data_heatmap4)[names(data_heatmap4) == "p.value.mod2"] <- "aprot_pval2"
names(data_heatmap4)[names(data_heatmap4) == "p.value.mod2.bon"] <- "aprot_pval2bon"

data_heatmap5 <- merge(data_heatmap4, prot_ratio[,c("Metabolite","Beta.mod2", "p.value.mod2", "p.value.mod2.bon")], by.x="Metabolite", by.y ="Metabolite", all.x=T)
colnames(data_heatmap5) 
names(data_heatmap5)[names(data_heatmap5) == "Beta.mod2"] <- "prot_ratio_coeff2"
names(data_heatmap5)[names(data_heatmap5) == "p.value.mod2"] <- "prot_ratio_pval2"
names(data_heatmap5)[names(data_heatmap5) == "p.value.mod2.bon"] <- "prot_ratio_pval2bon"

data_heatmap6 <- merge(data_heatmap5, legu[,c("Metabolite","Beta.mod2", "p.value.mod2", "p.value.mod2.bon")], by.x="Metabolite", by.y ="Metabolite", all.x=T)
colnames(data_heatmap6) 
names(data_heatmap6)[names(data_heatmap6) == "Beta.mod2"] <- "legu_coeff2"
names(data_heatmap6)[names(data_heatmap6) == "p.value.mod2"] <- "legu_pval2"
names(data_heatmap6)[names(data_heatmap6) == "p.value.mod2.bon"] <- "legu_pval2bon"

data_heatmap7 <- merge(data_heatmap6, nuts[,c("Metabolite","Beta.mod2", "p.value.mod2", "p.value.mod2.bon")], by.x="Metabolite", by.y ="Metabolite", all.x=T)
colnames(data_heatmap7) 
names(data_heatmap7)[names(data_heatmap7) == "Beta.mod2"] <- "nuts_coeff2"
names(data_heatmap7)[names(data_heatmap7) == "p.value.mod2"] <- "nuts_pval2"
names(data_heatmap7)[names(data_heatmap7) == "p.value.mod2.bon"] <- "nuts_pval2bon"

data_heatmap8 <- merge(data_heatmap7, rgrain[,c("Metabolite","Beta.mod2", "p.value.mod2", "p.value.mod2.bon")], by.x="Metabolite", by.y ="Metabolite", all.x=T)
colnames(data_heatmap8) 
names(data_heatmap8)[names(data_heatmap8) == "Beta.mod2"] <- "rgrain_coeff2"
names(data_heatmap8)[names(data_heatmap8) == "p.value.mod2"] <- "rgrain_pval2"
names(data_heatmap8)[names(data_heatmap8) == "p.value.mod2.bon"] <- "rgrain_pval2bon"

data_heatmap9 <- merge(data_heatmap8, wgrain[,c("Metabolite","Beta.mod2", "p.value.mod2", "p.value.mod2.bon")], by.x="Metabolite", by.y ="Metabolite", all.x=T)
colnames(data_heatmap9) 
names(data_heatmap9)[names(data_heatmap9) == "Beta.mod2"] <- "wgrain_coeff2"
names(data_heatmap9)[names(data_heatmap9) == "p.value.mod2"] <- "wgrain_pval2"
names(data_heatmap9)[names(data_heatmap9) == "p.value.mod2.bon"] <- "wgrain_pval2bon"

data_heatmap10 <- merge(data_heatmap9, procm[,c("Metabolite","Beta.mod2", "p.value.mod2", "p.value.mod2.bon")], by.x="Metabolite", by.y ="Metabolite", all.x=T)
colnames(data_heatmap10) 
names(data_heatmap10)[names(data_heatmap10) == "Beta.mod2"] <- "procm_coeff2"
names(data_heatmap10)[names(data_heatmap10) == "p.value.mod2"] <- "procm_pval2"
names(data_heatmap10)[names(data_heatmap10) == "p.value.mod2.bon"] <- "procm_pval2bon"

data_heatmap11 <- merge(data_heatmap10, rmeat[,c("Metabolite","Beta.mod2", "p.value.mod2", "p.value.mod2.bon")], by.x="Metabolite", by.y ="Metabolite", all.x=T)
colnames(data_heatmap11) 
names(data_heatmap11)[names(data_heatmap11) == "Beta.mod2"] <- "rmeat_coeff2"
names(data_heatmap11)[names(data_heatmap11) == "p.value.mod2"] <- "rmeat_pval2"
names(data_heatmap11)[names(data_heatmap11) == "p.value.mod2.bon"] <- "rmeat_pval2bon"

data_heatmap12 <- merge(data_heatmap11, poul[,c("Metabolite","Beta.mod2", "p.value.mod2", "p.value.mod2.bon")], by.x="Metabolite", by.y ="Metabolite", all.x=T)
colnames(data_heatmap12) 
names(data_heatmap12)[names(data_heatmap12) == "Beta.mod2"] <- "poul_coeff2"
names(data_heatmap12)[names(data_heatmap12) == "p.value.mod2"] <- "poul_pval2"
names(data_heatmap12)[names(data_heatmap12) == "p.value.mod2.bon"] <- "poul_pval2bon"

data_heatmap13 <- merge(data_heatmap12, fish[,c("Metabolite","Beta.mod2", "p.value.mod2", "p.value.mod2.bon")], by.x="Metabolite", by.y ="Metabolite", all.x=T)
colnames(data_heatmap13) 
names(data_heatmap13)[names(data_heatmap13) == "Beta.mod2"] <- "fish_coeff2"
names(data_heatmap13)[names(data_heatmap13) == "p.value.mod2"] <- "fish_pval2"
names(data_heatmap13)[names(data_heatmap13) == "p.value.mod2.bon"] <- "fish_pval2bon"

data_heatmap14 <- merge(data_heatmap13, egg[,c("Metabolite","Beta.mod2", "p.value.mod2", "p.value.mod2.bon")], by.x="Metabolite", by.y ="Metabolite", all.x=T)
colnames(data_heatmap14) 
names(data_heatmap14)[names(data_heatmap14) == "Beta.mod2"] <- "egg_coeff2"
names(data_heatmap14)[names(data_heatmap14) == "p.value.mod2"] <- "egg_pval2"
names(data_heatmap14)[names(data_heatmap14) == "p.value.mod2.bon"] <- "egg_pval2bon"

data_heatmap15 <- merge(data_heatmap14, dairy[,c("Metabolite","Beta.mod2", "p.value.mod2", "p.value.mod2.bon")], by.x="Metabolite", by.y ="Metabolite", all.x=T)
colnames(data_heatmap15) 
names(data_heatmap15)[names(data_heatmap15) == "Beta.mod2"] <- "dairy_coeff2"
names(data_heatmap15)[names(data_heatmap15) == "p.value.mod2"] <- "dairy_pval2"
names(data_heatmap15)[names(data_heatmap15) == "p.value.mod2.bon"] <- "dairy_pval2bon"


#Add endpoints to the heatmap
data_heatmap16 <- merge(data_heatmap15, mort[,c("Metabolite","HR.mod2", "p.value.mod2", "p.value.mod2.bon")], by.x="Metabolite", by.y ="Metabolite", all.x=T)
colnames(data_heatmap16) 
names(data_heatmap16)[names(data_heatmap16) == "HR.mod2"] <- "mort_coeff2"
names(data_heatmap16)[names(data_heatmap16) == "p.value.mod2"] <- "mort_pval2"
names(data_heatmap16)[names(data_heatmap16) == "p.value.mod2.bon"] <- "mort_pval2bon"

data_heatmap17 <- merge(data_heatmap16, longev[,c("Metabolite","HR.mod2", "p.value.mod2", "p.value.mod2.bon")], by.x="Metabolite", by.y ="Metabolite", all.x=T)
colnames(data_heatmap17) 
names(data_heatmap17)[names(data_heatmap17) == "HR.mod2"] <- "longev_coeff2"
names(data_heatmap17)[names(data_heatmap17) == "p.value.mod2"] <- "longev_pval2"
names(data_heatmap17)[names(data_heatmap17) == "p.value.mod2.bon"] <- "longev_pval2bon"

data_heatmap18 <- merge(data_heatmap17, db[,c("Metabolite","HR.mod2", "p.value.mod2", "p.value.mod2.bon")], by.x="Metabolite", by.y ="Metabolite", all.x=T)
colnames(data_heatmap18) 
names(data_heatmap18)[names(data_heatmap18) == "HR.mod2"] <- "db_coeff2"
names(data_heatmap18)[names(data_heatmap18) == "p.value.mod2"] <- "db_pval2"
names(data_heatmap18)[names(data_heatmap18) == "p.value.mod2.bon"] <- "db_pval2bon"


data_heatmap19 <- merge(data_heatmap18, cvd[,c("Metabolite","HR.mod2", "p.value.mod2", "p.value.mod2.bon")], by.x="Metabolite", by.y ="Metabolite", all.x=T)
colnames(data_heatmap19) 
names(data_heatmap19)[names(data_heatmap19) == "HR.mod2"] <- "cvd_coeff2"
names(data_heatmap19)[names(data_heatmap19) == "p.value.mod2"] <- "cvd_pval2"
names(data_heatmap19)[names(data_heatmap19) == "p.value.mod2.bon"] <- "cvd_pval2bon"



#### Merge with met names

data_heatmap_final <- merge(data_heatmap19, met_names[,c("HMDB", "rep", "name", "cate")], by.x="Metabolite", by.y ="HMDB", all.x=T)

# write.table(data_heatmap_final, file = "Results_heatmap.csv", sep=",", quote=F, row.names=F, col.names=T)
# current file named Results_heatmap.csv has the CATE that had comma and reads as several columns instead of one single column. 
# If you re-run this code it will be ok because I excluded cate and CATE

save.image(file=file.path(repo, "Heatmap/heatmap_multi-ep.RData"))


# vprot_p <- read.csv(file.path(repo, "Results_vprot_p.csv"), header=T) 
# aprot_p <- read.csv(file.path(repo, "Results_aprot_p.csv"), header=T) 
# prot_ratio <- read.csv(file.path(repo, "Results_prot_ratio.csv"), header=T) 
# legu <- read.csv(file.path(repo, "Results_legu.csv"), header=T) 
# nuts <- read.csv(file.path(repo, "Results_nuts.csv"), header=T) 
# rgrain <- read.csv(file.path(repo, "Results_rgrain.csv"), header=T) 
# wgrain <- read.csv(file.path(repo, "Results_wgrain.csv"), header=T) 
# procm <- read.csv(file.path(repo, "Results_procm.csv"), header=T) 
# rmeat <- read.csv(file.path(repo, "Results_rmeat.csv"), header=T) 
# poul <- read.csv(file.path(repo, "Results_poul.csv"), header=T) 
# fish <- read.csv(file.path(repo, "Results_fish.csv"), header=T) 
# egg <- read.csv(file.path(repo, "Results_egg.csv"), header=T) 
# dairy <- read.csv(file.path(repo, "Results_dairy.csv"), header=T) 


#### Heatmap using GGPLOT ####

# pval_heatmap <- data_heatmap_final[,c(1, 2, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16, 18, 19, 21, 22, 24, 25, 27, 28, 30, 31, 33, 34, 36, 37, 39, 40, 42, 43, 45)]
# pval_heatmap <- data_heatmap_final[,c(1, 2, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16, 18, 19, 21, 22, 24, 25, 27, 28, 30, 31, 33, 34, 36, 37, 39, 41)]
pval_heatmap <- data_heatmap_final[,c(1:55, 57)]


# put 1 to the non-corrected pvalues below 0.05
pval_heatmap$vprot_star_pval=ifelse(pval_heatmap$vprot_pval2<=0.05,1,0)
pval_heatmap$aprot_star_pval=ifelse(pval_heatmap$aprot_pval2<=0.05,1,0)
pval_heatmap$prot_ratio_star_pval=ifelse(pval_heatmap$prot_ratio_pval2<=0.05,1,0)
pval_heatmap$legu_star_pval=ifelse(pval_heatmap$legu_pval2<=0.05,1,0)
pval_heatmap$nuts_star_pval=ifelse(pval_heatmap$nuts_pval2<=0.05,1,0)
pval_heatmap$rgrain_star_pval=ifelse(pval_heatmap$rgrain_pval2<=0.05,1,0)
pval_heatmap$wgrain_star_pval=ifelse(pval_heatmap$wgrain_pval2<=0.05,1,0)
pval_heatmap$procm_star_pval=ifelse(pval_heatmap$procm_pval2<=0.05,1,0)
pval_heatmap$rmeat_star_pval=ifelse(pval_heatmap$rmeat_pval2<=0.05,1,0)
pval_heatmap$poul_star_pval=ifelse(pval_heatmap$poul_pval2<=0.05,1,0)
pval_heatmap$fish_star_pval=ifelse(pval_heatmap$fish_pval2<=0.05,1,0)
pval_heatmap$egg_star_pval=ifelse(pval_heatmap$egg_pval2<=0.05,1,0)
pval_heatmap$dairy_star_pval=ifelse(pval_heatmap$dairy_pval2<=0.05,1,0)

pval_heatmap$mort_star_pval=ifelse(pval_heatmap$mort_pval2<=0.05,1,0)
pval_heatmap$longev_star_pval=ifelse(pval_heatmap$longev_pval2<=0.05,1,0)
pval_heatmap$db_star_pval=ifelse(pval_heatmap$db_pval2<=0.05,1,0)
pval_heatmap$cvd_star_pval=ifelse(pval_heatmap$cvd_pval2<=0.05,1,0)


# put 1 to the Bonferroni corrected pvalues below 0.05
pval_heatmap$vprot_star_bon=ifelse(pval_heatmap$vprot_pval2bon<=0.05,1,0)
pval_heatmap$aprot_star_bon=ifelse(pval_heatmap$aprot_pval2bon<=0.05,1,0)
pval_heatmap$prot_ratio_star_bon=ifelse(pval_heatmap$prot_ratio_pval2bon<=0.05,1,0)
pval_heatmap$legu_star_bon=ifelse(pval_heatmap$legu_pval2bon<=0.05,1,0)
pval_heatmap$nuts_star_bon=ifelse(pval_heatmap$nuts_pval2bon<=0.05,1,0)
pval_heatmap$rgrain_star_bon=ifelse(pval_heatmap$rgrain_pval2bon<=0.05,1,0)
pval_heatmap$wgrain_star_bon=ifelse(pval_heatmap$wgrain_pval2bon<=0.05,1,0)
pval_heatmap$procm_star_bon=ifelse(pval_heatmap$procm_pval2bon<=0.05,1,0)
pval_heatmap$rmeat_star_bon=ifelse(pval_heatmap$rmeat_pval2bon<=0.05,1,0)
pval_heatmap$poul_star_bon=ifelse(pval_heatmap$poul_pval2bon<=0.05,1,0)
pval_heatmap$fish_star_bon=ifelse(pval_heatmap$fish_pval2bon<=0.05,1,0)
pval_heatmap$egg_star_bon=ifelse(pval_heatmap$egg_pval2bon<=0.05,1,0)
pval_heatmap$dairy_star_bon=ifelse(pval_heatmap$dairy_pval2bon<=0.05,1,0)

pval_heatmap$mort_star_bon=ifelse(pval_heatmap$mort_pval2bon<=0.05,1,0)
pval_heatmap$longev_star_bon=ifelse(pval_heatmap$longev_pval2bon<=0.05,1,0)
pval_heatmap$db_star_bon=ifelse(pval_heatmap$db_pval2bon<=0.05,1,0)
pval_heatmap$cvd_star_bon=ifelse(pval_heatmap$cvd_pval2bon<=0.05,1,0)


# add the two together to be able to identify those who have both pvalues significant (remained significant even after Bonferroni correction = **)
pval_heatmap$vprot_star=pval_heatmap$vprot_star_bon+pval_heatmap$vprot_star_pval
pval_heatmap$aprot_star=pval_heatmap$aprot_star_bon+pval_heatmap$aprot_star_pval
pval_heatmap$prot_ratio_star=pval_heatmap$prot_ratio_star_bon+pval_heatmap$prot_ratio_star_pval
pval_heatmap$legu_star=pval_heatmap$legu_star_bon+pval_heatmap$legu_star_pval
pval_heatmap$nuts_star=pval_heatmap$nuts_star_bon+pval_heatmap$nuts_star_pval
pval_heatmap$rgrain_star=pval_heatmap$rgrain_star_bon+pval_heatmap$rgrain_star_pval
pval_heatmap$wgrain_star=pval_heatmap$wgrain_star_bon+pval_heatmap$wgrain_star_pval
pval_heatmap$procm_star=pval_heatmap$procm_star_bon+pval_heatmap$procm_star_pval
pval_heatmap$rmeat_star=pval_heatmap$rmeat_star_bon+pval_heatmap$rmeat_star_pval
pval_heatmap$poul_star=pval_heatmap$poul_star_bon+pval_heatmap$poul_star_pval
pval_heatmap$fish_star=pval_heatmap$fish_star_bon+pval_heatmap$fish_star_pval
pval_heatmap$egg_star=pval_heatmap$egg_star_bon+pval_heatmap$egg_star_pval
pval_heatmap$dairy_star=pval_heatmap$dairy_star_bon+pval_heatmap$dairy_star_pval

pval_heatmap$mort_star=pval_heatmap$mort_star_bon+pval_heatmap$mort_star_pval
pval_heatmap$longev_star=pval_heatmap$longev_star_bon+pval_heatmap$longev_star_pval
pval_heatmap$db_star=pval_heatmap$db_star_bon+pval_heatmap$db_star_pval
pval_heatmap$cvd_star=pval_heatmap$cvd_star_bon+pval_heatmap$cvd_star_pval


pval_bon <- pval_heatmap[,c(56, 91:107)]
pval_bon$pvalwt_vprot=3000
pval_bon$pvalwt_aprot=3000
pval_bon$pvalwt_pr=3000
pval_bon <- pval_bon[,c(1, 19:21, 2:18)]
colnames(pval_bon) = c("Metabolite","PP signature", "AP signature", "PR signature", "PP", "AP", "PR", "Legumes", "Nuts", "Refined grains", "Whole grains", "Processed meat", "Red meat", "Poultry", "Fish", "Egg", "Dairy", "Total mortality", "Longevity", "Type 2 Diabetes", "CVD")


estimate <- data_heatmap_final[,c(57, 2:5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53)]
estimate$vprot_coef_numeric <- round(estimate$vprot_coef_numeric, 3)
estimate$aprot_coef_numeric <- round(estimate$aprot_coef_numeric, 3)
estimate$pr_coef_numeric <- round(estimate$pr_coef_numeric, 3)
colnames(estimate) = c("Metabolite","PP signature", "AP signature", "PR signature", "PP", "AP", "PR", "Legumes", "Nuts", "Refined grains", "Whole grains", "Processed meat", "Red meat", "Poultry", "Fish", "Egg", "Dairy", "Total mortality", "Longevity", "Type 2 Diabetes", "CVD")

Order <- order(estimate$`PP signature`)
estimate <- estimate[Order,]
pval_bon <- pval_bon[Order,]

# Extract the weights to keep them for labeling
enet_weights <- estimate[,c(1:4)]
colnames(enet_weights) = c("Metabolite","Unscaled_plant", "Unscaled_animal", "Unscaled_pr")

# Scale values to be able to apply color following the same scale
scaled_estimate <- scale(estimate[,c(2:21)]) # this is to exclude the enet weights in the scaling

full_estimate <- cbind(enet_weights, scaled_estimate)
full_estimate <- full_estimate[,c(1, 5:24)]
# colnames(full_estimate) = c("Metabolite","PP signature", "AP signature", "PR signature", "PP", "AP", "PR", "Legumes", "Nuts", "Refined grains", "Whole grains", "Processed meat", "Red meat", "Poultry", "Fish", "Egg", "Dairy", "Total mortality", "Cancer mortality", "CVD mortality", "Longevity")

# pval_long <- gather(pval_bon, Factor, Pvalue, Weight:Smoking, factor_key=TRUE)

library(reshape)
estimate_melt <- melt(full_estimate)
pval_bon_melt <- melt(pval_bon)
enet_weights_melt <- melt(enet_weights)


estimate_melt$stars = ""
estimate_melt$stars[which(pval_bon_melt$value==1)] = "*"
estimate_melt$stars[which(pval_bon_melt$value==2)] = "**"
estimate_melt$stars[which(pval_bon_melt$value>=3000)] = "" # no print of label/weight estimates
# estimate_melt$stars[which(pval_bon_melt$value>=3000)] = paste(enet_weights_melt$value)

## Create heatmap with ggplot
estimate_melt$Metabolite = factor(estimate_melt$Metabolite, levels = estimate$Metabolite)

range(estimate_melt$value, na.rm = TRUE) #-3.33 to 4.61




colours <- c("#0046FF", "#5876FF", "#7F97FF", "#A9B9FF", "#D4DBFF", "#FFFFFF", "#FFE9F3","#FFB0D6", "#FF7FBC","#FF5AA8","#FF2F92")


plot <- ggplot(data = estimate_melt, aes(x = variable, y = Metabolite, fill = value)) + #data and lables
                     geom_tile(color = "white", linewidth = 1) + 
                # facet_grid(cols=vars(Factor)) + # try to segment plot by column; change set for factor
                scale_fill_gradientn(limits=range(-4.4, 4.4), colours=colours) + #make it symmetrical on both sides for colors to match orientation of associations
                geom_text(aes(x = variable, y = Metabolite, label=stars), 
                              color="black", 
                              size=6, 
                              nudge_x = 0, 
                              nudge_y = 0) +
                theme_minimal()+ # minimal theme
                theme(text = element_text(family = "Arial", size=16)) +
                theme(axis.text.y = element_text(angle = 0, 
                                                vjust = 0.5, 
                                                size = 18, 
                     
                                                colour = "black"), 
                      axis.text.x = element_text(angle = 90, 
                                                vjust = 0.5, 
                                                hjust=0.05,
                                                size = 18, 
                                               
                                                colour = "black"),
                      axis.title = element_blank(),
                      legend.title = element_text(size=18), 
                      legend.text = element_blank(),
                      legend.key.height = unit(0.5, 'cm'),
                      legend.key.width = unit(1.5, 'cm')) +
                      scale_x_discrete(position = "top") +
                theme(strip.text.y = element_text(angle = 0, 
                                                  hjust = 0, 
                                                  vjust = 0.5, 
                                                  
                                                  size = 20), 
                                                  legend.position="bottom")+
                guides(fill=guide_colorbar(title = "Negative            Positive",
                                          title.position = "top",
                                          nrow = 2, 
                                          label.position = "bottom"))

ggsave(filename=file.path(file.path(repo, "Heatmap/Heatmap_ztrans_multi-ep.svg")), 
       plot = plot, dpi = 300, width = 15, height = 40)  






####################################################
#### VENN DIAGRAM — area-proportional (eulerr) ####
####################################################

# install.packages("eulerr")  # run once if not installed
library(eulerr)

# Validate counts from data_heatmap2
in_vprot <- !is.na(data_heatmap2$vprot_coef_numeric)
in_aprot <- !is.na(data_heatmap2$aprot_coef_numeric)
in_pr    <- !is.na(data_heatmap2$pr_coef_numeric)

cat("Plant only:          ", sum( in_vprot & !in_aprot & !in_pr), "\n")
cat("Animal only:         ", sum(!in_vprot &  in_aprot & !in_pr), "\n")
cat("PR only:             ", sum(!in_vprot & !in_aprot &  in_pr), "\n")
cat("Plant & Animal only: ", sum( in_vprot &  in_aprot & !in_pr), "\n")
cat("Plant & PR only:     ", sum( in_vprot & !in_aprot &  in_pr), "\n")
cat("Animal & PR only:    ", sum(!in_vprot &  in_aprot &  in_pr), "\n")
cat("All three:           ", sum( in_vprot &  in_aprot &  in_pr), "\n")
cat("Totals — Plant:", sum(in_vprot), "| Animal:", sum(in_aprot), "| PR:", sum(in_pr), "\n")

# Build area-proportional Euler diagram
fit <- euler(c(
  "Plant protein"                          = sum( in_vprot & !in_aprot & !in_pr),
  "Animal protein"                         = sum(!in_vprot &  in_aprot & !in_pr),
  "Protein ratio (P:A)"                    = sum(!in_vprot & !in_aprot &  in_pr),
  "Plant protein&Animal protein"           = sum( in_vprot &  in_aprot & !in_pr),
  "Plant protein&Protein ratio (P:A)"      = sum( in_vprot & !in_aprot &  in_pr),
  "Animal protein&Protein ratio (P:A)"     = sum(!in_vprot &  in_aprot &  in_pr),
  "Plant protein&Animal protein&Protein ratio (P:A)" = sum( in_vprot &  in_aprot &  in_pr)
))

venn_plot <- plot(fit,
     quantities = list(type = c("counts"), fontsize = 14, col = "black"),
     labels     = FALSE,
     fills      = list(fill = "transparent"),
     edges      = list(col = c("#FF1493", "#0000CD", "#1a2940"), lwd = 2))

# install.packages("ggforce")  # run once if needed
library(ggforce)

# Extract circle parameters from eulerr fit
e <- as.data.frame(fit$ellipses)
# eulerr: rows are in order of the named sets passed to euler()
# Row 1 = Plant protein, Row 2 = Animal protein, Row 3 = Protein ratio P:A
e$name  <- c("Plant protein", "Animal protein", "Protein ratio P:A")
e$total <- c(sum(in_vprot), sum(in_aprot), sum(in_pr))
e$col   <- c("#FF1493", "#0000CD", "#1a2940")

# Region count label positions — computed from exact ellipse coordinates
region_labels <- data.frame(
  x     = c(-4.0,  4.3, -0.43,  0.16, -2.5,  2.7,  0.0),
  y     = c(-2.3, -0.81,  4.1, -2.5,   0.8,  1.2, -0.2),
  label = as.character(c(
    sum( in_vprot & !in_aprot & !in_pr),   # Plant only
    sum(!in_vprot &  in_aprot & !in_pr),   # Animal only
    sum(!in_vprot & !in_aprot &  in_pr),   # PR only
    sum( in_vprot &  in_aprot & !in_pr),   # Plant & Animal only
    sum( in_vprot & !in_aprot &  in_pr),   # Plant & PR only
    sum(!in_vprot &  in_aprot &  in_pr),   # Animal & PR only
    sum( in_vprot &  in_aprot &  in_pr)    # All three
  ))
)

# Outside label positions — pushed outward from each circle center (scale = 1.2)
push_label <- function(cx, cy, r, scale = 1.2) {
  len <- sqrt(cx^2 + cy^2)
  if (len == 0) return(c(cx, cy + r * scale))
  c(cx + (cx / len) * r * scale, cy + (cy / len) * r * scale)
}

e$lx <- mapply(function(x, y, r) push_label(x, y, r)[1], e$h, e$k, e$a)
e$ly <- mapply(function(x, y, r) push_label(x, y, r)[2], e$h, e$k, e$a)

venn_gg <- ggplot() +
  geom_circle(data = e,
              aes(x0 = h, y0 = k, r = a, color = col),
              fill = "transparent", linewidth = 1.2) +
  scale_color_identity() +
  geom_text(data = region_labels,
            aes(x = x, y = y, label = label),
            size = 5) +
  geom_label(data = e,
             aes(x = lx, y = ly,
                 label = paste0(name, "\n(", total, ")")),
             fontface = "bold", size = 5,
             label.size = NA, fill = "transparent") +
  coord_equal(xlim = c(-8, 8), ylim = c(-7, 8), clip = "off") +
  theme_void()

ggsave(file.path(repo, "Heatmap/Venn_proportional.pdf"),
       plot = venn_gg, width = 10, height = 8, device = "pdf")

