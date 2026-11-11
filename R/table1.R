
########################################
### Table 1 and Supplemental Table 1 ###
########################################


#load the imputed and ztransformed metabolites
load("/udd/nhtes/git/protein-metab-mortality/imp_ztran_met_PR_multi_ep.RData") 
gdata::keep(AllCohorts_clean, sure=T)

AllCohorts_clean_ztransformed <- AllCohorts_clean

## convert categorical variables to factors
factor=c("white", "fast", "multivit", "htn", "hchol", "diab", "hbpmed", "cholmed", "mifh", "dbfh",
         "cacoall", "labcode", "cohort", "sex") 
for (i in 1:length(factor)){
  AllCohorts_clean_ztransformed[, factor[i]]=as.factor(AllCohorts_clean_ztransformed[, factor[i]])
}


library(tableone)
table1_cohort <- CreateTableOne(data = AllCohorts_clean_ztransformed,
                                vars = c("ageyr","sex","white", "fast", "multivit", "diab", "htn", "hchol",
                                         "hbpmed", "cholmed", "calor_avg", "smoke_avg", "bmi_avg", "th_avg",
                                         "alco_avg", "AHEI_mod_avg", "totprot_avg", "totprot_avg_eadj", "totprot_avg_p", "aprot_avg", "aprot_avg_p", "aprot_avg_eadj", "vprot_avg", "vprot_avg_p", "vprot_avg_eadj", "protein_ratio"),
                                factorVars = c("sex","white", "fast", "multivit", "diab", "htn", "hchol", "hbpmed", "cholmed", "smoke_avg"),
                                strata = "cohort", includeNA = T, addOverall = T)
csv_cohort <- print(table1_cohort,
                    printToggle = FALSE)

# table1_healthylifestyle <- CreateTableOne(data = AllCohorts_clean_ztransformed,
#                                 vars = c("ageyr","sex","white", "fast", "multivit", "diab", "htn", "hchol",
#                                          "hbpmed", "cholmed", "calor_avg", "smoke_avg", "bmi_avg", "th_avg",
#                                          "alco_avg", "AHEI_mod_avg"),
#                                 factorVars = c("sex","white", "fast", "multivit", "diab", "htn", "hchol", "hbpmed", "cholmed", "healthy_lifestyle_score_AHEI", "smoke_avg"),
#                                 strata = "healthy_lifestyle_score_AHEI", includeNA = T, addOverall = T)
# csv_healthylifestyle<- print(table1_healthylifestyle,
#                     printToggle = FALSE)

library(openxlsx)
write.xlsx(list("bycohort"=csv_cohort),
           rowNames = T,
           file = "/udd/nhtes/git/protein-metab-mortality/Tables/table1.xlsx")


###################
### Table 1 WHI ###
###################


#load the imputed and ztransformed metabolites
WHI <- read.csv("/udd/nhtes/git/protein-metab-mortality/WHI-replication/WHI_replic.csv")

## convert categorical variables to factors
factor=c("ETHNIC", "SMOKING") 
for (i in 1:length(factor)){
  WHI[, factor[i]]=as.factor(WHI[, factor[i]])
}


library(tableone)
table1_cohort <- CreateTableOne(data = WHI,
                                vars = c("AGE","ETHNIC", "BMI", "SMOKING", "F60ENRGY", "MSMINWK",
                                          "ahei_mod", "totprot_avg_p", "vprot_avg_p", "aprot_avg_p",  "protein_ratio"),
                                factorVars = c("ETHNIC", "SMOKING"))



# n                                       2092         
#   AGE (mean (SD))                        67.08 (6.83)  
#   ETHNIC (%)                                           
#      American Indian or Alaskan Native      13 ( 0.6)  
#      Asian or Pacific Islander              33 ( 1.6)  
#      Black or African-American             242 (11.6)  
#      Hispanic/Latino                        58 ( 2.8)  
#      Other                                  49 ( 2.3)  
#      White (not of Hispanic origin)       1697 (81.1)  
#   BMI (mean (SD))                        28.60 (5.98)  
#   SMOKING (%)                                          
#                                             24 ( 1.1)  
#      Current Smoker                        240 (11.5)  
#      Never Smoked                          993 (47.5)  
#      Past Smoker                           835 (39.9)  
#   F60ENRGY (mean (SD))                 1574.08 (579.74)
#   MSMINWK (mean (SD))                    84.08 (130.98)
#   ahei_mod (mean (SD))                   38.74 (8.67)  
#   totprot_avg_p (mean (SD))              16.72 (3.12)  
#   vprot_avg_p (mean (SD))                 4.99 (1.19)  
#   aprot_avg_p (mean (SD))                11.71 (3.11)  
#   protein_ratio (mean (SD))               0.47 (0.20) 