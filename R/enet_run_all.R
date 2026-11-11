#### Run full elastic-net signature pipeline (plant, animal, ratio) for discovery cohort.
#### setwd to R/; edit 00_paths.R first.

source("00_paths.R")
source("enet_signature_config.R")
source("enet_step1_train_test_100reps.R")
source("enet_step2_leave_one_out_training.R")
source("enet_step3_merge_final_scores.R")

cat("=== Step 1: train/test, 100 repetitions ===\n")
enet_run_all_step1()

cat("=== Step 2: leave-one-out scores (training set) ===\n")
enet_run_all_step2()

cat("=== Step 3: merge scores to full cohort ===\n")
enet_run_all_step3()

cat("Done. Merge files are under repo/Enet/<Plant|Animal|Protein_ratio>/.\n")
