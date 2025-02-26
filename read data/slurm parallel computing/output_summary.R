# output_summary.R
# This script reads all the RDS files produced from the Multinomial Log‐Linear pipeline
# (stored in the "data_CE" directory) and summarizes both the misclassification rates
# and the joint cross-entropies (Joint_CE) by (p, n, Model).

library(dplyr)

# Directory of RDS files
directory_path <- "/home/shenx/zhao1118/AD-real-data/data_CE"

# List all RDS files in the directory
rds_files <- list.files(path = directory_path, pattern = "\\.RDS$", full.names = TRUE)

# Initialize a list to store each file's results
all_results_list <- list()

for (file in rds_files) {
  # Read the RDS object (each is a list with Misclassification_rate, Joint_CE, Model, p, n, rep_ID)
  result_obj <- readRDS(file)
  
  # 'Misclassification_rate' is a 1 x 13 matrix; convert it to a numeric vector
  mis_vec <- as.numeric(result_obj$Misclassification_rate)
  # 'Joint_CE' is also a 1 x 13 matrix; convert it to a numeric vector
  joint_vec <- as.numeric(result_obj$Joint_CE)
  
  # Create a data frame row with labeled columns for misclassification and joint CE
  row_df <- data.frame(
    p       = result_obj$p,
    n       = result_obj$n,
    Model   = result_obj$Model,
    rate=result_obj$rate,
    rep_ID  = result_obj$rep_ID,
    d3_LMult         = mis_vec[1],
    d3_LPois         = mis_vec[2],
    d3_GMult         = mis_vec[3],
    d3_GPois         = mis_vec[4],
    d2_LMult         = mis_vec[5],
    d2_LPois         = mis_vec[6],
    d2_GMult         = mis_vec[7],
    d2_GPois         = mis_vec[8],
    d1_LMult         = mis_vec[9],
    d1_LPois         = mis_vec[10],
    d1_GMult         = mis_vec[11],
    d1_GPois         = mis_vec[12],
    Theta_GroupLasso = mis_vec[13],
    d3_LMult_joint_CE         = joint_vec[1],
    d3_LPois_joint_CE         = joint_vec[2],
    d3_GMult_joint_CE         = joint_vec[3],
    d3_GPois_joint_CE         = joint_vec[4],
    d2_LMult_joint_CE         = joint_vec[5],
    d2_LPois_joint_CE         = joint_vec[6],
    d2_GMult_joint_CE         = joint_vec[7],
    d2_GPois_joint_CE         = joint_vec[8],
    d1_LMult_joint_CE         = joint_vec[9],
    d1_LPois_joint_CE         = joint_vec[10],
    d1_GMult_joint_CE         = joint_vec[11],
    d1_GPois_joint_CE         = joint_vec[12],
    Theta_GroupLasso_joint_CE = joint_vec[13]
  )
  
  all_results_list[[length(all_results_list) + 1]] <- row_df
}

# Combine all individual rows into one master data frame
all_results_df <- bind_rows(all_results_list)



# Step 1. For each combination of p, n, Model, and rep_ID,
# compute the minimum value of each target variable over the different rate values.
min_results <- all_results_df %>%
  group_by(p, n, Model, rep_ID) %>%
  summarize(
    d3_LMult = min(d3_LMult, na.rm = TRUE),
    d3_LPois = min(d3_LPois, na.rm = TRUE),
    d3_GMult = min(d3_GMult, na.rm = TRUE),
    d3_GPois = min(d3_GPois, na.rm = TRUE),
    d2_LMult = min(d2_LMult, na.rm = TRUE),
    d2_LPois = min(d2_LPois, na.rm = TRUE),
    d2_GMult = min(d2_GMult, na.rm = TRUE),
    d2_GPois = min(d2_GPois, na.rm = TRUE),
    d1_LMult = min(d1_LMult, na.rm = TRUE),
    d1_LPois = min(d1_LPois, na.rm = TRUE),
    d1_GMult = min(d1_GMult, na.rm = TRUE),
    d1_GPois = min(d1_GPois, na.rm = TRUE),
    Theta_GroupLasso = min(Theta_GroupLasso, na.rm = TRUE),
    d3_LMult_joint_CE = min(d3_LMult_joint_CE, na.rm = TRUE),
    d3_LPois_joint_CE = min(d3_LPois_joint_CE, na.rm = TRUE),
    d3_GMult_joint_CE = min(d3_GMult_joint_CE, na.rm = TRUE),
    d3_GPois_joint_CE = min(d3_GPois_joint_CE, na.rm = TRUE),
    d2_LMult_joint_CE = min(d2_LMult_joint_CE, na.rm = TRUE),
    d2_LPois_joint_CE = min(d2_LPois_joint_CE, na.rm = TRUE),
    d2_GMult_joint_CE = min(d2_GMult_joint_CE, na.rm = TRUE),
    d2_GPois_joint_CE = min(d2_GPois_joint_CE, na.rm = TRUE),
    d1_LMult_joint_CE = min(d1_LMult_joint_CE, na.rm = TRUE),
    d1_LPois_joint_CE = min(d1_LPois_joint_CE, na.rm = TRUE),
    d1_GMult_joint_CE = min(d1_GMult_joint_CE, na.rm = TRUE),
    d1_GPois_joint_CE = min(d1_GPois_joint_CE, na.rm = TRUE),
    Theta_GroupLasso_joint_CE = min(Theta_GroupLasso_joint_CE, na.rm = TRUE)
  ) %>%
  ungroup()

# Step 2. For each combination of p, n, and Model, compute the sample mean and SD across rep_ID.
final_summary <- min_results %>%
  group_by(p, n, Model) %>%
  summarize(
    d3_LMult_mean               = mean(d3_LMult, na.rm = TRUE),
    d3_LMult_sd                 = sd(d3_LMult, na.rm = TRUE),
    d3_LPois_mean               = mean(d3_LPois, na.rm = TRUE),
    d3_LPois_sd                 = sd(d3_LPois, na.rm = TRUE),
    d3_GMult_mean               = mean(d3_GMult, na.rm = TRUE),
    d3_GMult_sd                 = sd(d3_GMult, na.rm = TRUE),
    d3_GPois_mean               = mean(d3_GPois, na.rm = TRUE),
    d3_GPois_sd                 = sd(d3_GPois, na.rm = TRUE),
    
    d2_LMult_mean               = mean(d2_LMult, na.rm = TRUE),
    d2_LMult_sd                 = sd(d2_LMult, na.rm = TRUE),
    d2_LPois_mean               = mean(d2_LPois, na.rm = TRUE),
    d2_LPois_sd                 = sd(d2_LPois, na.rm = TRUE),
    d2_GMult_mean               = mean(d2_GMult, na.rm = TRUE),
    d2_GMult_sd                 = sd(d2_GMult, na.rm = TRUE),
    d2_GPois_mean               = mean(d2_GPois, na.rm = TRUE),
    d2_GPois_sd                 = sd(d2_GPois, na.rm = TRUE),
    
    d1_LMult_mean               = mean(d1_LMult, na.rm = TRUE),
    d1_LMult_sd                 = sd(d1_LMult, na.rm = TRUE),
    d1_LPois_mean               = mean(d1_LPois, na.rm = TRUE),
    d1_LPois_sd                 = sd(d1_LPois, na.rm = TRUE),
    d1_GMult_mean               = mean(d1_GMult, na.rm = TRUE),
    d1_GMult_sd                 = sd(d1_GMult, na.rm = TRUE),
    d1_GPois_mean               = mean(d1_GPois, na.rm = TRUE),
    d1_GPois_sd                 = sd(d1_GPois, na.rm = TRUE),
    
    Theta_GroupLasso_mean       = mean(Theta_GroupLasso, na.rm = TRUE),
    Theta_GroupLasso_sd         = sd(Theta_GroupLasso, na.rm = TRUE),
    
    d3_LMult_joint_CE_mean      = mean(d3_LMult_joint_CE, na.rm = TRUE),
    d3_LMult_joint_CE_sd        = sd(d3_LMult_joint_CE, na.rm = TRUE),
    d3_LPois_joint_CE_mean      = mean(d3_LPois_joint_CE, na.rm = TRUE),
    d3_LPois_joint_CE_sd        = sd(d3_LPois_joint_CE, na.rm = TRUE),
    d3_GMult_joint_CE_mean      = mean(d3_GMult_joint_CE, na.rm = TRUE),
    d3_GMult_joint_CE_sd        = sd(d3_GMult_joint_CE, na.rm = TRUE),
    d3_GPois_joint_CE_mean      = mean(d3_GPois_joint_CE, na.rm = TRUE),
    d3_GPois_joint_CE_sd        = sd(d3_GPois_joint_CE, na.rm = TRUE),
    
    d2_LMult_joint_CE_mean      = mean(d2_LMult_joint_CE, na.rm = TRUE),
    d2_LMult_joint_CE_sd        = sd(d2_LMult_joint_CE, na.rm = TRUE),
    d2_LPois_joint_CE_mean      = mean(d2_LPois_joint_CE, na.rm = TRUE),
    d2_LPois_joint_CE_sd        = sd(d2_LPois_joint_CE, na.rm = TRUE),
    d2_GMult_joint_CE_mean      = mean(d2_GMult_joint_CE, na.rm = TRUE),
    d2_GMult_joint_CE_sd        = sd(d2_GMult_joint_CE, na.rm = TRUE),
    d2_GPois_joint_CE_mean      = mean(d2_GPois_joint_CE, na.rm = TRUE),
    d2_GPois_joint_CE_sd        = sd(d2_GPois_joint_CE, na.rm = TRUE),
    
    d1_LMult_joint_CE_mean      = mean(d1_LMult_joint_CE, na.rm = TRUE),
    d1_LMult_joint_CE_sd        = sd(d1_LMult_joint_CE, na.rm = TRUE),
    d1_LPois_joint_CE_mean      = mean(d1_LPois_joint_CE, na.rm = TRUE),
    d1_LPois_joint_CE_sd        = sd(d1_LPois_joint_CE, na.rm = TRUE),
    d1_GMult_joint_CE_mean      = mean(d1_GMult_joint_CE, na.rm = TRUE),
    d1_GMult_joint_CE_sd        = sd(d1_GMult_joint_CE, na.rm = TRUE),
    d1_GPois_joint_CE_mean      = mean(d1_GPois_joint_CE, na.rm = TRUE),
    d1_GPois_joint_CE_sd        = sd(d1_GPois_joint_CE, na.rm = TRUE),
    
    Theta_GroupLasso_joint_CE_mean = mean(Theta_GroupLasso_joint_CE, na.rm = TRUE),
    Theta_GroupLasso_joint_CE_sd   = sd(Theta_GroupLasso_joint_CE, na.rm = TRUE),
    
    .groups = "drop"
  )

# Save the summarized results to an RDS file
saveRDS(final_summary, file = "/home/shenx/zhao1118/AD-real-data/Cross_entropy_summary_1000rep.rds")

# Optional: Write to CSV if desired
# write.csv(final_summary, file = "/home/shenx/zhao1118/AD-real-data/Cross_entropy_summary_1000rep.csv", row.names = FALSE)

