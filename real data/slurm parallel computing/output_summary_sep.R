# output_summary_sep.R
# This script reads in all the RDS result files produced by the separate real data pipeline
# (stored in the "data_sep" directory) and summarizes the misclassification rates and joint
# cross-entropies by (p, n, Model). The final summary is saved as an RDS file.

library(dplyr)

# Directory containing the separate approach results
directory_path <- "/home/shenx/zhao1118/AD-real-data/data_sep"

# List all RDS files in the directory
rds_files <- list.files(path = directory_path, pattern = "\\.RDS$", full.names = TRUE)

# Initialize a list to store each file's results
all_results_list <- list()

# Loop over each RDS file and extract the relevant information
for (file in rds_files) {
  # Read the RDS object (each is a list with misclassification, joint_CE, Model, p, n, rep_ID)
  result_obj <- readRDS(file)
  
  # Create a data frame row from the result
  row_df <- data.frame(
    p = result_obj$p,
    n = result_obj$n,
    Model = result_obj$Model,
    rep_ID = result_obj$rep_ID,
    misclassification = result_obj$misclassification,
    joint_CE = result_obj$joint_CE
  )
  
  # Append this row to our list
  all_results_list[[length(all_results_list) + 1]] <- row_df
}

# Combine all individual rows into one master data frame
all_results_df <- bind_rows(all_results_list)

# Summarize the results by (p, n, Model): compute mean and standard deviation for misclassification and joint_CE
final_summary <- all_results_df %>%
  group_by(p, n, Model) %>%
  summarize(
    misclassification_mean = mean(misclassification, na.rm = TRUE),
    misclassification_sd   = sd(misclassification, na.rm = TRUE),
    joint_CE_mean          = mean(joint_CE, na.rm = TRUE),
    joint_CE_sd            = sd(joint_CE, na.rm = TRUE),
    .groups = "drop"
  )

# Save the summarized results to an RDS file
saveRDS(final_summary, file = "/home/shenx/zhao1118/AD-real-data/sep_summary_1000rep.rds")

# Optional: Write to CSV if desired
# write.csv(final_summary, file = "/home/shenx/zhao1118/AD-real-data/sep_summary_1000rep.csv", row.names = FALSE)
