# output_summary_overlap.R
# This script reads all the RDS files produced from the Overlapping Group Lasso pipeline
# (stored in the "data_overlap" directory) and summarizes both the misclassification rates
# and the joint cross-entropies (Joint_CE) by (p, n, Model).

library(dplyr)

# Specify directory containing overlapping approach results
directory_path <- "/home/shenx/zhao1118/AD-real-data/data_overlap"

# List all RDS files in the directory
rds_files <- list.files(path = directory_path, pattern = "\\.RDS$", full.names = TRUE)

all_rows <- list()  # Will store a row of data for each RDS file

for (file in rds_files) {
  # Load the RDS object
  result_obj <- readRDS(file)
  
  # mis_class is a 1 x 6 matrix; convert it to a numeric vector
  mis_vec <- as.numeric(result_obj$Misclassification_rate)
  # Joint_CE is a 1 x 6 matrix; convert it to a numeric vector
  joint_vec <- as.numeric(result_obj$Joint_CE)
  
  # Create a data frame with the (p, n, Model, rep_ID) and both sets of metrics
  row_df <- data.frame(
    p        = result_obj$p,
    n        = result_obj$n,
    Model    = result_obj$Model,
    rate     = result_obj$rate,
    rep_ID   = result_obj$rep_ID,
    d3_OMult = mis_vec[1],
    d3_OPois = mis_vec[2],
    d2_OMult = mis_vec[3],
    d2_OPois = mis_vec[4],
    d1_OMult = mis_vec[5],
    d1_OPois = mis_vec[6],
    d3_OMult_joint_CE = joint_vec[1],
    d3_OPois_joint_CE = joint_vec[2],
    d2_OMult_joint_CE = joint_vec[3],
    d2_OPois_joint_CE = joint_vec[4],
    d1_OMult_joint_CE = joint_vec[5],
    d1_OPois_joint_CE = joint_vec[6]
  )
  
  all_rows[[length(all_rows) + 1]] <- row_df
}

# Combine all rows into one data frame
all_results_df <- bind_rows(all_rows)


# Step 1. For each replicate (p, n, Model, rep_ID), compute the minimum value of each metric.
min_results <- all_results_df %>%
  group_by(p, n, Model, rep_ID) %>%
  summarize(
    d3_OMult            = min(d3_OMult, na.rm = TRUE),
    d3_OPois            = min(d3_OPois, na.rm = TRUE),
    d2_OMult            = min(d2_OMult, na.rm = TRUE),
    d2_OPois            = min(d2_OPois, na.rm = TRUE),
    d1_OMult            = min(d1_OMult, na.rm = TRUE),
    d1_OPois            = min(d1_OPois, na.rm = TRUE),
    d3_OMult_joint_CE   = min(d3_OMult_joint_CE, na.rm = TRUE),
    d3_OPois_joint_CE   = min(d3_OPois_joint_CE, na.rm = TRUE),
    d2_OMult_joint_CE   = min(d2_OMult_joint_CE, na.rm = TRUE),
    d2_OPois_joint_CE   = min(d2_OPois_joint_CE, na.rm = TRUE),
    d1_OMult_joint_CE   = min(d1_OMult_joint_CE, na.rm = TRUE),
    d1_OPois_joint_CE   = min(d1_OPois_joint_CE, na.rm = TRUE)
  ) %>%
  ungroup()

# Step 2. Now summarize by (p, n, Model): compute sample mean and sample standard deviation.
final_summary <- min_results %>%
  group_by(p, n, Model) %>%
  summarize(
    d3_OMult_mean            = mean(d3_OMult, na.rm = TRUE),
    d3_OMult_sd              = sd(d3_OMult, na.rm = TRUE),
    d3_OPois_mean            = mean(d3_OPois, na.rm = TRUE),
    d3_OPois_sd              = sd(d3_OPois, na.rm = TRUE),
    d2_OMult_mean            = mean(d2_OMult, na.rm = TRUE),
    d2_OMult_sd              = sd(d2_OMult, na.rm = TRUE),
    d2_OPois_mean            = mean(d2_OPois, na.rm = TRUE),
    d2_OPois_sd              = sd(d2_OPois, na.rm = TRUE),
    d1_OMult_mean            = mean(d1_OMult, na.rm = TRUE),
    d1_OMult_sd              = sd(d1_OMult, na.rm = TRUE),
    d1_OPois_mean            = mean(d1_OPois, na.rm = TRUE),
    d1_OPois_sd              = sd(d1_OPois, na.rm = TRUE),
    
    d3_OMult_joint_CE_mean   = mean(d3_OMult_joint_CE, na.rm = TRUE),
    d3_OMult_joint_CE_sd     = sd(d3_OMult_joint_CE, na.rm = TRUE),
    d3_OPois_joint_CE_mean   = mean(d3_OPois_joint_CE, na.rm = TRUE),
    d3_OPois_joint_CE_sd     = sd(d3_OPois_joint_CE, na.rm = TRUE),
    d2_OMult_joint_CE_mean   = mean(d2_OMult_joint_CE, na.rm = TRUE),
    d2_OMult_joint_CE_sd     = sd(d2_OMult_joint_CE, na.rm = TRUE),
    d2_OPois_joint_CE_mean   = mean(d2_OPois_joint_CE, na.rm = TRUE),
    d2_OPois_joint_CE_sd     = sd(d2_OPois_joint_CE, na.rm = TRUE),
    d1_OMult_joint_CE_mean   = mean(d1_OMult_joint_CE, na.rm = TRUE),
    d1_OMult_joint_CE_sd     = sd(d1_OMult_joint_CE, na.rm = TRUE),
    d1_OPois_joint_CE_mean   = mean(d1_OPois_joint_CE, na.rm = TRUE),
    d1_OPois_joint_CE_sd     = sd(d1_OPois_joint_CE, na.rm = TRUE),
    .groups = "drop"
  )
# Save the summarized results to an RDS file
saveRDS(final_summary, file = "/home/shenx/zhao1118/AD-real-data/overlap_summary_1000rep.rds")

# Optional: Write to CSV if desired
# write.csv(final_summary, file = "/home/shenx/zhao1118/AD-real-data/overlap_summary_100rep.csv", row.names = FALSE)

