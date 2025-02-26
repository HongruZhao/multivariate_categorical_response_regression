# Load necessary library
library(dplyr)

# Specify the directory containing RDS files
directory_path <- "/home/shenx/zhao1118/AD-simulation/structure_learning_2/data_overlap"

# List all RDS files in the directory
rds_files <- list.files(path = directory_path, pattern = "\\.RDS$", full.names = TRUE)

# Initialize an empty list to store data frames
df_list <- list()

# Loop over each RDS file, read it, and store the data frame in the list
for (file in rds_files) {
  df <- readRDS(file)
  df_list[[length(df_list) + 1]] <- df
}

# Initialize a data frame to store unique combinations of (p, n, Model)
unique_combinations <- unique(do.call(rbind, lapply(df_list, function(x) data.frame(p = x$p, n = x$n, Model = x$Model))))

# Initialize 'combined_data' as an empty list here before the loop
combined_data <- list()

# Iterate over unique combinations
for (i in seq_len(nrow(unique_combinations))) {
  combo <- unique_combinations[i, ]
  
  # Initialize combined_data[[i]] as a list
  combined_data[[i]] <- list()
  
  # Filter df_list for elements matching the current combination
  matching_dfs <- lapply(df_list, function(df) {
    if (df$p == combo$p && df$n == combo$n && df$Model == combo$Model) return(df)
    else return(NULL)
  })
  matching_dfs <- matching_dfs[!sapply(matching_dfs, is.null)]
  
  # Combine data for matching elements
  combined_data[[i]]$Hellinger_dist <- do.call(rbind, lapply(matching_dfs, `[[`, "Hellinger_dist"))
  combined_data[[i]]$Misclassification_rate <- do.call(rbind, lapply(matching_dfs, `[[`, "Misclassification_rate"))
  combined_data[[i]]$F_measure <- do.call(rbind, lapply(matching_dfs, `[[`, "F_measure"))
  combined_data[[i]]$G_measure <- do.call(rbind, lapply(matching_dfs, `[[`, "G_measure"))
  combined_data[[i]]$FPR_mat <- do.call(rbind, lapply(matching_dfs, `[[`, "FPR_mat"))
  combined_data[[i]]$FNR_mat <- do.call(rbind, lapply(matching_dfs, `[[`, "FNR_mat"))
  combined_data[[i]]$rep_ID <- do.call(rbind, lapply(matching_dfs, `[[`, "rep_ID"))
  
  # Store the (p, n, Model) combination for reference
  combined_data[[i]]$p <- combo$p
  combined_data[[i]]$n <- combo$n
  combined_data[[i]]$Model <- combo$Model
}


setwd("/home/shenx/zhao1118/AD-simulation/structure_learning_2")
saveRDS(combined_data, file = "overlap_summary_100rep.rds")
