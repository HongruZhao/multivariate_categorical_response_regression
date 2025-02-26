# Load the ggplot2 package for plotting
library(ggplot2)

# Read in the cross entropy summary and the overlap summary data from RDS files
CE_summary <- readRDS("Cross_entropy_summary_100rep.rds")
overlap_summary <- readRDS("overlap_summary_100rep.rds")

# Extract unique settings (p, n, Model) from the CE_summary list
unique_setting <- unique(do.call(rbind, lapply(CE_summary, function(x) data.frame(p = x$p, n = x$n, Model = x$Model))))

# Extract unique settings without the 'n' variable (p and Model only)
unique_setting_without_n <- unique(do.call(rbind, lapply(CE_summary, function(x) data.frame(p = x$p,   Model = x$Model))))

# Function to convert two summary lists into a combined data frame for a given metric
list_to_df_2 <- function(summary_list, summary_list_2, metric) {
  # Initialize a list to store processed elements for each repetition
  processed_list <- list()
  
  # Loop through each summary element in the summary list
  for (i in seq_along(summary_list)) {
    summary_element <- summary_list[[i]]
    summary_element_2 <- summary_list_2[[i]]
    
    # Create a new list to store processed data for the current repetition
    list_new <- list()
    # Combine the metric data from both summary elements column-wise
    list_new[[metric]] <- cbind(summary_element[[metric]], summary_element_2[[metric]])
    
    # Define a full vector of estimator names
    names_vector_full <- c("O-Mult", "O-Pois", "L-Mult", "L-Pois", "G-Mult", "G-Pois", "G-Mult-theta", "Sep-Mult", "Oracle")
    # Use only as many names as there are columns in the combined metric data
    names_vector <- names_vector_full[1:ncol(list_new[[metric]])]
    num_rows <- nrow(list_new[[metric]])
    
    # Create a matrix of estimator names to assign to each row of the metric data
    list_new$nam <- matrix(names_vector, nrow = num_rows, ncol = length(names_vector), byrow = TRUE)
    
    # Get the dimensions of the metric matrix
    dim_metric <- dim(list_new[[metric]])
    # Create matrices for n, Model, and p values, matching the dimensions of the metric data
    list_new$n_mat <- matrix(rep(summary_element$n, prod(dim_metric)), nrow = dim_metric[1], ncol = dim_metric[2])
    list_new$Model_mat <- matrix(rep(summary_element$Model, prod(dim_metric)), nrow = dim_metric[1], ncol = dim_metric[2])
    list_new$p_mat <- matrix(rep(summary_element$p, prod(dim_metric)), nrow = dim_metric[1], ncol = dim_metric[2])
    
    # Store the processed list for the current repetition
    processed_list[[i]] <- list_new
  }
  
  # Combine the processed elements from all repetitions into one result list
  combined_results <- list()
  elements_to_combine <- names(processed_list[[1]])
  
  # For each element (metric, names, n, Model, p), combine all repetitions row-wise
  for (element_name in elements_to_combine) {
    element_list <- lapply(processed_list, function(x) x[[element_name]])
    combined_matrix <- do.call(rbind, element_list)
    combined_results[[element_name]] <- combined_matrix
  }
  
  # Convert the combined matrices to vectors and store them in a final data frame
  final_df <- data.frame(
    Hellinger_dist = as.vector(combined_results[[metric]]),
    nam = as.vector(combined_results$nam),
    n_mat = as.vector(combined_results$n_mat),
    Model_mat = as.vector(combined_results$Model_mat),
    p_mat = as.vector(combined_results$p_mat)
  )
  
  return(final_df)
}

# Define color mapping for the different estimators using hex color codes
color_mapping <- c(
  "O-Mult"         = "#ADD8E6",  # Light Blue
  "L-Mult"         = "#0000CD",  # Medium Blue
  "G-Mult"         = "#00008B",  # Dark Blue
  "O-Pois"         = "#FFDAB9",  # Light Orange (Peach Puff)
  "L-Pois"         = "#FFA500",  # Orange
  "G-Pois"         = "#FF8C00",  # Dark Orange
  "G-Mult-theta"   = "#800080",  # Purple
  "Sep-Mult"       = "#FFC0CB",  # Pink
  "Oracle"         = "#A9A9A9"   # Dark Grey
)

# Print the names of the elements in the first CE_summary list element
names(CE_summary[[1]])
# Expected output (example):
# "Hellinger_dist"         "Misclassification_rate" "F_measure"              
# "G_measure"              "FPR_mat"                "FNR_mat"                
# "rep_ID"                 "p"                      "n"                      
# "Model"                  

# Set the metric name for the y-axis; you can switch this to "Misclassification_rate" if needed
name_of_ylab = "Hellinger_dist"
#### To plot misclassification rate, uncomment the next line:
#### name_of_ylab = "Misclassification_rate"

# Generate the final data frame from the overlap and cross entropy summaries using the specified metric
final_df <- list_to_df_2(overlap_summary, CE_summary, name_of_ylab)

# Display the unique settings without n (p and Model)
unique_setting_without_n

# Loop over the first 6 unique settings to create and save boxplots for each setting
for (ijjkkk in 1:6) {
  # Create an index for rows in final_df that match the current setting (p and Model)
  iinndd = final_df$p_mat == unique_setting_without_n[ijjkkk, 1] & final_df$Model_mat == unique_setting_without_n[ijjkkk, 2]
  dat = final_df[iinndd, ]
  # Convert n_mat to character type for proper factor ordering in plots
  dat$n_mat = as.character(dat$n_mat)
  # Set the factor levels for n_mat to ensure they appear in the desired order
  dat$n_mat <- factor(dat$n_mat, levels = c("100", "300", "500", "1000", "2000"))
  
  # Set the factor levels for the estimator names based on the chosen metric
  if (name_of_ylab == "Misclassification_rate") {
    dat$nam = factor(dat$nam, levels = c("O-Mult", "O-Pois", "L-Mult", "L-Pois", "G-Mult", "G-Pois", "G-Mult-theta", "Sep-Mult", "Oracle"))
  } else {
    dat$nam = factor(dat$nam, levels = c("O-Mult", "O-Pois", "L-Mult", "L-Pois", "G-Mult", "G-Pois", "G-Mult-theta", "Sep-Mult"))
  }
  
  # Set the y-axis label based on the metric being plotted
  if (name_of_ylab == "Misclassification_rate") {
    ylab_lab = "Misclassification Rate"
  } else {
    ylab_lab = "Hellinger Distance"
  }
  
  # Construct the plot title dynamically using the current p and Model (Scheme) values
  plot_title <- paste("p =", unique_setting_without_n[ijjkkk, 1], "and Scheme =", unique_setting_without_n[ijjkkk, 2])
  
  # Create a boxplot with ggplot2
  p1 <- ggplot(dat, aes(y = Hellinger_dist, x = n_mat, fill = nam)) +
    geom_boxplot(outlier.shape = NA, outlier.size = 0, lwd = 0.05) +
    theme_bw() +
    xlab("n") +  # Set the x-axis label
    ylab(ylab_lab) +   # Set the y-axis label based on the metric
    labs(fill = "Estimators") +  # Set the legend title for the estimator names
    ggtitle(plot_title) +  # Set the plot title dynamically
    scale_fill_manual(values = color_mapping) +  # Use custom colors for the fill based on estimator
    theme(
      plot.title = element_text(hjust = 0.5),  # Center the plot title
      # legend.position = "none"  # Option to remove the legend
      legend.position = "bottom"  # Place the legend below the plot
      # legend.position = "none"  # Alternative: remove the legend
    )
  
  # Define the filename for saving the plot (as a PDF)
  filenames <- paste("plot/", ylab_lab, "_p", unique_setting_without_n[ijjkkk, 1], "_scheme", unique_setting_without_n[ijjkkk, 2], ".pdf", sep = "")
  
  # Save the generated plot to the specified file with given width and height
  ggsave(filenames, plot = p1, width = 6, height = 4)
}
