library(dplyr)

# Helper function to read a summary file and drop the first three columns
read_and_clean <- function(file) {
  data <- readRDS(file)
  data[ , -c(1:3)]
}

# Read the summary files
overlap_summary_new <- read_and_clean(file.path(directory_path, "overlap_summary_1000rep.rds"))
penalty_summary_new <- read_and_clean(file.path(directory_path, "Cross_entropy_summary_1000rep.rds"))
sep_summary_new     <- read_and_clean(file.path(directory_path, "sep_summary_1000rep.rds"))

# Function to extract misclassification and CE parts given the number of methods and the offset for CE columns
extract_summary <- function(data, n_methods, ce_offset) {
  list(
    mis_mean = data[, (1:n_methods) * 2 - 1],
    mis_sd   = data[, (1:n_methods) * 2],
    ce_mean  = data[, ((1:n_methods) * 2 - 1) + ce_offset],
    ce_sd    = data[, ((1:n_methods) * 2) + ce_offset]
  )
}

# Extract the matrices:
# For overlap (6 methods; mis stats are in columns 1–12, and CE stats start at column 13)
overlap <- extract_summary(overlap_summary_new, 6, ce_offset = 12)

# For penalty (13 methods; mis stats are in columns 1–26, and CE stats start at column 27)
penalty <- extract_summary(penalty_summary_new, 13, ce_offset = 26)

# For sep (assumed to be 1 method; columns are fixed)
sep <- list(
  mis_mean = sep_summary_new[, 1],
  mis_sd   = sep_summary_new[, 2],
  ce_mean  = sep_summary_new[, 3],
  ce_sd    = sep_summary_new[, 4]
)

# Combine the matrices column-wise
combined_mis_mean <- cbind(penalty$mis_mean, overlap$mis_mean, sep$mis_mean)
combined_mis_sd   <- cbind(penalty$mis_sd,   overlap$mis_sd,   sep$mis_sd)
combined_CE_mean  <- cbind(penalty$ce_mean,  overlap$ce_mean,  sep$ce_mean)
combined_CE_sd    <- cbind(penalty$ce_sd,    overlap$ce_sd,    sep$ce_sd)

# Define method names (adjust as needed)
name_vec <- c("d3_LMult", "d3_LPois", "d3_GMult", "d3_GPois",
              "d2_LMult", "d2_LPois", "d2_GMult", "d2_GPois",
              "d1_LMult", "d1_LPois", "d1_GMult", "d1_GPois",
              "Theta_GroupLasso", "d3_OMult", "d3_OPois",
              "d2_OMult", "d2_OPois", "d1_OMult", "d1_OPois",
              "sep")

colnames(combined_mis_mean) <- name_vec
rownames(combined_mis_mean) <- c("p=100", "p=200", "p=300", "p=400", "p=500")

colnames(combined_mis_sd)   <- name_vec
rownames(combined_mis_sd) <- c("p=100", "p=200", "p=300", "p=400", "p=500")

colnames(combined_CE_mean)  <- name_vec
rownames(combined_CE_mean) <- c("p=100", "p=200", "p=300", "p=400", "p=500")

colnames(combined_CE_sd)    <- name_vec
rownames(combined_CE_sd) <- c("p=100", "p=200", "p=300", "p=400", "p=500")



library(xtable)

# 1) Transpose the matrix
mat_t <- t(combined_mis_mean)
round(mat_t, 4)
# 2) Create the xtable object
mis_tab <- xtable(
  mat_t,
  caption = "Misclassification rates",
  label = "tab:mis",
  digits=4
)


# 3) Print with desired alignment and digits
print(
  mis_tab,
  align = c("l", rep("r", ncol(mat_t))),
  include.rownames = TRUE,
  sanitize.text.function = identity,
  type = "latex",
  digits = 4  # Add digits=4 here to control decimal places
)



# 1) Transpose the matrix
mat_t <- t(combined_CE_mean)

# 2) Create the xtable object
mis_tab <- xtable(
  mat_t,
  caption = "Empirical cross entropy",
  label = "tab:ce",
  digits = 4
)


# 3) Print with desired alignment and digits
print(
  mis_tab,
  align = c("l", rep("r", ncol(mat_t))),
  include.rownames = TRUE,
  sanitize.text.function = identity,
  type = "latex",
  digits = 4
)

