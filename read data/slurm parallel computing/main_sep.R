###############################################
## Real Data Pipeline (Separate Z1, Z2, Z3)
## Using glmnet with binomial & multinomial
## Self-contained version
###############################################

# (A) Slurm or local parallel indexing
ID_index <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(ID_index)) ID_index <- 1  # fallback if not in Slurm
num_parallel <- 50  # e.g., #SBATCH --array=1-50

# If glmnet might be missing, optionally:
#if (!requireNamespace("glmnet", quietly = TRUE)) {
#  install.packages("glmnet", repos="https://cloud.r-project.org/")
#}
.libPaths(c("/panfs/jay/groups/27/shenx/zhao1118/AD-real-data/library", .libPaths()))
library(Matrix)
library(glmnet)

#########################
# (B) Model parameters
#########################
nreps <- 1000  # total replications
# must be divisible by num_parallel in typical Slurm usage

# Real data scenario
p_vec     <- c(100,200,300,400,500)   # number of columns in X (beyond the first 3 for Y)
Model_vec <- c(1)     # scheme=1
n_vec     <- c(100)   # training size(s)

params <- expand.grid("p" = p_vec, "Model" = Model_vec, "n" = n_vec)
# The dataset has total size n=202 => we typically do:
#  - train=100, validation=50, test=52

##########################################################
## (C) read the final real dataset with 202 rows
##     and 3 factor columns (Z1, Z2, Z3), plus many genes
##########################################################
dataset  <- read.csv("/home/shenx/zhao1118/AD-real-data/dataset.csv")
dfReal   <- dataset
# Convert first 3 columns to factor
dfReal[[1]] <- as.factor(dfReal[[1]])
dfReal[[2]] <- as.factor(dfReal[[2]])
dfReal[[3]] <- as.factor(dfReal[[3]])

##########################################################
## (D) Utility functions: map_cats_to_j123, one_hot_tensor
##     to build or decode Y
##########################################################
map_cats_to_j123 <- function(z1, z2, z3, lev1, lev2, lev3) {
  j1 <- match(z1, lev1)  # 1..(length(lev1))
  j2 <- match(z2, lev2)  # 1..(length(lev2))
  j3 <- match(z3, lev3)  # 1..(length(lev3))
  c(j1, j2, j3)
}

one_hot_tensor <- function(j1, j2, j3, J) {
  arr <- array(0, dim=J)  # J = c(2,2,4)
  arr[j1, j2, j3] <- 1
  as.vector(arr)
}

##########################################################
## (E) convert_response => one-hot encode (Z1,Z2,Z3)
##########################################################
convert_response <- function(dataset, J=c(2,2,4)) {
  stopifnot(ncol(dataset) >= 3)
  lev1 <- levels(dataset[[1]])
  lev2 <- levels(dataset[[2]])
  lev3 <- levels(dataset[[3]])
  nSamples <- nrow(dataset)
  
  Y <- matrix(0, nrow=prod(J), ncol=nSamples)
  for (i in seq_len(nSamples)) {
    z1 <- dataset[i,1]
    z2 <- dataset[i,2]
    z3 <- dataset[i,3]
    j123 <- map_cats_to_j123(z1, z2, z3, lev1, lev2, lev3)
    Y[, i] <- one_hot_tensor(j123[1], j123[2], j123[3], J)
  }
  Y
}

##########################################################
## (F) decode => from one-hot to (z1,z2,z3) in {1,2}x{1,2}x{1..4}
##########################################################
decode_Y_to_Z123 <- function(Y_mat, J=c(2,2,4)) {
  # Y_mat: shape = 16 x #samples if J=c(2,2,4)
  # We decode the index to z1,z2,z3 in {1,2} for z1,z2 and {1..4} for z3
  nSamples <- ncol(Y_mat)
  z1_vec <- integer(nSamples)
  z2_vec <- integer(nSamples)
  z3_vec <- integer(nSamples)
  # total states = 2*2*4 = 16
  # Suppose the row "pos" in [1..16], we want to decode that => z1,z2,z3
  # We'll do zero-based then decode
  for (col in seq_len(nSamples)) {
    pos_1 <- which(Y_mat[, col] == 1L)
    # If there's exactly one "1" per column => pos_1 is the row index in 1..16
    zero_idx <- pos_1 - 1  # now 0..15
    # decode
    # The dimension c(2,2,4) => fastest index is z1 in 0..1, next z2 in 0..1, next z3 in 0..3
    # i.e. zero_idx = z1_0 + 2*(z2_0 + 2*z3_0)
    # We'll do integer arithmetic
    z3_0 <- zero_idx %/% 4
    remainder <- zero_idx %% 4
    z2_0 <- remainder %/% 2
    z1_0 <- remainder %% 2
    # convert to 1-based
    z1_vec[col] <- z1_0 + 1
    z2_vec[col] <- z2_0 + 1
    z3_vec[col] <- z3_0 + 1
  }
  list(z1=z1_vec, z2=z2_vec, z3=z3_vec)
}

##########################################################
## (G) cross-entropy for binomial
##########################################################
compute_CE_binomial <- function(fitobj, X_val, y_val) {
  lam_seq <- fitobj$lambda
  nlam <- length(lam_seq)
  ce_vals <- numeric(nlam)
  # We interpret y_val in {1,2}, so let's map "2 => 1" and "1 => 0"
  y_num <- ifelse(y_val == 2, 1, 0)
  
  for (i in seq_len(nlam)) {
    probPos <- predict(fitobj, newx=X_val, s=lam_seq[i], type="response")[,1]
    # clamp
    probPos <- pmin(pmax(probPos, 1e-15), 1 - 1e-15)
    ce_vals[i] <- -mean( y_num*log(probPos) + (1 - y_num)*log(1 - probPos) )
  }
  ce_vals
}

##########################################################
## (H) cross-entropy for multinomial
##########################################################
compute_CE_multinomial <- function(fitobj, X_val, y_val) {
  lam_seq <- fitobj$lambda
  nlam <- length(lam_seq)
  ce_vals <- numeric(nlam)
  K <- length(unique(y_val))
  n_val <- length(y_val)
  
  for(i in seq_len(nlam)) {
    prob_array <- predict(fitobj, newx=X_val, s=lam_seq[i], type="response")
    # If 3D => prob_array has shape [n_val, K, 1]
    # If 2D => prob_array has shape [n_val, K]
    if(length(dim(prob_array)) == 3) {
      prob_mat <- prob_array[,,1]
    } else {
      prob_mat <- prob_array
    }
    # y_val in {1..K}
    chosen_prob <- numeric(n_val)
    for(j in seq_len(n_val)) {
      chosen_prob[j] <- max(prob_mat[j, y_val[j]], 1e-15)
    }
    ce_vals[i] <- -mean(log(chosen_prob))
  }
  ce_vals
}

##########################################################
## (I) make_X to create design matrix with intercept
##     (We do random column sub-sampling for demonstration)
##########################################################
make_X_old <- function(dataset, addIntercept=TRUE, colSample=10, seed=123) {
  # dataset: 1..3 are factor => columns 4.. are numeric genes
  X0full <- as.matrix(dataset[, -c(1,2,3), drop=FALSE])
  nColsFull <- ncol(X0full)
  # colSample minus 1 for intercept
  colSample <- colSample - 1
  if(nColsFull > colSample) {
    set.seed(seed)
    pick_cols <- sample(seq_len(nColsFull), colSample)
    X0full <- X0full[, pick_cols, drop=FALSE]
    cat(sprintf("   Subsampling %d of %d columns...\n", colSample, nColsFull))
  }
  X0t <- t(X0full)
  if(addIntercept) {
    X1 <- rbind(1, X0t)  # final shape => (#features) x (#samples)
    return(X1)
  } else {
    return(X0t)
  }
}

make_X <- function(dataset, addIntercept = TRUE, colSample, seed = 123) {
  # This function assumes the input 'dataset' is a data frame created by your script,
  # with the first three columns named "disease.state", "other", and "tissue" (the responses),
  # and the remaining columns being predictors.
  
  # Define the response column names.
  response_cols <- c("disease.state", "other", "tissue")
  
  # Check that the required response columns exist.
  if (!all(response_cols %in% colnames(dataset))) {
    stop("The dataset must contain the response columns: 'disease.state', 'other', and 'tissue'.")
  }
  
  # ----------------------------
  # Part 1: Separate Predictors and Responses
  # ----------------------------
  
  # Extract predictor columns (all columns except the responses).
  predictor_data <- dataset[, !(colnames(dataset) %in% response_cols), drop = FALSE]
  Xall <- as.matrix(predictor_data)
  
  # Order predictors by their Median Absolute Deviation (MAD)
  mad_values <- apply(Xall, 2, mad)
  selected_cols <- order(mad_values, decreasing = TRUE)
  
  # Retain only the top 3000 predictors (or fewer if not available)
  num_to_keep <- min(3000, ncol(Xall))
  Xall <- Xall[, selected_cols[1:num_to_keep], drop = FALSE]
  
  # ----------------------------
  # Part 2: Predictor Selection using SNP_Prune
  # ----------------------------
  
  # Define the SNP_Prune function.
  # This function iteratively selects the first predictor and removes
  # any predictors that have an absolute correlation >= cor.cutoff with it.
  SNP_Prune <- function(X, cor.cutoff = 0.7) {
    keep <- NULL
    keepnames <- NULL
    Xsorted <- X
    
    while (ncol(Xsorted) > 1) {
      # Select the first column in the current sorted matrix.
      keep <- cbind(keep, Xsorted[, 1])
      keepnames <- c(keepnames, colnames(Xsorted)[1])
      
      # If only one predictor remains, exit the loop.
      if (ncol(Xsorted) == 1) break
      
      # Compute absolute correlations between the first predictor and the remaining predictors.
      temp <- abs(cor(Xsorted[, 1], Xsorted[, -1]))
      
      # Identify predictors with correlation >= cor.cutoff.
      rm <- which(temp >= cor.cutoff)
      
      # Remove the first predictor and any predictors that are too highly correlated with it.
      Xsorted <- Xsorted[, -c(1, if (length(rm) > 0) rm + 1), drop = FALSE]
    }
    
    # If one predictor remains, add it.
    if (ncol(Xsorted) == 1) {
      keep <- cbind(keep, Xsorted)
      keepnames <- c(keepnames, colnames(Xsorted))
    }
    
    return(list("X" = keep, "snpnames" = keepnames))
  }
  
  # Determine the number of predictors to select.
  # If an intercept is to be added, we need to select one fewer predictor.
  num_to_select <- if (addIntercept) colSample - 1 else colSample
  
  # For reproducibility in SNP selection (if needed), set the seed.
  set.seed(seed)
  
  # Apply SNP_Prune with a correlation cutoff of 0.7.
  X_temp <- SNP_Prune(Xall, cor.cutoff = 0.7)
  
  # Check how many predictors remain after pruning.
  num_selected <- ncol(X_temp$X)
  if (num_selected < num_to_select) {
    warning(paste("Only", num_selected, "predictors are available after pruning; returning all available predictors."))
    num_to_select <- num_selected
  }
  
  # Select the first 'num_to_select' predictors.
  X_selected <- X_temp$X[, 1:num_to_select, drop = FALSE]
  colnames(X_selected) <- X_temp$snpnames[1:num_to_select]
  
  # If addIntercept is TRUE, add an intercept column (a column of ones) as the first predictor.
  if (addIntercept) {
    intercept <- rep(1, nrow(X_selected))
    X_selected <- cbind(Intercept = intercept, X_selected)
  }
  
  # ----------------------------
  # Part 3: Return the Predictor Matrix
  # ----------------------------
  
  final_dataset <- t(X_selected)
  return(final_dataset)
}



###################################
## (J) Main loop
###################################
for(ind_loop in seq_len(nrow(params))) {
  p     <- params[ind_loop, "p"]
  Model <- params[ind_loop, "Model"]
  n     <- params[ind_loop, "n"]
  
  for(rep_i in ID_index + (0:(nreps/num_parallel - 1))*num_parallel) {
    cat(sprintf(">>> Real data separate approach: p=%d, Model=%d, n=%d, rep_i=%d\n", 
                p, Model, n, rep_i))
    
    # 1) Convert entire data => Y (16 x 202)
    Y_full <- convert_response(dfReal, J=c(2,2,4))
    # 2) Build X
    X_full <- make_X(dfReal, addIntercept=TRUE, colSample=p, seed=2024)
    # total samples = 202
    n_train      <- n
    n_validation <- 150 - n
    n_test       <- 52
    stopifnot(n_train + n_validation + n_test == 202)
    
    # 3) Shuffle
    set.seed(2024 + rep_i)
    perm <- sample(202)
    X_full_perm <- X_full[, perm]
    Y_full_perm <- Y_full[, perm]
    
    X_train      <- X_full_perm[, 1:n_train]
    Y_train      <- Y_full_perm[, 1:n_train]
    X_validation <- X_full_perm[, (n_train+1):(n_train+n_validation)]
    Y_validation <- Y_full_perm[, (n_train+1):(n_train+n_validation)]
    X_test       <- X_full_perm[, (n_train+n_validation+1):(n_train+n_validation+n_test)]
    Y_test       <- Y_full_perm[, (n_train+n_validation+1):(n_train+n_validation+n_test)]
    
    # 4) decode => (z1,z2,z3) for train, valid, test
    trainZ <- decode_Y_to_Z123(Y_train, J=c(2,2,4))
    z1_train <- trainZ$z1  # in {1,2}
    z2_train <- trainZ$z2  # in {1,2}
    z3_train <- trainZ$z3  # in {1,2,3,4}
    
    validZ <- decode_Y_to_Z123(Y_validation, J=c(2,2,4))
    z1_valid <- validZ$z1
    z2_valid <- validZ$z2
    z3_valid <- validZ$z3
    
    testZ <- decode_Y_to_Z123(Y_test, J=c(2,2,4))
    z1_test <- testZ$z1
    z2_test <- testZ$z2
    z3_test <- testZ$z3
    
    # 5) Fit 3 separate glmnet models
    # remove intercept row => shape = n_train x (p-1)
    X_train_glm <- t(X_train[-1, ])
    z1_factor_train <- factor(z1_train, levels=c(1,2))
    z2_factor_train <- factor(z2_train, levels=c(1,2))
    z3_factor_train <- factor(z3_train, levels=c(1,2,3,4))
    
    fit_z1 <- glmnet(x=X_train_glm, y=z1_factor_train, family="binomial")
    fit_z2 <- glmnet(x=X_train_glm, y=z2_factor_train, family="binomial")
    fit_z3 <- glmnet(x=X_train_glm, y=z3_factor_train,
                     family="multinomial", type.multinomial="grouped")
    
    # 6) pick lambda by cross-entropy on validation
    X_val_glm <- t(X_validation[-1, ])
    ce_z1 <- compute_CE_binomial(fit_z1, X_val_glm, z1_valid)
    best_idx_z1 <- which.min(ce_z1)
    best_lambda_z1 <- fit_z1$lambda[best_idx_z1]
    
    ce_z2 <- compute_CE_binomial(fit_z2, X_val_glm, z2_valid)
    best_idx_z2 <- which.min(ce_z2)
    best_lambda_z2 <- fit_z2$lambda[best_idx_z2]
    
    ce_z3 <- compute_CE_multinomial(fit_z3, X_val_glm, z3_valid)
    best_idx_z3 <- which.min(ce_z3)
    best_lambda_z3 <- fit_z3$lambda[best_idx_z3]
    
    # 7) compute test misclassification (joint)
    X_test_glm <- t(X_test[-1, ])
    n_test_size <- ncol(X_test)
    
    # predict z1
    prob_z1_test <- predict(fit_z1, newx=X_test_glm, s=best_lambda_z1, type="response")[,1]
    z1hat <- ifelse(prob_z1_test >= 0.5, 2L, 1L)
    
    # predict z2
    prob_z2_test <- predict(fit_z2, newx=X_test_glm, s=best_lambda_z2, type="response")[,1]
    z2hat <- ifelse(prob_z2_test >= 0.5, 2L, 1L)
    
    # predict z3
    prob_z3_test_array <- predict(fit_z3, newx=X_test_glm, s=best_lambda_z3, type="response")
    if(length(dim(prob_z3_test_array)) == 3) {
      prob_z3_test_mat <- prob_z3_test_array[,,1]  # shape (# test obs) x # classes
    } else {
      prob_z3_test_mat <- prob_z3_test_array
    }
    z3hat <- apply(prob_z3_test_mat, 1, which.max)  # 1..4
    
    # check correctness
    correct_flags <- (z1hat == z1_test) & (z2hat == z2_test) & (z3hat == z3_test)
    num_correct <- sum(correct_flags)
    joint_misclassification_rate <- 1 - (num_correct / length(correct_flags))
    cat(sprintf("   joint misclassification = %.4f\n", joint_misclassification_rate))
    
    
    # For z1 (binomial)
    y1_true <- ifelse(z1_test == 2, 1, 0)
    CE_z1 <- -mean( y1_true * log(prob_z1_test) + (1 - y1_true) * log(1 - prob_z1_test) )
    
    # For z2 (binomial)
    y2_true <- ifelse(z2_test == 2, 1, 0)
    CE_z2 <- -mean( y2_true * log(prob_z2_test) + (1 - y2_true) * log(1 - prob_z2_test) )
    
    # For z3 (multinomial)
    # For each test sample i, take the predicted probability assigned to the true class.
    # z3_test is assumed to be in {1,2,3,4}.
    p_true_z3 <- sapply(seq_len(nrow(prob_z3_test_mat)), function(i) {
      prob_z3_test_mat[i, z3_test[i]]
    })
    CE_z3 <- -mean(log(p_true_z3))
    
    # Joint cross entropy: add the three separate cross entropies together.
    joint_CE <- CE_z1 + CE_z2 + CE_z3
    
    
    # 8) Save results => data_sep
    savename <- paste0("/home/shenx/zhao1118/AD-real-data/data_sep/Model",
                       Model, "_p", p, "_n", n, "_repID", rep_i, ".RDS")
    Results <- list(
      misclassification = joint_misclassification_rate,
      joint_CE=joint_CE,
      Model = Model,
      p = p,
      n = n,
      rep_ID = rep_i
    )
    saveRDS(Results, file=savename)
  } # end for rep_i
} # end for(ind_loop in seq_len(nrow(params)))

