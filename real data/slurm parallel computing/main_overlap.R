############################################################
# Overlapping Group Lasso on Real Data (6 Cases) 
#
# Cases:
#   1) d=3, O-Mult (overlap + multinomial)
#   2) d=3, O-Pois (overlap + Poisson)
#   3) d=2, O-Mult
#   4) d=2, O-Pois
#   5) d=1, O-Mult
#   6) d=1, O-Pois
#
# We only store mis_class (1 x 6 matrix).
# nlam=200 for each path.
############################################################

#########################
# Slurm or local array
#########################
ID_index <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if(is.na(ID_index)) ID_index <- 1
num_parallel <- 1000  # e.g. #SBATCH --array=1-250

.libPaths(c("/panfs/jay/groups/27/shenx/zhao1118/AD-real-data/library", .libPaths()))
library(Matrix)
library(glmnet)


#########################
# Real data parameters
#########################
nreps <- 1000
# nreps should be exactly divisible by num_parallel in a typical Slurm scenario.

p_vec = c(100,200,300,400,500)
Model_vec = c(1)  # scheme=1
n_vec = c(100)    # training=100
rate_vec  = c(1,2,4,6,10)          # NEW: rate values

# The total dataset size is 202, so:
#   - training: n in {50, 100}
#   - validation: 150 - n  (i.e., 100 or 50)
#   - test: 52

params <- expand.grid("p" = p_vec , "Model" = Model_vec, "n" = n_vec, "rate" = rate_vec)

############################################################
# 1) Basic design functions: design_U, design_H, etc.
############################################################
design_U <- function(m) {
  U = matrix(0, nrow=m, ncol=m-1)
  for(j in 1:(m-1)){
    U[1:(j+1), j] = c(rep(1,j), -j)/sqrt(j*(j+1))
  }
  return(U)
}

design_H_k <- function(J,k) {
  q=length(J)
  for_ind=combn(q, k)
  rep_full=ncol(for_ind)
  ind_end=rep(1,rep_full)
  for(j in 1:rep_full){
    current=1
    k_vec=for_ind[,j]
    ind_end[j]=prod(J[k_vec]-1)
    if(ind_end[j]==0) next
    for(i in 1:q){
      if(i %in% k_vec){
        current = kronecker(design_U(J[i]), current)
      } else {
        current = kronecker(rep(1,J[i])/sqrt(J[i]), current)
      }
    }
    if(j==1){
      H=current
    } else {
      H=cbind(H,current)
    }
  }
  return(list(H=H, ind_end=ind_end))
}

design_H <- function(J,d) {
  ind_end=1
  H = rep(1, prod(J))/sqrt(prod(J))
  if(d==0){
    return(list(H=t(t(H)), ind_start=1, ind_end=1))
  } else {
    for(k in 1:d){
      output=design_H_k(J,k)
      H=cbind(H, output$H)
      ind_end=c(ind_end, output$ind_end)
    }
    ind_sta=rep(1,length(ind_end))
    for(i in 2:length(ind_end)){
      ind_end[i] = ind_end[i] + ind_end[i-1]
      ind_sta[i] = ind_end[i-1] + 1
    }
    colnames(H)<-NULL
    return(list(H=H, ind_start=ind_sta, ind_end=ind_end))
  }
}

############################################################
# 2) Generating minimal "dummy" data just to get H
#    We'll store: 
#       generate_multinomial_final(...) 
############################################################
cov_decay <- function(p, r){
  mat=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      mat[i,j] = r^abs(i-j)
    }
  }
  return(mat)
}

multinomial_Y_generate <- function(H,X,beta){
  u_mat=H%*%(beta%*%X)
  p_mat=exp(u_mat)/ (rep(1,nrow(u_mat))%*%t(colSums(exp(u_mat))))
  Y=rmultinom(1, size=1, prob=p_mat[,1])
  for(j in 2:ncol(u_mat)){
    Y=cbind(Y, rmultinom(1, size=1, prob=p_mat[,j]))
  }
  return(Y)
}

beta_generate <- function(L_sum,p,ind_end, s, unif_range=c(-2,2), spectral_normalize=5){
  beta=runif(L_sum*p, min=unif_range[1], max=unif_range[2])
  dim(beta)=c(L_sum,p)
  if(s+1<length(ind_end)){
    beta[(ind_end[s+1]+1):L_sum, ]=0
  }
  beta=beta/norm(beta,"2")*spectral_normalize
  return(beta)
}

generate_multinomial_final <- function(J,d,p,n_sample, ind_start_X, ind_end_X,
                                       r=1/2, s=length(J), beta=NULL, R=NULL){
  obj=design_H(J,d)
  H=obj$H
  ind_start_H=obj$ind_start
  ind_end_H=obj$ind_end
  L_sum=ncol(H)
  if(is.null(beta)){
    beta=beta_generate(L_sum, p, ind_end=ind_end_H, s=s)
  }
  if(is.null(R)){
    R=chol(cov_decay(p-1,r))
  }
  X=t(R)%*%matrix(rnorm((p-1)*n_sample), nrow=(p-1))
  X=rbind(rep(1,n_sample),X)
  Y=multinomial_Y_generate(H,X,beta)
  return(list(
    X=X, Y=Y, beta=beta, R=R, 
    H=H, ind_start_H=ind_start_H, ind_end_H=ind_end_H
  ))
}

############################################################
# 3) Real Data: convert_response, make_X
############################################################
one_hot_tensor <- function(j1,j2,j3,J){
  arr=array(0,dim=J)
  arr[j1,j2,j3]=1
  as.vector(arr)
}

map_cats_to_j123 <- function(z1,z2,z3,lev1,lev2,lev3){
  j1=match(z1,lev1)
  j2=match(z2,lev2)
  j3=match(z3,lev3)
  c(j1,j2,j3)
}

convert_response <- function(dataset, J=c(2,2,4)){
  stopifnot(ncol(dataset)>=3)
  lev1=levels(dataset[[1]])
  lev2=levels(dataset[[2]])
  lev3=levels(dataset[[3]])
  nSamples=nrow(dataset)
  
  Y=matrix(0,nrow=prod(J),ncol=nSamples)
  for(i in seq_len(nSamples)){
    z1=dataset[i,1]
    z2=dataset[i,2]
    z3=dataset[i,3]
    j123=map_cats_to_j123(z1,z2,z3,lev1,lev2,lev3)
    Y[,i]=one_hot_tensor(j123[1], j123[2], j123[3], J)
  }
  Y
}

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


############################################################
# 4) Overlapping group-lasso code (the key part)
############################################################
generate_combinations <- function(J, d) {
  out_list=list()
  for(k in 1:d){
    comb_k = combn(length(J), k)
    for(cc in 1:ncol(comb_k)){
      out_list[[length(out_list)+1]] = comb_k[,cc]
    }
  }
  return(out_list)
}

is_subset <- function(set1, set2) {
  all(set1 %in% set2)
}

group_hierarchy <- function(J, d, ind_start_H, ind_end_H){
  H_list_ind = generate_combinations(J, d)
  maxrow = max(ind_end_H)
  group_mat = matrix(FALSE, nrow = maxrow, ncol = length(H_list_ind)+1)
  group_mat[,1] = TRUE
  for(g in 2:ncol(group_mat)){
    start_ = ind_start_H[g]
    end_   = ind_end_H[g]
    group_mat[start_:end_, g] = TRUE
  }
  return(group_mat)
}

project_ball <- function(vec, radius){
  norm_vec = norm(as.matrix(vec), type="f")
  if(norm_vec > radius){
    vec = vec / norm_vec * radius
  }
  return(vec)
}

proximal_new <- function(beta_u, group_H_index, weight, lambda=1, iter_max=50, epsilon=1e-8){
  weight = weight * lambda
  ngroup = ncol(group_H_index)
  group = array(0, dim=c(dim(beta_u), ngroup))
  beta = beta_u
  feasible = rep(1, ngroup)
  for(j in 1:iter_max){
    for(g in 1:ngroup){
      rows_ = which(group_H_index[,g])
      b_sub = beta[rows_, , drop=FALSE]
      b_sub_proj = project_ball(b_sub, weight[g])
      dif = group[rows_, , g] - b_sub_proj
      feasible[g] = sum(abs(dif))
      beta[rows_, ] = beta[rows_, ] + group[rows_, , g]
      group[rows_, , g] = b_sub_proj
      beta[rows_, ] = beta[rows_, ] - group[rows_, , g]
    }
    if(sum(feasible) < epsilon) break
  }
  return(beta)
}

overlapping_group_lasso_proximal <- function(beta_u, group_H_index,
                                             ind_start_X, ind_end_X,
                                             weight, lambda=1,
                                             iter_max=50, epsilon=1e-8){
  for_max = length(ind_start_X)
  final = matrix(0, nrow=nrow(beta_u), ncol=ncol(beta_u))
  for(j in 1:for_max){
    cstart = ind_start_X[j]
    cend   = ind_end_X[j]
    final[, cstart:cend] = proximal_new(
      beta_u = as.matrix(beta_u[, cstart:cend]),
      group_H_index = group_H_index,
      weight = weight[, j],
      lambda = lambda,
      iter_max = iter_max,
      epsilon = epsilon
    )
  }
  return(final)
}

lambda_max_mult <- function(H,X,Y,ind_start_H,ind_end_H,ind_start_X,ind_end_X, w=1){
  beta_zero = matrix(0, nrow=ncol(H), ncol=nrow(X))
  val = grad_ll_mult(beta_zero, H, X, Y)
  lambda_max = max( group_lasso_matrix(val, ind_start_H, ind_end_H, ind_start_X, ind_end_X )/w )
  return(lambda_max)
}

lambda_max_pois <- function(H,X,Y,ind_start_H,ind_end_H,ind_start_X,ind_end_X, w=1){
  beta_zero = matrix(0, nrow=ncol(H), ncol=nrow(X))
  val = grad_ll_pois(beta_zero, H, X, Y)
  lambda_max = max( group_lasso_matrix(val, ind_start_H, ind_end_H, ind_start_X, ind_end_X )/w )
  return(lambda_max)
}

lambda_seq <- function(lambda_max, lambda.min.ratio=0.618, nlambda=100){
  exp(seq(from = log(lambda_max), by = log(lambda.min.ratio), length.out = nlambda))
}

group_lasso_matrix <- function(beta, ind_start_H, ind_end_H, ind_start_X, ind_end_X){
  row_indices <- mapply(seq, ind_start_H, ind_end_H, SIMPLIFY = FALSE)
  col_indices <- mapply(seq, ind_start_X, ind_end_X, SIMPLIFY = FALSE)
  get_norm <- function(r, c) {
    sub_matrix <- beta[r, c, drop = FALSE]
    norm(sub_matrix, type = "F")
  }
  norms <- outer(row_indices, col_indices, Vectorize(get_norm))
  return(norms)
}

grad_ll_mult <- function(beta, H, X, Y){
  expu_mat = exp(H %*% (beta %*% X))
  n_sa = ncol(Y)
  term_mat = rep(1, nrow(expu_mat)) %*% t(colSums(expu_mat))
  p_mat = expu_mat / term_mat
  dis_mat = p_mat - Y
  term = 0
  for(j in 1:n_sa){
    term = term + dis_mat[,j] %*% t(X[,j])
  }
  return(t(H) %*% term / n_sa)
}

ll_mult <- function(beta, H, X, Y){
  n_sa = ncol(Y)
  u_mat = H %*% (beta %*% X)
  C_vec = log(colSums(exp(u_mat)))
  C = sum(colSums(Y) * C_vec)
  return((-sum(Y * u_mat) + C) / n_sa)
}

grad_ll_pois <- function(beta, H, X, Y){
  mu_mat = exp(H %*% (beta %*% X))
  n_sa = ncol(Y)
  dis_mat = mu_mat - Y
  term = 0
  for(j in 1:n_sa){
    term = term + dis_mat[,j] %*% t(X[,j])
  }
  return(t(H) %*% term / n_sa)
}

ll_pois <- function(beta, H, X, Y){
  n_sa = ncol(Y)
  u_mat = H %*% (beta %*% X)
  C = sum(exp(u_mat))
  return((-sum(Y * u_mat) + C) / n_sa)
}

calculate_l2_norms_weight <- function(beta, group_H_index, ind_start_X, ind_end_X, w){
  total_norm = 0
  beta2 = beta^2
  ngroups = ncol(group_H_index)
  nblocks = length(ind_start_X)
  for(i in seq_len(ngroups)){
    rows_ = which(group_H_index[,i])
    for(b in seq_len(nblocks)){
      cstart = ind_start_X[b]
      cend   = ind_end_X[b]
      seg = beta2[rows_, cstart:cend, drop=FALSE]
      total_norm = total_norm + sqrt(sum(seg)) * w[i,b]
    }
  }
  return(total_norm)
}

backtracking_overlap_PGD_single <- function(gradient, loss, beta_init=NULL,
                                            H, X, Y, step_size,
                                            ind_start_H, ind_end_H,
                                            ind_start_X, ind_end_X,
                                            lambda, group_H_index,
                                            shrinking_gamma=0.618,
                                            epsilon=1e-8, max_iter=1000,
                                            min_iter=10, w=1){
  gamma = shrinking_gamma
  eta = step_size
  if(is.null(beta_init)){
    beta_loop = matrix(0, nrow=ncol(H), ncol=nrow(X))
  } else {
    beta_loop = beta_init
  }
  cost = rep(0, max_iter)
  z = beta_loop
  log_loss_z = loss(z, H, X, Y)
  grad = gradient(z, H, X, Y)
  
  # initial attempt
  repeat{
    beta_loop_new = overlapping_group_lasso_proximal(
      z - eta * grad,
      group_H_index,
      ind_start_X, ind_end_X,
      w, lambda = lambda * eta,
      iter_max = 50, epsilon = 1e-8
    )
    D_loop = beta_loop_new - z
    log_loss_beta = loss(beta_loop_new, H, X, Y)
    if(log_loss_beta <= log_loss_z + sum(grad * D_loop) + sum(D_loop^2) / (2 * eta)){
      break
    }
    eta = gamma * eta
  }
  cst_old = log_loss_beta + lambda * calculate_l2_norms_weight(beta_loop_new, group_H_index, ind_start_X, ind_end_X, w)
  for(i in 1:max_iter){
    z = beta_loop_new + i/(i+3) * (beta_loop_new - beta_loop)
    beta_loop = beta_loop_new
    grad = gradient(z, H, X, Y)
    log_loss_z = loss(z, H, X, Y)
    repeat{
      beta_loop_new = overlapping_group_lasso_proximal(
        z - eta * grad,
        group_H_index,
        ind_start_X, ind_end_X,
        w, lambda = lambda * eta,
        iter_max = 50, epsilon = 1e-8
      )
      D_loop = beta_loop_new - z
      log_loss_beta = loss(beta_loop_new, H, X, Y)
      if(log_loss_beta <= log_loss_z + sum(grad * D_loop) + sum(D_loop^2) / (2 * eta)){
        break
      }
      eta = gamma * eta
    }
    cost[i] = log_loss_beta + lambda * calculate_l2_norms_weight(beta_loop_new, group_H_index, ind_start_X, ind_end_X, w)
    cst_new = cost[i]
    if(cst_old - cst_new < epsilon && sum(D_loop^2) < epsilon){
      if(i > min_iter){
        beta_final = beta_loop_new
        break
      }
    }
    if(i == max_iter){
      beta_final = beta_loop_new
    }
    cst_old = cst_new
  }
  training_loss = cost
  return(list(beta_final = beta_final, training_loss = training_loss[1:i], eta_end = eta))
}

backtracking_overlap_PGD_path <- function(gradient, loss, H, X, Y,
                                          step_size, ind_start_H, ind_end_H,
                                          ind_start_X, ind_end_X,
                                          lam_seq, shrinking_gamma=0.618,
                                          epsilon=1e-8, group_H_index,
                                          max_iter=100, min_iter=2, w=1){
  gamma = shrinking_gamma
  nlambda = length(lam_seq)
  # first
  ot = backtracking_overlap_PGD_single(
    gradient, loss, beta_init = NULL,
    H, X, Y, step_size,
    ind_start_H, ind_end_H,
    ind_start_X, ind_end_X,
    lam_seq[1], group_H_index,
    shrinking_gamma = gamma,
    epsilon = epsilon,
    max_iter = max_iter,
    min_iter = min_iter,
    w = w
  )
  eta = ot$eta_end
  beta_tensor = array(0, dim = c(ncol(H), nrow(X), nlambda))
  beta_tensor[,,1] = ot$beta_final
  
  lam_ind = 1
  for(lambda_cur in lam_seq[-1]){
    lam_ind = lam_ind + 1
    ot = backtracking_overlap_PGD_single(
      gradient, loss,
      beta_init = beta_tensor[,,lam_ind-1],
      H, X, Y, eta / gamma,
      ind_start_H, ind_end_H,
      ind_start_X, ind_end_X,
      lambda_cur, group_H_index,
      shrinking_gamma = gamma,
      epsilon = epsilon,
      max_iter = max_iter,
      min_iter = min_iter,
      w = w
    )
    eta = ot$eta_end
    beta_tensor[,,lam_ind] = ot$beta_final
  }
  return(list(beta_tensor = beta_tensor))
}

############################################################
# 5) Validation + misclassification
############################################################
beta_validation <- function(beta_tensor, lam_seq, loss, H,
                            X_validation, Y_validation,
                            n_train, plot_TF=FALSE){
  nlambda = dim(beta_tensor)[3]
  cross_entropy = rep(0, nlambda)
  for(i in 1:nlambda){
    cross_entropy[i] = loss(beta_tensor[,,i], H, X_validation, Y_validation)
  }
  min_cross_entropy = min(cross_entropy)
  min_index = which.min(cross_entropy)
  # minimal approach => no 1SE calc
  se_val = 0
  cross_entropy_threshold = min_cross_entropy + se_val
  valid_lams = lam_seq[cross_entropy <= cross_entropy_threshold]
  lambda_1se = max(valid_lams)
  lambda_1se_index = which(lam_seq == lambda_1se)[1]
  if(plot_TF){
    plot(log(lam_seq), cross_entropy, type='l', xlab='log(lambda)', ylab='CE')
    abline(h=cross_entropy_threshold, col='red', lty=2)
    points(log(lam_seq[min_index]), cross_entropy[min_index], col='darkgreen', pch=19)
  }
  return(list(lambda_min_index = min_index, lambda_1se_index = lambda_1se_index))
}

calculate_misclassification_rate <- function(X_test, Y_test, beta_est, H, theta=NULL){
  if(is.null(theta)){
    exp_mat = exp(H %*% (beta_est %*% X_test))
  } else {
    exp_mat = exp(theta %*% X_test)
  }
  predicted_prob = sweep(exp_mat, 2, colSums(exp_mat), FUN='/')
  predicted_classes = max.col(t(predicted_prob))
  true_classes = max.col(t(Y_test))
  return(sum(predicted_classes != true_classes) / ncol(Y_test))
}

weight_generate <- function(J, d, ind_end_H, ind_end_X, weight_deep=NULL, weight_groupX=NULL, scheme_X=0){
  if(is.null(weight_groupX)){
    weight_groupX = rep(1, length(ind_end_X))
  }
  if(is.null(weight_deep)){
    weight_deep_vec = rep(1, length(ind_end_H))
  } else {
    q = length(J)
    weight_deep_vec = weight_deep[1]
    for(i in 1:d){
      weight_deep_vec = c(weight_deep_vec, rep(weight_deep[i+1], choose(q, i)))
    }
  }
  weight_mat = outer(weight_deep_vec, weight_groupX, "*")
  if(scheme_X == 1){
    weight_mat[1:ind_end_H[1], 1] = 0
  }
  if(scheme_X == 2){
    weight_mat[ind_end_H[length(J)+1], 1] = 0
  }
  return(weight_mat)
}


custom_weight_generate <- function(J, d, ind_end_H, ind_end_X, rate, weight_groupX = NULL) {
  # Create weight_deep vector: length = d+1, with elements 1, rate, rate^2, ..., rate^d.
  weight_deep <- c(1, rate^(0:(d-1)))
  
  # Call the original weight_generate function with our weight_deep.
  weight_mat <- weight_generate(J, d, ind_end_H, ind_end_X, weight_deep = weight_deep, weight_groupX = weight_groupX)
  
  return(weight_mat)
}
############################################################
# 6) Main Overlapping G-Lasso Real Data Loop
############################################################
dataset = read.csv("/panfs/jay/groups/27/shenx/zhao1118/AD-real-data/dataset.csv")
dfReal = dataset
dfReal[[1]] = as.factor(dfReal[[1]])
dfReal[[2]] = as.factor(dfReal[[2]])
dfReal[[3]] = as.factor(dfReal[[3]])

for(ind_loop in seq_len(nrow(params))){
  p     <- params[ind_loop, "p"]
  Model <- params[ind_loop, "Model"]
  n     <- params[ind_loop, "n"]
  rate  <- params[ind_loop, "rate"]   # NEW: retrieve the rate
  
  # For each replication
  for(rep_i in ID_index + (0:(nreps/num_parallel - 1)) * num_parallel){
    cat(sprintf(">>> p=%d, Model=%d, n=%d, rep_i=%d\n", p, Model, n, rep_i))
    
    # 1) Convert to Y, X
    J = c(2,2,4)
    Y_full = convert_response(dfReal, J)
    X_full = make_X(dfReal, addIntercept=TRUE, colSample=p, seed=2024)
    
    n_train = n
    n_validation = 150 - n
    n_test = 52
    stopifnot(n_train + n_validation + n_test == 202)
    
    # Shuffle
    set.seed(2024 + rep_i)
    perm = sample(202)
    X_full_perm = X_full[, perm]
    Y_full_perm = Y_full[, perm]
    
    X_train = X_full_perm[, 1:n_train]
    Y_train = Y_full_perm[, 1:n_train]
    X_validation = X_full_perm[, (n_train + 1):(n_train + n_validation)]
    Y_validation = Y_full_perm[, (n_train + 1):(n_train + n_validation)]
    X_test = X_full_perm[, (n_train + n_validation + 1):(n_train + n_validation + n_test)]
    Y_test = Y_full_perm[, (n_train + n_validation + 1):(n_train + n_validation + n_test)]
    
    # We'll do 6 cases => mis_class is 1 x 6
    mis_class = matrix(0, nrow=1, ncol=6)
    # --- NEW: Initialize a matrix to hold the joint cross entropy for each case ---
    joint_CE_mat = matrix(0, nrow=1, ncol=6)
    
    # We define nlam = 200
    nlam = 200
    epsilon = 1e-4
    
    # Cases:
    # 1 => d=3, O-Mult
    # 2 => d=3, O-Pois
    # 3 => d=2, O-Mult
    # 4 => d=2, O-Pois
    # 5 => d=1, O-Mult
    # 6 => d=1, O-Pois
    for(caseID in 1:6){
      if(caseID %in% c(1,2)){
        d_val = 3
      } else if(caseID %in% c(3,4)){
        d_val = 2
      } else {
        d_val = 1
      }
      
      # "Dummy" to get H_, ind_start_H_, etc.
      dummy = generate_multinomial_final(J=J, d=d_val, p=p,
                                         n_sample=2, ind_start_X=1:2, ind_end_X=1:2,
                                         r=1/2, s=length(J), beta=NULL, R=NULL)
      H_ = dummy$H
      ind_start_H_ = dummy$ind_start_H
      ind_end_H_ = dummy$ind_end_H
      
      # Build group hierarchy
      group_H_idx = group_hierarchy(J, d_val, ind_start_H_, ind_end_H_)
      
      # Define the column block as c(1, p)
      ind_end_X_ = c(1, p)
      ind_start_X_ = c(1, 1 + ind_end_X_[-length(ind_end_X_)])
      
      # *** NEW: Use custom_weight_generate instead of manual weight generation ***
      weight_mat <- custom_weight_generate(J, d_val, ind_end_H_, ind_end_X_, rate = rate)
      
      if(caseID %in% c(1,3,5)){
        # Overlapping + Mult
        cat("Case", caseID, ": d=", d_val, ", Overlapping Mult\n")
        lam_max_val = lambda_max_mult(H_, X_train, Y_train,
                                      ind_start_H_, ind_end_H_,
                                      ind_start_X_, ind_end_X_, w=weight_mat)
        lam_seq = lambda_seq(lam_max_val, lambda.min.ratio=0.90, nlambda=nlam)
        eta_ = prod(J) * n_train / norm(X_train, "2")^2
        out_path = backtracking_overlap_PGD_path(
          gradient = grad_ll_mult, loss = ll_mult,
          H = H_, X = X_train, Y = Y_train,
          step_size = eta_,
          ind_start_H = ind_start_H_,
          ind_end_H = ind_end_H_,
          ind_start_X = ind_start_X_,
          ind_end_X = ind_end_X_,
          lam_seq = lam_seq,
          shrinking_gamma = 0.618,
          epsilon = epsilon,
          group_H_index = group_H_idx,
          max_iter = 100,
          min_iter = 2,
          w = weight_mat
        )
        beta_tensor = out_path$beta_tensor
        val_out = beta_validation(beta_tensor, lam_seq, ll_mult, H_,
                                  X_validation, Y_validation,
                                  n_train = n, plot_TF = FALSE)
        best_idx = val_out$lambda_min_index
        mis_class[1, caseID] = calculate_misclassification_rate(
          X_test, Y_test, beta_est = beta_tensor[,,best_idx], H = H_
        )
        # --- NEW: Compute joint CE on test set using ll_mult ---
        joint_CE_mat[1, caseID] = ll_mult(beta_tensor[,,best_idx], H_, X_test, Y_test)
        
      } else {
        # Overlapping + Poisson (cases 2,4,6)
        cat("Case", caseID, ": d=", d_val, ", Overlapping Pois\n")
        lam_max_val = lambda_max_pois(H_, X_train, Y_train,
                                      ind_start_H_, ind_end_H_,
                                      ind_start_X_, ind_end_X_, w=weight_mat)
        lam_seq = lambda_seq(lam_max_val, lambda.min.ratio=0.85, nlambda=nlam)
        eta_ = n_train / norm(X_train, "2")^2
        out_path = backtracking_overlap_PGD_path(
          gradient = grad_ll_pois, loss = ll_pois,
          H = H_, X = X_train, Y = Y_train,
          step_size = eta_,
          ind_start_H = ind_start_H_,
          ind_end_H = ind_end_H_,
          ind_start_X = ind_start_X_,
          ind_end_X = ind_end_X_,
          lam_seq = lam_seq,
          shrinking_gamma = 0.618,
          epsilon = epsilon,
          group_H_index = group_H_idx,
          max_iter = 100,
          min_iter = 2,
          w = weight_mat
        )
        beta_tensor = out_path$beta_tensor
        val_out = beta_validation(beta_tensor, lam_seq, ll_mult, H_,
                                  X_validation, Y_validation,
                                  n_train = n, plot_TF = FALSE)
        best_idx = val_out$lambda_min_index
        mis_class[1, caseID] = calculate_misclassification_rate(
          X_test, Y_test, beta_est = beta_tensor[,,best_idx], H = H_
        )
        # --- NEW: Compute joint CE on test set using ll_pois ---
        joint_CE_mat[1, caseID] = ll_pois(beta_tensor[,,best_idx], H_, X_test, Y_test)
      }
    } # end for caseID in 1:6
    
    # Save results => now we store mis_class and joint_CE for 6 cases
    savename <- paste0("/panfs/jay/groups/27/shenx/zhao1118/AD-real-data/data_overlap/Model",
                       Model, "_p", p, "_n", n, "_rate", rate, "_repID", rep_i, ".RDS")
    Results <- list(
      Misclassification_rate = mis_class,  # 1 x 6 matrix
      Joint_CE = joint_CE_mat,             # --- NEW: Joint CE matrix, 1 x 6
      Model = Model,
      p = p,
      n = n,
      rate = rate,
      rep_ID = rep_i
    )
    saveRDS(Results, file = savename)
    
  } # end for rep_i
} # end for ind_loop in seq_len(nrow(params))

