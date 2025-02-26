###############################################
#  Real Data Pipeline: Multinomial Log-Linear
###############################################

# -------------------------------
# Slurm-related (if needed)
# -------------------------------
ID_index <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(is.na(ID_index)) ID_index <- 1  # Fallback if not running on Slurm
num_parallel=1000  # e.g., if you used #SBATCH --array=1-25

## At the start of your R script:

#if (!requireNamespace("glmnet", quietly = TRUE)) {
#  install.packages("glmnet", repos = "https://cloud.r-project.org/")
#}
.libPaths(c("/panfs/jay/groups/27/shenx/zhao1118/AD-real-data/library", .libPaths()))
library(Matrix)
library(glmnet)


# --------------------------------------------
# Set model parameters 
# --------------------------------------------
nreps <- 1000
# nreps should be exactly divisible by num_parallel in a typical Slurm scenario.

# Real data scenario:
p_vec      =  c(100,200,300,400,500)    # columns (after the first 3 are used for Y)
Model_vec  = c(1)        # only scheme=1
n_vec      = c(100)   # training sizes we want to try
rate_vec  = c(1,2,4,6,10)          # NEW: rate values

# The total dataset size is 202, so:
#   - training: n in {50, 100}
#   - validation: 150 - n  (i.e., 100 or 50)
#   - test: 52

params <- expand.grid("p" = p_vec , "Model" = Model_vec, "n" = n_vec, "rate" = rate_vec)

######### Basic functions (Design U, H_k, etc.)
design_U <- function(m) {
  U = matrix(0, nrow = m, ncol = m-1 )
  for(j in 1:m-1){
    U[1:(j+1),j]=c(rep(1,j), - j )/sqrt(j*(j+1))
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
    if( ind_end[j]==0 ){next}
    for(i in 1:q ){
      if( i %in% k_vec ){
        current = kronecker( design_U(J[i]) , current )
      }else{
        current = kronecker( rep(1,J[i])/sqrt(J[i]) , current )
      }
    }
    if(j==1){
      H=current
    }else{
      H=cbind(H,current)
    }
  }
  return( list(H=H, ind_end=ind_end ) )
}

design_H <- function(J,d) {
  ind_end=1
  H=rep(1,prod(J))/sqrt(prod(J))
  if(d==0){
    return( list(H=t(t(H)), ind_start=1, ind_end=1 ) )
  }else{
    for(k in 1:d ){
      output=design_H_k(J, k )  
      H=cbind(H, output$H )
      ind_end=c(ind_end, output$ind_end )
    }
    ind_sta=rep(1,length(ind_end))
    for(i in 2:length(ind_end)){
      ind_end[i] = ind_end[i] + ind_end[i-1] 
      ind_sta[i] = ind_end[i-1] + 1
    }
    colnames(H)<-NULL
    return( list(H=H, ind_start=ind_sta, ind_end=ind_end ) )
  }
}

beta_generate <- function(L_sum,p,ind_end,s,unif_range=c(-2,2), spectral_normalize=5){
  beta=runif(L_sum*p, min=unif_range[1], max=unif_range[2])
  dim(beta)=c(L_sum,p)
  if( s+1< length(ind_end) ){
    beta[(ind_end[s+1]+1):L_sum,]=0
  }
  beta=beta/norm(beta, type = "2")*spectral_normalize
  return(beta)
}

beta_generate_scheme <- function(J,d,L_sum,p,ind_start_H,ind_end_H,unif_range=c(-2,2),scheme=1 ){
  beta=matrix(0,nrow=L_sum, ncol=p  )
  q=length(J)
  col_ind=c(1,sort(sample(2:p,2))) 
  ii=1
  if(scheme==1 ){
    # just fill in first-order terms (example)
    for_ind=combn(q, 1)
    rep_full=ncol(for_ind)
    for(j in 1:rep_full){
      ii=ii+1
      rows_sel <- ind_start_H[ii]:ind_end_H[ii]
      beta[ rows_sel, col_ind ] = matrix(
        runif(length(rows_sel)*length(col_ind), min=unif_range[1], max=unif_range[2]),
        ncol=length(col_ind)
      )
    }
  }
  return(beta)
}

cov_decay <- function(p , r) {
  matrix <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p ) {
      matrix[i, j] <- r^abs(i - j)
    }
  }
  return(matrix)
}

poisson_Y_generate <- function(H, X, beta) {
  u_mat = H %*% beta %*% X
  Y=rpois(length(u_mat), exp(u_mat))
  dim(Y)=dim(u_mat)
  return(Y)
}

multinomial_Y_generate <- function(H, X, beta) {
  u_mat = H %*% beta %*% X
  p_mat=exp(u_mat) / (rep(1,nrow(u_mat))%*%t(colSums(exp(u_mat))) )
  Y=rmultinom( 1, size=1, prob=p_mat[,1])
  for(j in 2:ncol(u_mat) ){
    Y=cbind(Y,rmultinom(1, size=1, prob=p_mat[,j]))
  }
  return(Y)
}

ll_pois <- function(beta,H,X,Y){
  n_sa=ncol(Y)
  u_mat=H%*% (beta %*% X)
  C=sum(exp(u_mat))
  return( (-sum(Y*u_mat) + C)/n_sa)
}

grad_ll_pois <- function(beta,H,X,Y){
  mu_mat=exp(H%*% (beta %*% X) )
  n_sa=ncol(Y)
  dis_mat=mu_mat - Y
  term=0
  for(j in 1:n_sa){
    term=term+dis_mat[,j]%*%t(X[,j])
  }
  return(t(H)%*%term/n_sa)
}

ll_mult <- function(beta,H,X,Y){
  n_sa=ncol(Y)
  u_mat=H%*% (beta %*% X)
  C_vec=log(colSums(exp(u_mat)))
  C=sum(colSums(Y)*C_vec)
  return( (-sum(Y*u_mat) + C)/n_sa )
}

ll_mult_se <- function(beta,H,X,Y){
  n_sa=ncol(Y)
  u_mat=H%*% (beta %*% X)
  C_vec=log(colSums(exp(u_mat)))
  C=sum(colSums(Y)*C_vec)
  ll_vec=-colSums(Y*u_mat) +colSums(Y)*C_vec
  return( list(mean_ll=mean(ll_vec), se_ll=sd(ll_vec)  ) )
}

grad_ll_mult <- function(beta,H,X,Y){
  expu_mat=exp(H%*% (beta %*% X) )
  n_sa=ncol(Y)
  term_mat=(rep(1,nrow(expu_mat))%*%t(colSums(expu_mat)))
  p_mat=expu_mat/term_mat
  dis_mat=p_mat-Y
  term=0
  for(j in 1:n_sa){
    term=term+dis_mat[,j]%*%t(X[,j])
  }
  return(t(H)%*%term/n_sa)
}

group_lasso <- function(beta, ind_start_H,ind_end_H,ind_start_X,ind_end_X, w=1){
  row_indices <- mapply(seq, ind_start_H, ind_end_H, SIMPLIFY = FALSE)
  col_indices <- mapply(seq, ind_start_X, ind_end_X, SIMPLIFY = FALSE)
  get_norm <- function(r, c) {
    sub_matrix <- beta[r, c, drop = FALSE]
    norm(sub_matrix, type = "F")
  }
  norms <- outer(row_indices, col_indices, Vectorize(get_norm))
  return( sum(norms*w) )
}

group_lasso_proximal <- function(beta_z, w, lambda, ind_start_H, ind_end_H, ind_start_X, ind_end_X) {
  row_indices <- mapply(seq, ind_start_H, ind_end_H, SIMPLIFY = FALSE)
  col_indices <- mapply(seq, ind_start_X, ind_end_X, SIMPLIFY = FALSE)
  if (length(w) == 1){
    ww=rep(w,length(row_indices)*length(col_indices))
    dim(ww)=c(length(row_indices), length(col_indices))
    w=ww
  }
  process_submatrix <- function(r, c, i, j) {
    b_sub <- beta_z[r, c]
    frob_norm <- sqrt(sum(b_sub^2))
    if (frob_norm > lambda * w[i, j]) {
      return((1 - lambda * w[i, j] / frob_norm) * b_sub)
    } else {
      return(matrix(0, nrow = length(r), ncol = length(c)))
    }
  }
  final <- matrix(0, nrow = nrow(beta_z), ncol = ncol(beta_z))
  for (i in 1:length(row_indices)) {
    for (j in 1:length(col_indices)) {
      final[row_indices[[i]], col_indices[[j]]] <- 
        process_submatrix(row_indices[[i]], col_indices[[j]], i, j)
    }
  }
  return(final)
}

generate_multinomial_final<-function(J,d,p,n_sample, ind_start_X,ind_end_X,r=1/2,s=length(J),beta=NULL,R=NULL){
  obj=design_H(J,d)
  H=obj$H
  ind_start_H=obj$ind_start
  ind_end_H=obj$ind_end
  L_sum=ncol(H)
  if (is.null(beta)) {
    beta=beta_generate(L_sum,p,ind_end=ind_end_H,s)
  }
  if (is.null(R)) {
    R = chol(cov_decay(p-1, r))
  }
  X = t(R) %*% matrix(rnorm((p-1)*n_sample), nrow=p-1 )
  X = rbind(rep(1,n_sample), X)
  Y=multinomial_Y_generate(H, X, beta)
  return(list(X=X, Y=Y, beta=beta, R=R, H=H, ind_start_H=ind_start_H, ind_end_H=ind_end_H ))
}

group_lasso_matrix <- function(beta, ind_start_H,ind_end_H,ind_start_X,ind_end_X ){
  row_indices <- mapply(seq, ind_start_H, ind_end_H, SIMPLIFY = FALSE)
  col_indices <- mapply(seq, ind_start_X, ind_end_X, SIMPLIFY = FALSE)
  get_norm <- function(r, c) {
    sub_matrix <- beta[r, c, drop = FALSE]
    norm(sub_matrix, type = "F")
  }
  norms <- outer(row_indices, col_indices, Vectorize(get_norm))
  return(norms)
}

lambda_max_mult <- function(H,X,Y,ind_start_H,ind_end_H,ind_start_X,ind_end_X, w=1){
  beta_zero = matrix(0, nrow=ncol(H), ncol=nrow(X))
  val = grad_ll_mult(beta_zero,H,X,Y)
  lambda_max = max( group_lasso_matrix(val,ind_start_H,ind_end_H,ind_start_X,ind_end_X )/w )
  return(lambda_max)
}

lambda_max_pois <- function(H,X,Y,ind_start_H,ind_end_H,ind_start_X,ind_end_X, w=1){
  beta_zero = matrix(0, nrow=ncol(H), ncol=nrow(X))
  val = grad_ll_pois(beta_zero,H,X,Y)
  lambda_max = max( group_lasso_matrix(val,ind_start_H,ind_end_H,ind_start_X,ind_end_X )/w )
  return(lambda_max)
}

lambda_seq<-function(lambda_max,lambda.min.ratio=0.618,nlambda=100){
  exp(seq(from = log(lambda_max),  by = log(lambda.min.ratio), length.out = nlambda ))
}

weight_generate <- function( J, d, ind_end_H, ind_end_X, weight_deep=NULL, weight_groupX=NULL ){
  if ( is.null(weight_groupX) ){
    weight_groupX=rep(1,length(ind_end_X))
  }
  if ( is.null(weight_deep) ){
    weight_deep_vec=rep(1,length(ind_end_H))
  }else{
    q=length(J)
    weight_deep_vec=weight_deep[1]
    for(i in 1:d){
      weight_deep_vec=c(weight_deep_vec, rep(weight_deep[i+1], choose(q, i)) )
    }
  }
  weight_mat = outer(weight_deep_vec, weight_groupX, "*")
  return(weight_mat)
}

custom_weight_generate <- function(J, d, ind_end_H, ind_end_X, rate, weight_groupX = NULL) {
  # Create weight_deep vector: length = d+1, with elements 1, rate, rate^2, ..., rate^d.
  weight_deep <- c(1, rate^(0:(d-1)))
  
  # Call the original weight_generate function with our weight_deep.
  weight_mat <- weight_generate(J, d, ind_end_H, ind_end_X, weight_deep = weight_deep, weight_groupX = weight_groupX)
  
  return(weight_mat)
}


backtracking_PGD_single<-function(gradient,loss,beta_init=NULL,H,X,Y,
                                  step_size,ind_start_H,ind_end_H,
                                  ind_start_X,ind_end_X,lambda,
                                  shrinking_gamma=0.618,epsilon=1e-8,
                                  max_iter=1000,min_iter=10, w=1){
  gamma=shrinking_gamma
  eta=step_size
  if(is.null(beta_init)){
    beta_loop=matrix(0, nrow=ncol(H), ncol=nrow(X))
  } else {
    beta_loop=beta_init
  }
  cost=rep(0,max_iter)
  z=beta_loop
  log_loss_z=loss(z, H,X,Y)
  grad=gradient(z, H,X,Y)
  repeat{
    beta_loop_new=group_lasso_proximal(z-eta*grad , w=w, lambda*eta,
                                       ind_start_H, ind_end_H,
                                       ind_start_X, ind_end_X)
    D_loop = beta_loop_new - z
    log_loss_beta =loss(beta_loop_new, H, X, Y)
    if ( log_loss_beta <= log_loss_z + sum(grad * D_loop)+sum(D_loop^2) /(eta * 2 ) ){
      break
    }
    eta = gamma * eta
  }
  cst_old=log_loss_beta + lambda*group_lasso(beta_loop_new ,ind_start_H,ind_end_H,ind_start_X,ind_end_X, w=w) 
  for (i in 1:max_iter) {
    z = beta_loop_new + i/(i+3)* ( beta_loop_new-beta_loop)
    beta_loop = beta_loop_new 
    grad=gradient(z, H,X,Y)
    log_loss_z=loss(z, H, X, Y)
    repeat{
      beta_loop_new=group_lasso_proximal(z-eta*grad, w=w, lambda*eta,
                                         ind_start_H, ind_end_H, ind_start_X, ind_end_X)
      D_loop = beta_loop_new - z
      log_loss_beta =loss(beta_loop_new, H, X, Y)
      if ( log_loss_beta <= log_loss_z + sum(grad * D_loop)+sum(D_loop^2) /(eta * 2 ) ) {
        break
      }
      eta = gamma * eta
    }
    cost[i]=log_loss_beta+lambda*group_lasso(beta_loop_new, ind_start_H,ind_end_H,ind_start_X,ind_end_X, w=w)
    cst_new=cost[i]
    if (cst_old - cst_new < epsilon & sum(D_loop^2) < epsilon ){
      if (i>min_iter){
        beta_final=beta_loop_new
        break
      }
    }
    if (i == max_iter){
      beta_final=beta_loop_new
    }
    cst_old=cst_new
  }
  training_loss=cost[1:i]
  return(list(beta_final=beta_final,training_loss=training_loss, eta_end=eta))
}

backtracking_PGD_path<-function(gradient,loss,H,X,Y,step_size,ind_start_H,ind_end_H,
                                ind_start_X,ind_end_X,lam_seq,
                                shrinking_gamma=0.618,epsilon=1e-8,
                                max_iter=100,min_iter=10,w=1){
  gamma=shrinking_gamma
  nlambda=length(lam_seq)
  ot=backtracking_PGD_single(gradient,loss,beta_init=NULL,H,X,Y,
                             step_size,ind_start_H,ind_end_H,ind_start_X,ind_end_X,
                             lam_seq[1],gamma,epsilon,max_iter,min_iter,w)  
  eta=ot$eta_end
  beta_tensor=array(0, dim=c(ncol(H), nrow(X), nlambda))
  beta_tensor[,,1]=ot$beta_final
  for(lam_ind in 2:nlambda){
    ot=backtracking_PGD_single(gradient,loss,
                               beta_init=beta_tensor[,,lam_ind-1],
                               H,X,Y,eta/gamma,
                               ind_start_H,ind_end_H,ind_start_X,ind_end_X,
                               lam_seq[lam_ind],gamma,epsilon,max_iter,min_iter,w)  
    eta=ot$eta_end
    beta_tensor[,,lam_ind]=ot$beta_final
  }
  return( list( beta_tensor=beta_tensor ) )
}

beta_validation<-function(beta_tensor,lam_seq,loss,H,X_validation,Y_validation,n_train,plot_TF=FALSE){
  nlambda=dim(beta_tensor)[3]
  cross_entropy=rep(0, nlambda)
  for(i in 1:nlambda){
    cross_entropy[i]=loss(beta_tensor[,,i],H,X_validation,Y_validation)
  }
  min_cross_entropy <- min(cross_entropy)
  min_index <- which.min(cross_entropy)
  se_val <- ll_mult_se(beta_tensor[,,min_index], H, X_validation, Y_validation)$se_ll / sqrt(n_train)
  cross_entropy_threshold <- min_cross_entropy + se_val
  valid_lambdas <- lam_seq[ cross_entropy <= cross_entropy_threshold ]
  lambda_1se <- max(valid_lambdas)
  lambda_1se_index = min(which( cross_entropy <= cross_entropy_threshold ))
  if(plot_TF == TRUE){
    plot(log(lam_seq), cross_entropy, type = "l", xlab = "log(λ)", ylab = "Cross Entropy")
    abline(h = cross_entropy_threshold, col = "red", lty = 2)
    points(log(lambda_1se), min(cross_entropy[lam_seq == lambda_1se]), col = "blue", pch = 3)
    text(log(lambda_1se), min(cross_entropy[lam_seq == lambda_1se]), labels ="λ.1se" , pos =3, col = "blue")
    points(log(lam_seq[min_index]), min_cross_entropy, col = "darkgreen", pch =4)
    text(log(lam_seq[min_index]), min_cross_entropy, labels = "λ.min " , pos = 2,col = "darkgreen")
  }
  return(list(lambda_min_index=min_index, lambda_1se_index=lambda_1se_index ))
}

calculate_misclassification_rate <- function(X_test, Y_test, beta_est, H, theta=NULL) {
  # If 'theta' is NULL, use beta_est & H. Otherwise use 'theta' alone.
  if( is.null(theta) ){
    exp_mat <- exp(H %*% beta_est %*% X_test)
  }else{
    exp_mat <- exp(theta %*% X_test)
  }
  predicted_probabilities <- sweep(exp_mat, 2, colSums(exp_mat), FUN="/")
  predicted_classes <- max.col(t(predicted_probabilities))
  true_classes <- max.col(t(Y_test))
  num_misclassified <- sum(predicted_classes != true_classes)
  total_observations <- ncol(Y_test)
  misclassification_rate <- num_misclassified / total_observations
  return(misclassification_rate)
}

signal_strengthen <- function(beta, delta=2, sigma=2) {
  result <- matrix(0, nrow = nrow(beta), ncol = ncol(beta))
  for (i in 1:nrow(beta)) {
    for (j in 1:ncol(beta)) {
      if (beta[i, j] > 0) {
        result[i, j] <- (beta[i, j] + delta) / sigma
      } else if (beta[i, j] < 0) {
        result[i, j] <- (beta[i, j] - delta) / sigma
      } else {
        result[i, j] <- 0
      }
    }
  }
  return(result)
}

#### For real data conversion
one_hot_tensor <- function(j1, j2, j3, J) {
  arr <- array(0, dim = J)
  arr[j1, j2, j3] <- 1
  as.vector(arr)
}

map_cats_to_j123 <- function(z1, z2, z3, lev1, lev2, lev3) {
  j1 <- match(z1, lev1)
  j2 <- match(z2, lev2)
  j3 <- match(z3, lev3)
  c(j1, j2, j3)
}

convert_response <- function(dataset, J = c(2,2,4)) {
  stopifnot(ncol(dataset) >= 3)
  lev1 <- levels(dataset[[1]])
  lev2 <- levels(dataset[[2]])
  lev3 <- levels(dataset[[3]])
  nSamples <- nrow(dataset)
  
  Y <- matrix(0, nrow=prod(J), ncol=nSamples)
  for(i in seq_len(nSamples)) {
    z1 <- dataset[i,1]
    z2 <- dataset[i,2]
    z3 <- dataset[i,3]
    j123 <- map_cats_to_j123(z1, z2, z3, lev1, lev2, lev3)
    Y[, i] <- one_hot_tensor(j123[1], j123[2], j123[3], J)
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


# Additional THETA-versions for methods (4) and (5)
### old
ll_mult_theta <- function(theta, X, Y){
  n_sa = ncol(Y)
  u_mat = theta %*% X
  C_vec = log(colSums(exp(u_mat)))
  C = sum(colSums(Y)*C_vec)
  return((-sum(Y*u_mat) + C)/n_sa)
}
### new
debug_ll_mult_theta <- function(theta, X, Y) {
  cat("\n------------------\n")
  cat("DEBUG in ll_mult_theta:\n")
  cat("dim(theta) =", dim(theta), "\n")
  cat("dim(X)     =", dim(X), "\n")
  cat("dim(Y)     =", dim(Y), "\n")
  
  # 1) Check if the columns of 'theta' match the rows of 'X'
  #    =>  (nrow(theta) might be #classes, ncol(theta) = #features)
  #    =>  X should be (#features) x (#samples)
  if(ncol(theta) != nrow(X)) {
    stop(paste0(
      "Non-conformable: ncol(theta)=", ncol(theta), 
      " != nrow(X)=", nrow(X)
    ))
  }
  
  # 2) Check if the rows of 'Y' match the #classes in 'theta'
  if(nrow(Y) != nrow(theta)) {
    stop(paste0(
      "Non-conformable: nrow(Y)=", nrow(Y),
      " != nrow(theta)=", nrow(theta)
    ))
  }
  
  # 3) Check if the columns of X, Y match
  #    => both should be #samples
  if(ncol(X) != ncol(Y)) {
    stop(paste0(
      "X, Y have different #samples: ncol(X)=", ncol(X), 
      " != ncol(Y)=", ncol(Y)
    ))
  }
  
  # If all checks pass, we do the log-likelihood
  n_sa = ncol(Y)
  u_mat = theta %*% X   # => (#classes) x (#samples)
  # Then we sum exponentiated 'u_mat' to get normalizing constants
  C_vec = log(colSums(exp(u_mat))) 
  # Weighted by how many samples in each column (or Y colSums)
  C = sum(colSums(Y) * C_vec)
  
  # final negative log-likelihood
  ll_val = (-sum(Y * u_mat) + C) / n_sa
  
  cat("DEBUG: -sum(Y*u_mat) =", -sum(Y*u_mat), 
      ", C =", C, 
      ", n_sa =", n_sa, 
      "\n")
  cat("DEBUG: returning ll_val =", ll_val, "\n")
  cat("------------------\n\n")
  
  return(ll_val)
}



grad_ll_mult_theta <- function(theta, X, Y){
  expu_mat = exp(theta %*% X)
  n_sa=ncol(Y)
  term_mat = (rep(1,nrow(expu_mat))%*%t(colSums(expu_mat)))
  p_mat = expu_mat/term_mat
  dis_mat = p_mat - Y
  term=0
  for(j in 1:n_sa){
    term=term + dis_mat[,j] %*% t(X[,j])
  }
  return(term/n_sa)
}

# Prox for row-lasso or row-group-lasso on THETA
lasso_proximal_theta <- function(Z, lambda) {
  # Soft threshold each element
  theta <- sign(Z)*pmax(abs(Z)-lambda, 0)
  return(theta)
}

group_lasso_proximal_theta <- function(theta_z, lambda) {
  frob_norms <- sqrt(rowSums(theta_z^2))
  scaling_factors <- pmax(1 - lambda/frob_norms, 0)
  result <- sweep(theta_z, 1, scaling_factors, '*')
  return(result)
}

# We'll define separate functions to run APGD for THETA
accelerated_proximal_gradient_descent_fixed_eta <- function(X, Y, lambda, step_size, epsilon=1e-8,
                                                            max_iter=1000, min_iter=10, theta_init=NULL, L=1) {
  # L=1 => standard lasso in THETA
  # L=2 => group-lasso in THETA
  if(L==1){
    proximal_solver = lasso_proximal_theta
  } else {
    proximal_solver = group_lasso_proximal_theta
  }
  eta = step_size
  if(is.null(theta_init)){
    theta = matrix(0, nrow=nrow(Y), ncol=nrow(X))
  } else {
    theta=theta_init
  }
  theta_old = theta
  t_old = 1
  cost = numeric(max_iter)
  for(i in 1:max_iter){
    z = theta + ((t_old - 1)/(t_old+2))*(theta - theta_old)
    grad = grad_ll_mult_theta(z, X, Y)
    # Single step, no backtracking (or you could do that if needed)
    theta_new = proximal_solver(z - eta*grad, lambda*eta)
    D_loop = theta_new - z
    log_loss_new = ll_mult_theta(theta_new, X, Y)
    if(L==1){
      cost[i] = log_loss_new + lambda*sum(abs(theta_new))
    } else {
      cost[i] = log_loss_new + lambda*sum(sqrt(rowSums(theta_new^2)))
    }
    if(i>1){
      if(abs(cost[i]-cost[i-1]) < epsilon && sum(D_loop^2)<epsilon && i>min_iter){
        break
      }
    }
    theta_old = theta
    theta = theta_new
    t_old = t_old + 1
  }
  return(list(theta_final=theta, cost=cost[1:i]))
}

accelerated_proximal_gradient_descent_path_fixed_eta <- function(X, Y, lam_seq, step_size, epsilon=1e-8,
                                                                 max_iter=100, min_iter=10, L=1) {
  nlambda = length(lam_seq)
  #  theta has dimension (#categories, #features), i.e. nrow(Y) x nrow(X)
  nfeatures = nrow(Y)
  nresponses= nrow(X)
  theta_tensor = array(0, dim=c(nfeatures, nresponses, nlambda))
  theta_init = NULL
  for(lam_ind in 1:nlambda){
    lambda = lam_seq[lam_ind]
    out_ = accelerated_proximal_gradient_descent_fixed_eta(X, Y, lambda, 
                                                           step_size, epsilon, max_iter, 
                                                           min_iter, theta_init, L)
    theta_tensor[,,lam_ind] = out_$theta_final
    theta_init = out_$theta_final
  }
  return(theta_tensor)
}

beta_validation_theta <- function(theta_tensor, lam_seq, X_validation, Y_validation, plot_FT=FALSE){
  nlambda = dim(theta_tensor)[3]
  cross_entropy = rep(0, nlambda)
  for(i in 1:nlambda){
    cross_entropy[i] = ll_mult_theta(theta_tensor[,,i], X_validation, Y_validation)
  }
  min_index = which.min(cross_entropy)
  if(plot_FT){
    plot(log(lam_seq), cross_entropy, type='l', xlab='log(lambda)', ylab='Neg LogLik')
    points(log(lam_seq[min_index]), cross_entropy[min_index], pch=19, col='red')
  }
  return(min_index)
}



tensor_marginal <- function(J, Y, marginal_int){
  nn = ncol(Y)
  output_mar = matrix(0, nrow=J[marginal_int], ncol=nn)
  for(i in 1:nn){
    y_tensor <- array(Y[,i], dim=J)
    output_mar[,i] = apply(y_tensor, marginal_int, sum)
  }
  return(output_mar)
}
### old
beta_validation_glmnet_theta <- function(glmnet_output, X_validation, Y_validation, plot_TF=FALSE){
  # Evaluate ll_mult_theta 
  lam_seq = glmnet_output$lambda
  cross_entropy=rep(0, length(lam_seq))
  for(i in seq_along(lam_seq)){
    coef_list = coef(glmnet_output, s=lam_seq[i])
    # Combine them to a matrix => rows=categories, columns=features
    # For a simple example with 'family="multinomial"', we might do:
    # If the glmnet output is K classes => we get K coef() calls.
    theta_est = t(do.call(cbind, lapply(coef_list, as.matrix)))
    cross_entropy[i] = ll_mult_theta(theta_est, X_validation, Y_validation)
  }
  min_idx = which.min(cross_entropy)
  if(plot_TF){
    plot(log(lam_seq), cross_entropy, type='l')
    points(log(lam_seq[min_idx]), cross_entropy[min_idx], col='red', pch=19)
  }
  theta.min_list <- coef(glmnet_output, s=lam_seq[min_idx])
  theta.min = t(do.call(cbind, lapply(theta.min_list, as.matrix)))
  return(list(lambda_min_index=min_idx, theta.min=theta.min))
}
### new
beta_validation_glmnet_theta <- function(glmnet_output, X_validation, Y_validation, plot_TF=FALSE){
  # -----------------------------------------
  # 1) Quick checks
  # -----------------------------------------
  stopifnot(is.numeric(glmnet_output$lambda))                    # glmnet must have a numeric lambda sequence
  stopifnot(nrow(X_validation) >= 1, ncol(X_validation) >= 1)    # X_validation must be non-empty
  stopifnot(nrow(Y_validation) >= 1, ncol(Y_validation) >= 1)    # Y_validation must be non-empty
  # Typically we want the same # of samples in X_validation and Y_validation:
  stopifnot(ncol(Y_validation) == ncol(X_validation))
  
  # The sequence of lambdas from glmnet
  lam_seq <- glmnet_output$lambda
  # We'll store a cross-entropy (or negative log-likelihood) for each lambda
  cross_entropy <- numeric(length(lam_seq))
  
  # -----------------------------------------
  # 2) Loop over each lambda
  # -----------------------------------------
  for(i in seq_along(lam_seq)) {
    cat("\n-----------------------------\n")
    cat("DEBUG: i =", i, 
        ", lambda =", lam_seq[i],
        ", length(glmnet_output$beta) might vary if family='multinomial'\n")
    
    # coef_list is typically a list of length = #classes for multinomial
    # or length=1 for binomial.
    coef_list <- coef(glmnet_output, s = lam_seq[i])
    cat("DEBUG: length(coef_list) =", length(coef_list), "\n")
    
    # Convert each item to a matrix, e.g. (p+1) x 1, then cbind them
    mat_list <- lapply(coef_list, as.matrix)
    for(ci in seq_along(mat_list)) {
      cat("   Class", ci, ", dim(mat_list[[ci]]) =", dim(mat_list[[ci]]), "\n")
    }
    
    # After do.call(...), we get a (p+1) x K matrix for K classes, then transpose => K x (p+1)
    theta_est <- t(do.call(cbind, mat_list))
    cat("DEBUG: After cbind & transpose, dim(theta_est) =", dim(theta_est), "\n")
    
    # Show dimension info for X_validation, Y_validation
    cat("DEBUG: dim(X_validation) =", dim(X_validation),
        ", dim(Y_validation) =", dim(Y_validation), "\n")
    
    # -----------------------------------------
    # 3) Evaluate cross-entropy or negative log-likelihood
    # -----------------------------------------
    cat("DEBUG: calling ll_mult_theta...\n")
    cross_entropy[i] <- ll_mult_theta(theta_est, X_validation, Y_validation)
    cat("DEBUG: cross_entropy[", i, "] =", cross_entropy[i], "\n")
  }
  
  # -----------------------------------------
  # 4) Pick best lambda
  # -----------------------------------------
  min_idx <- which.min(cross_entropy)
  cat("\n====================================\n")
  cat("DEBUG: best lambda index =", min_idx, ", => lambda =", lam_seq[min_idx], "\n")
  cat("DEBUG: cross_entropy =", cross_entropy, "\n")
  
  # (Optional) Plot the cross-entropy curve
  if(plot_TF){
    plot(log(lam_seq), cross_entropy, type='l',
         xlab="log(lambda)", ylab="Cross Entropy (neg LL)",
         main="Cross Entropy vs. Lambda")
    points(log(lam_seq[min_idx]), cross_entropy[min_idx], col='red', pch=19)
  }
  
  # -----------------------------------------
  # 5) Extract final theta at best lambda
  # -----------------------------------------
  theta.min_list <- coef(glmnet_output, s = lam_seq[min_idx])
  # again, might be length=#classes or length=1 for binomial
  mat_final <- lapply(theta.min_list, as.matrix)
  theta.min <- t(do.call(cbind, mat_final))
  cat("DEBUG: final best model => dim(theta.min) =", dim(theta.min), "\n")
  
  # return everything useful
  return(list(
    lambda_min_index  = min_idx,
    theta.min         = theta.min,
    lam_seq           = lam_seq,
    cross_entropy     = cross_entropy
  ))
}



tensor_marginal <- function(J, Y, marginal_int) {
  # J is c(2,2,4)
  # Y is prod(J) x #samples
  nn = ncol(Y)
  output_mar = matrix(0, nrow=J[marginal_int], ncol=nn)
  for(i in 1:nn){
    y_vec = Y[,i]  # shape = prod(J)
    y_arr = array(y_vec, dim=J)  # turn into a 3D array of shape (2,2,4)
    # Summation along all dimensions except marginal_int
    output_mar[,i] = apply(y_arr, marginal_int, sum)
  }
  return(output_mar)
}





#-----------------------------------------------
# Now load your real data set
#-----------------------------------------------
dataset = read.csv("/panfs/jay/groups/27/shenx/zhao1118/AD-real-data/dataset.csv")
dfReal  = dataset
dfReal[[1]] = as.factor(dfReal[[1]])
dfReal[[2]] = as.factor(dfReal[[2]])
dfReal[[3]] = as.factor(dfReal[[3]])

make_dummy_for_d <- function(d, p, Model){
  # For weighting, we do:
  ind_end_X_dummy <- 1:p  # default
  # We build J=c(2,2,4)
  J <- c(2,2,4)
  # We do a dummy call:
  out_dummy <- generate_multinomial_final(J, d, p, n_sample=2,
                                          ind_start_X=1:10,  # doesn't matter
                                          ind_end_X=1:10,
                                          r=1/2,
                                          s=length(J),
                                          beta=NULL, R=NULL)
  H_ <- out_dummy$H
  ind_start_H_ <- out_dummy$ind_start_H
  ind_end_H_   <- out_dummy$ind_end_H
  L_sum_ <- ncol(H_)
  
  # We can also do a dummy beta:
  set.seed(7777)
  beta_sim_ <- beta_generate_scheme(J,d,L_sum_,p,
                                    ind_start_H_,ind_end_H_,
                                    unif_range=c(-2,2),
                                    scheme=Model)
  beta_sim_ <- signal_strengthen(beta_sim_, delta=2, sigma=2)
  
  # We'll define a dummy "lambda.min.ratio_seq":
  lambda.min.ratio_seq_ <- c(0.6, 0.6, 0.6, 0.6, 0.6)
  
  return(list(H=H_,
              ind_start_H=ind_start_H_,
              ind_end_H=ind_end_H_,
              L_sum=L_sum_,
              beta_sim=beta_sim_,
              lambda.min.ratio_seq=lambda.min.ratio_seq_))
}




for(ind_loop in seq_len(nrow(params))) {
  p     <- params[ind_loop, "p"]
  Model <- params[ind_loop, "Model"]
  n     <- params[ind_loop, "n"]
  rate  <- params[ind_loop, "rate"]   # NEW: retrieve the rate
  
  # Define tensor dimensions (assumed constant)
  J <- c(2, 2, 4)
  
  # STEP 2: Load the real data (dfReal is defined above)
  for(i in ID_index + (0:(nreps/num_parallel - 1))*num_parallel) {
    
    # Build Y, X
    Y_full <- convert_response(dfReal, J = J)
    X_full <- make_X(dfReal, addIntercept = TRUE, colSample = p, seed = 2024)
    
    n_train      <- n
    n_validation <- 150 - n
    n_test       <- 52
    stopifnot(n_train + n_validation + n_test == 202)
    
    # Optionally shuffle per replication
    set.seed(2024 + i)
    perm <- sample(202)
    X_full_per <- X_full[, perm]
    Y_full_per <- Y_full[, perm]
    
    X_train <- X_full_per[, 1:n_train]
    Y_train <- Y_full_per[, 1:n_train]
    
    X_validation <- X_full_per[, (n_train + 1):(n_train + n_validation)]
    Y_validation <- Y_full_per[, (n_train + 1):(n_train + n_validation)]
    
    X_test <- X_full_per[, (n_train + n_validation + 1):(n_train + n_validation + n_test)]
    Y_test <- Y_full_per[, (n_train + n_validation + 1):(n_train + n_validation + n_test)]
    
    # 13 methods:
    mis_class   = matrix(0, nrow = 1, ncol = 13)
    joint_CE_mat = matrix(0, nrow = 1, ncol = 13)
    nlam        <- 200
    
    # Loop over methods (cases 1 to 13)
    for(caseID in 1:13) {
      # Special case for theta-space group-lasso (caseID 13)
      if(caseID == 13) {
        cat("Case 13 => group-lasso in theta space (old #5)\n")
        lambda_max_G_lasso <- function(X, Y, w = 1) {
          theta0 = matrix(0, nrow = nrow(Y), ncol = nrow(X))
          grad0  = grad_ll_mult_theta(theta0, X, Y)
          rownorms = sqrt(rowSums(grad0^2))
          return(max(rownorms))
        }
        lam_max_13 = lambda_max_G_lasso(X_train, Y_train, w = 1)
        lam_seq_13 = lambda_seq(lam_max_13, lambda.min.ratio = 0.85, nlambda = nlam)
        step_theta_13 = 2 * n / norm(X_train, "2")^2
        theta_tensor_13 = accelerated_proximal_gradient_descent_path_fixed_eta(
          X_train, Y_train,
          lam_seq_13,
          step_size = step_theta_13,
          epsilon = 1e-4,
          max_iter = 1000,
          min_iter = 5,
          L = 2
        )
        min_idx_13 = beta_validation_theta(theta_tensor_13, lam_seq_13, X_validation, Y_validation)
        theta_star_13 = theta_tensor_13[,,min_idx_13]
        mis_class[1, caseID] = calculate_misclassification_rate(
          X_test, Y_test, beta_est = NULL, H = NULL, theta = theta_star_13
        )
        joint_CE_mat[1, caseID] = ll_mult_theta(theta_star_13, X_test, Y_test)
        
      } else {
        # For cases 1 to 12, determine d_val based on caseID:
        if(caseID >= 1 && caseID <= 4) {
          d_val = 3
        } else if(caseID >= 5 && caseID <= 8) {
          d_val = 2
        } else if(caseID >= 9 && caseID <= 12) {
          d_val = 1
        }
        
        # Retrieve design parameters via a dummy call
        dummyObj <- make_dummy_for_d(d_val, p, Model)
        H_           <- dummyObj$H
        ind_start_H_ <- dummyObj$ind_start_H
        ind_end_H_   <- dummyObj$ind_end_H
        lambda.min.ratio_seq_ <- dummyObj$lambda.min.ratio_seq
        
        # Set predictor grouping depending on method type.
        if(caseID %in% c(1,2,5,6,9,10)) {  # L- methods
          ind_end_X_ = 1:p
          ind_start_X_ = c(1, 1 + ind_end_X_[-length(ind_end_X_)])
        } else if(caseID %in% c(3,4,7,8,11,12)) {  # G- methods
          ind_end_X_ = c(1, p)
          ind_start_X_ = c(1, 2)
        }
        
        # Choose method branch based on caseID modulo 4:
        if(caseID %% 4 == 1) {
          # => L-Mult branch
          cat(sprintf("Case %d => d=%d => L-Mult...\n", caseID, d_val))
          lambda.min.ratio = lambda.min.ratio_seq_[1]
          w_mat <- custom_weight_generate(J, d_val, ind_end_H_, ind_end_X_, rate = rate)
          lambda_mult_max = lambda_max_mult(H_, X_train, Y_train,
                                            ind_start_H_, ind_end_H_,
                                            ind_start_X_, ind_end_X_, w = w_mat)
          lam_seq_mult = lambda_seq(lambda_mult_max, lambda.min.ratio, nlambda = nlam)
          eta_mult_ = prod(c(2,2,4)) * n / norm(X_train, "2")^2
          bpath = backtracking_PGD_path(grad_ll_mult, ll_mult,
                                        H_, X_train, Y_train,
                                        step_size = eta_mult_,
                                        ind_start_H_, ind_end_H_,
                                        ind_start_X_, ind_end_X_,
                                        lam_seq_mult,
                                        shrinking_gamma = 0.618,
                                        epsilon = 1e-4,
                                        max_iter = 1000,
                                        w = w_mat)
          bTensor = bpath$beta_tensor
          val_out = beta_validation(bTensor, lam_seq_mult, ll_mult,
                                    H_, X_validation, Y_validation,
                                    n_train = n, plot_TF = FALSE)
          best_ind = val_out$lambda_min_index
          mis_class[1, caseID] = calculate_misclassification_rate(
            X_test, Y_test, bTensor[,,best_ind], H_
          )
          joint_CE_mat[1, caseID] = ll_mult(bTensor[,,best_ind], H_, X_test, Y_test)
          
        } else if(caseID %% 4 == 2) {
          # => L-Pois branch
          cat(sprintf("Case %d => d=%d => L-Pois...\n", caseID, d_val))
          lambda.min.ratio = lambda.min.ratio_seq_[2]
          w_mat <- custom_weight_generate(J, d_val, ind_end_H_, ind_end_X_, rate = rate)
          lam_pois_max = lambda_max_pois(H_, X_train, Y_train,
                                         ind_start_H_, ind_end_H_,
                                         ind_start_X_, ind_end_X_, w = w_mat)
          lam_seq_pois = lambda_seq(lam_pois_max, lambda.min.ratio, nlam)
          eta_pois_ = n / norm(X_train, "2")^2
          bpath = backtracking_PGD_path(grad_ll_pois, ll_pois,
                                        H_, X_train, Y_train,
                                        step_size = eta_pois_,
                                        ind_start_H_, ind_end_H_,
                                        ind_start_X_, ind_end_X_,
                                        lam_seq_pois,
                                        shrinking_gamma = 0.618,
                                        epsilon = 1e-4,
                                        max_iter = 1000,
                                        w = w_mat)
          bTensor = bpath$beta_tensor
          val_out = beta_validation(bTensor, lam_seq_pois, ll_mult,
                                    H_, X_validation, Y_validation,
                                    n_train = n, plot_TF = FALSE)
          best_ind = val_out$lambda_min_index
          mis_class[1, caseID] = calculate_misclassification_rate(
            X_test, Y_test, bTensor[,,best_ind], H_
          )
          joint_CE_mat[1, caseID] = ll_pois(bTensor[,,best_ind], H_, X_test, Y_test)
          
        } else if(caseID %% 4 == 3) {
          # => G-Mult branch
          cat(sprintf("Case %d => d=%d => G-Mult...\n", caseID, d_val))
          lambda.min.ratio = lambda.min.ratio_seq_[3]
          w_mat <- custom_weight_generate(J, d_val, ind_end_H_, ind_end_X_, rate = rate)
          lam_mult_max = lambda_max_mult(H_, X_train, Y_train,
                                         ind_start_H_, ind_end_H_,
                                         ind_start_X_, ind_end_X_, w = w_mat)
          lam_seq_mult = lambda_seq(lam_mult_max, lambda.min.ratio, nlam)
          eta_mult_ = prod(c(2,2,4)) * n / norm(X_train, "2")^2
          bpath = backtracking_PGD_path(grad_ll_mult, ll_mult,
                                        H_, X_train, Y_train,
                                        eta_mult_,
                                        ind_start_H_, ind_end_H_,
                                        ind_start_X_, ind_end_X_,
                                        lam_seq_mult,
                                        shrinking_gamma = 0.618,
                                        epsilon = 1e-4,
                                        max_iter = 1000,
                                        w = w_mat)
          bTensor = bpath$beta_tensor
          val_out = beta_validation(bTensor, lam_seq_mult, ll_mult,
                                    H_, X_validation, Y_validation,
                                    n_train = n, plot_TF = FALSE)
          best_ind = val_out$lambda_min_index
          mis_class[1, caseID] = calculate_misclassification_rate(
            X_test, Y_test, bTensor[,,best_ind], H_
          )
          joint_CE_mat[1, caseID] = ll_mult(bTensor[,,best_ind], H_, X_test, Y_test)
          
        } else if(caseID %% 4 == 0) {
          # => G-Pois branch
          cat(sprintf("Case %d => d=%d => G-Pois...\n", caseID, d_val))
          lambda.min.ratio = lambda.min.ratio_seq_[4]
          w_mat <- custom_weight_generate(J, d_val, ind_end_H_, ind_end_X_, rate = rate)
          lam_pois_max = lambda_max_pois(H_, X_train, Y_train,
                                         ind_start_H_, ind_end_H_,
                                         ind_start_X_, ind_end_X_, w = w_mat)
          lam_seq_pois = lambda_seq(lam_pois_max, lambda.min.ratio, nlam)
          eta_pois_ = n / norm(X_train, "2")^2
          bpath = backtracking_PGD_path(grad_ll_pois, ll_pois,
                                        H_, X_train, Y_train,
                                        step_size = eta_pois_,
                                        ind_start_H_, ind_end_H_,
                                        ind_start_X_, ind_end_X_,
                                        lam_seq_pois,
                                        shrinking_gamma = 0.618,
                                        epsilon = 1e-4,
                                        max_iter = 1000,
                                        w = w_mat)
          bTensor = bpath$beta_tensor
          val_out = beta_validation(bTensor, lam_seq_pois, ll_mult,
                                    H_, X_validation, Y_validation,
                                    n_train = n, plot_TF = FALSE)
          best_ind = val_out$lambda_min_index
          mis_class[1, caseID] = calculate_misclassification_rate(
            X_test, Y_test, bTensor[,,best_ind], H_
          )
          joint_CE_mat[1, caseID] = ll_pois(bTensor[,,best_ind], H_, X_test, Y_test)
        } # end if for modulo branch
      } # end if else for caseID==13 vs. others
    } # end for(caseID in 1:13)
    
    # Save results
    savename <- paste0("/home/shenx/zhao1118/AD-real-data/data_CE/Model", Model, 
                       "_p", p, "_n", n,"_rate", rate,"_repID", i, ".RDS")
    Results <- list(
      Misclassification_rate = mis_class,
      Joint_CE = joint_CE_mat,
      Model = Model,
      p = p,
      n = n,
      rate = rate,
      rep_ID = i
    )
    saveRDS(Results, file = savename)
  } # end replication loop
} # end for(ind_loop in seq_len(nrow(params)))



