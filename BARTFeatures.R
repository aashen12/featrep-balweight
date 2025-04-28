# Fit BART features

library(Matrix)
library(stochtree)


extract_bart_features <- function(X_train, y_train, data, covs, n_components = 100, num_trees = 100, seed = 122345, verbose = FALSE) {
  
  X_test <- data %>% dplyr::select(all_of(covs))
  sigma_leaf <- 1/num_trees
  mean_forest_params <- list(num_trees=num_trees, sigma2_leaf_init=sigma_leaf)
  bart_model <- bart(X_train=X_train, y_train=y_train, mean_forest_params = mean_forest_params,
                     num_burnin=200, num_mcmc=1000, general_params=list(random_seed=seed))
  
  nkeep <- n_components
  nsamp = bart_model$model_params$num_samples
  wts0 = matrix(0, nrow=nrow(X_train), ncol=nkeep)
  wts1 <- matrix(0, nrow=nrow(X_test), ncol=nkeep)
  ix = round(seq(1, bart_model$model_params$num_samples, length.out=nkeep))
  i = 1
  for(ii in ix) {
    leaf_mat_train <- computeForestLeafIndices(bart_model, X_train, forest_type = "mean", 
                                               forest_inds = ii - 1)
    leaf_mat_test <- computeForestLeafIndices(bart_model, X_test, forest_type = "mean", 
                                              forest_inds = ii - 1)
    
    var_leaf = bart_model$sigma2_leaf_samples[ii]
    sigma2  = bart_model$sigma2_global_samples[ii]
    W_train <- sparseMatrix(i=rep(1:nrow(X_train),num_trees), j=leaf_mat_train + 1, x=1)*sqrt(var_leaf)
    W_test  <- sparseMatrix(i=rep(1:nrow(X_test),num_trees), j=leaf_mat_test + 1, x=1, dims=c(nrow(X_test), ncol(W_train)))*sqrt(var_leaf)
    
    #  for(i in round(seq(1, nrow(X_test), length.out=9))) {
    #  wts[,i] = as.matrix(solve(W_train%*%t(W_train) + diag(sigma2, nrow(W_train)), W_train%*%colMeans(W_test)))
    wts0[,i] = as.matrix(W_train%*%as.matrix(solve(t(W_train)%*%W_train + diag(sigma2, ncol(W_train)), colMeans(W_test))))
    wts1[,i] = as.matrix(W_test%*%as.matrix(solve(t(W_train)%*%W_train + diag(sigma2, ncol(W_train)), colMeans(W_test))))
    i = i+1
    # K = K + W_train%*%t(W_train)
    
  }
  wts0 <- data.frame(wts0)
  wts1 <- data.frame(wts1)

  bart_features <- wts1
  names(bart_features) <- paste0("BART_", seq_len(ncol(bart_features)))
  data_bart <- cbind(data, bart_features)
  out_list <-  list(
    features = bart_features,
    data_bart = data_bart
  )
}






