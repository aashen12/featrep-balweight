## this script contains functions for random forest feature representation
library(tidyverse)
library(randomForest)
library(Matrix)
library(irlba)  # For efficient truncated SVD


extract_rf_features <- function(rf_model, data, n_components=10, verbose = FALSE) {
  
  # This function extract 1,0's for RF feature representation
  
  # Get leaf node IDs for each tree
  n_components <- min(n_components, nrow(data))
  n_trees <- rf_model[["ntree"]]
  #print(n_trees)
  n_samples <- nrow(data)
  
  
  # Initialize sparse matrix list
  preds <- predict(rf_model, data, nodes=TRUE)
  leaf_ids <- attributes(preds)$nodes 
  encoding_list <- lapply(1:n_trees, function(i){
    #print(i)
    tree_leaf_ids <- leaf_ids[, i]
    
    # Create sparse one-hot encoding
    unique_leaves <- unique(tree_leaf_ids)
    one_hot <- sparseMatrix(
      i = 1:n_samples,
      j = match(tree_leaf_ids, unique_leaves),
      x = 1,
      dims = c(n_samples, length(unique_leaves))
    )
    one_hot
  })
  
  # Combine all encodings
  leaf_encoding <- do.call(cbind, encoding_list)
  
  # Perform truncated SVD
  if (verbose) print(paste("Finished creating matrix. Running SVD."))
  svd_result <- irlba(leaf_encoding, nv=n_components, maxit = 2000, verbose = F)
  
  # Scale by square root of singular values
  features <- svd_result$u %*% diag(sqrt(svd_result$d)) %>% data.frame()
  names(features) <- paste0("PC", 1:n_components)
  data_rf <- cbind(data, features)
  
  return(list(
    #leaf_encoding = leaf_encoding,
    data_rf = data_rf,
    features = features,
    explained_variance = svd_result$d^2 / sum(svd_result$d^2)
  ))
}







