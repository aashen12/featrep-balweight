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
  
  leaf_encoding <- scale(leaf_encoding, center = TRUE, scale = TRUE)
  
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


randomForestRegression <- function(data, true_att, feat_rep, ncomp = 5, verbose = FALSE) {
  #print(paste("True ATT:", round(true_att, 3)))
  data <- data %>% 
    dplyr::rename(treat = Z) #%>% dplyr::select(-Y0, -Y1)
  
  covs <- names(data)[!names(data) %in% c("Y", "treat")]
  covs <- c(covs, "-1")
  X <- scale(model.matrix(reformulate(covs), data))
  trt <- data$treat
  n <- nrow(data)
  
  data0 <- data %>% dplyr::filter(treat == 0)
  data1 <- data %>% dplyr::filter(treat == 1)
  
  form <- reformulate(covs, response = "Y")
  rfmod <- randomForest(form, data = data0, ntree = 100)
  rfpreds_1 <- predict(rfmod, newdata = data1, type = "response")
  
  # Next, extract feat reps
  rf_obj <- extract_rf_features(rfmod, data1, n_components = ncomp)
  rf_feats <- rf_obj$features # U matrix
  data_rf_feats <- data.frame(Y = data1$Y, rf_feats)
  #return(data_rf_feats)
  lmmod <- lm(Y ~ ., data = data_rf_feats)
  lm_rf_preds <- lmmod$fitted.values #predict(lmmod, newdata = data1, type = "response")
  
  att.raw.rf <- mean(data1$Y) - mean(rfpreds_1)
  att.ols.rf <- mean(data1$Y) - mean(lm_rf_preds)
  
  bias.raw.rf <- true_att - att.raw.rf
  bias.ols.rf <- true_att - att.ols.rf
  
  df_raw_rf <- data.frame(est = "raw.rf", "dgp" = dgp, "feat_rep" = feat_rep, 
                          "bias" = bias.raw.rf, "cvg" = NA, pbr = NA, ess = NA,
                          "true.att" = true_att, "est.att" = att.raw.rf)
  
  df_ols_rf <- data.frame(est = "ols.rf", "dgp" = dgp, "feat_rep" = feat_rep, 
                          "bias" = bias.ols.rf, "cvg" = NA, pbr = NA, ess = NA,
                          "true.att" = true_att, "est.att" = att.ols.rf)
  data.frame(bind_rows(df_raw_rf, df_ols_rf))
}



