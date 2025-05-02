eval_data <- function(dat, pilot.dat, treat.true = 5, verbose = FALSE) {
  
  data.c <- dat %>% filter(Z==0)
  data.t <- dat %>% filter(Z==1)
  
  # set.seed(1031)
  # train <- sample(nrow(data.c), round(.1*nrow(data.c)))
  # length(train)
  # data.c.train <- data.c[train, ]
  # data.c.test <- data.c[-train, ]
  data.c.train <- pilot.dat
  data <- rbind(data.t, data.c)
  # gets rid of the training sample used to select interactions & nonlinear covariates
  
  covs <- c("X1","X2", "X3", "X4", "X5")
  covs_gam_mod <- covs
  gam_fit <- fit_gam(data.c.train, covs_gam = covs_gam_mod, cutoff = 3)
  nonlinear_covs <- gam_fit$nonlinear_covs
  if (verbose) {
    print("gam fitted")
  }
  ############################################################################
  
  
  ### Generate splines ###
  # nc_spline <- c(3, 5, 10, 15)
  # spline.scenarios <- expand.grid(fr = "spline", ncomp = nc_spline)
  # 
  # spline_dfs <- lapply(nc_spline, function(i) {
  #   spline_obj <- generate_splines(data, nonlinear_covs = nonlinear_covs, df = i)
  #   data_spline <- spline_obj$data_spline
  #   data_spline
  # })
  # names(spline_dfs) <- paste0("spline_", nc_spline)
  # if (verbose) {
  #   print("splines generated")
  # }
  ############################################################################
  
  
  #### Generate Random Forest features #####
  nc_rf <- c(10, 20, 50, 100)
  rf.scenarios <- expand.grid(fr = c("rf_only", "rf_plus"), ncomp = nc_rf)
  ## First, run random forest
  
  form <- reformulate(covs, response = "Y")
  rfmod <- randomForest(form, data = data.c.train, ntree = 100)
  
  # Next, extract feat reps
  all_rf <- lapply(nc_rf, function(nc) {
    rf_obj <- extract_rf_features(rf_model = rfmod, data = data, n_components = nc, verbose = TRUE)
    data_rf_only <- rf_obj$features %>% 
      dplyr::mutate(Y = data$Y, treat = data$Z)
    var_exp <- rf_obj$explained_variance
    data_rf_plus <- rf_obj$data_rf
    ol <- list(
      data_rf_only = data_rf_only,
      data_rf_plus = data_rf_plus
    )
    ol
  })
  names(all_rf) <- paste0("rf_", nc_rf)
  ############################################################################
  
  ############################# Extract Gaussian Kernel Features #################################
  kbal.scenarios <- expand.grid(fr = c("kbal_only", "kbal_plus"), ncomp = nc_rf)
  X_unscaled <- model.matrix(reformulate(c(covs, "-1")), data)
  # scale to have variance 1 but not mean 0
  X <- scale(X_unscaled)
  kbal_objs <- lapply(1:length(nc_rf), function(i) {
    print(i)
    nc <- nc_rf[i]
    kbal_obj <- kbal::kbal(X, treatment = data$Z, numdims = nc)
    data_kbal_only <-  data.frame(
      kbal_obj$svdK$u[, 1:nc] %*% diag(sqrt(kbal_obj$svdK$d[1:nc]))
    )
    names(data_kbal_only) <- paste0("PC", 1:nc)
    var_exp <- kbal_obj$explained_variance
    data_kbal_plus <- cbind(data, data_kbal_only)
    ol <- list(
      data_kbal_only = data_kbal_only %>% dplyr::mutate(treat = data$Z, Y = data$Y),
      data_kbal_plus = data_kbal_plus
    )
    ol
  }
  )
  names(kbal_objs) <- paste0("kbal_", nc_rf)
  ###########################################################################
  
  raw.scenarios <- expand.grid(fr = "raw", ncomp = 0)
  ### Run the actual simulation ###
  scenarios <- rbind(raw.scenarios, rf.scenarios, kbal.scenarios)
  
  out <- lapply(1:nrow(scenarios), function(i){
    fr <- scenarios[i,1]
    nc <- scenarios[i,2]
    feat_rep <- paste(fr, nc, sep = "_")
    if (fr == "raw") {
      dataset <- data
    } else if (fr == "rf_plus") {
      dataset <- all_rf[[paste0("rf_", nc)]]$data_rf_plus
    } else if (fr == "rf_only") {
      dataset <- all_rf[[paste0("rf_", nc)]]$data_rf_only
    } else if (fr == "spline") {
      dataset <- spline_dfs[[paste0("spline_", nc)]]
    } else if (fr == "kbal_only") {
      dataset <- kbal_objs[[paste0("kbal_", nc)]]$data_kbal_only
    } else if (fr == "kbal_plus") {
      dataset <- kbal_objs[[paste0("kbal_", nc)]]$data_kbal_plus
    }
     
    # compute ATT, balance metrics, error metrics, and so forth for each estimator 
    bw_l2 <- balancingWeights(data = dataset, true_att = treat.true, feat_rep = feat_rep, type = "l2", verbose = FALSE)
    bw_inf <- balancingWeights(data = dataset, true_att = treat.true, feat_rep = feat_rep, type = "inf", verbose = FALSE)
    ipw <- logisticIPW(data = dataset, true_att = treat.true, feat_rep = feat_rep, verbose = FALSE)
    aug.l2 <- augmentedBalWeights(data = dataset, true_att = treat.true, feat_rep = feat_rep, out.mod = "ridge", ps.mod = "l2", verbose = FALSE)
    aug.inf <- augmentedBalWeights(data = dataset, true_att = treat.true, feat_rep = feat_rep, out.mod = "lasso", ps.mod = "inf", verbose = FALSE)
    aug.vanilla <- augmentedBalWeights(data = dataset, true_att = treat.true, feat_rep = feat_rep, out.mod = "ols", ps.mod = "ipw", verbose = FALSE)
    list(bw_l2 = bw_l2, bw_inf = bw_inf, ipw = ipw, 
         aug.l2 = aug.l2, aug.inf = aug.inf, aug.vanilla = aug.vanilla)#, var_exp = var_exp)
  })
  out
}