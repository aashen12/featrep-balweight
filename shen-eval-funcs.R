eval_data <- function(dat, treat.true = 5, verbose = FALSE) {
  
  data.c <- dat %>% filter(Z==0)
  data.t <- dat %>% filter(Z==1)
  
  set.seed(1031)
  train <- sample(nrow(data.c), round(.1*nrow(data.c)))
  length(train)
  data.c.train <- data.c[train, ]
  data.c.test <- data.c[-train, ]
  
  data <- rbind(data.t, data.c.test)
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
  nc_spline <- c(3, 5, 10, 15)
  spline.scenarios <- expand.grid(fr = "spline", ncomp = nc_spline)
  
  spline_dfs <- lapply(nc_spline, function(i) {
    spline_obj <- generate_splines(data, nonlinear_covs = nonlinear_covs, df = i)
    data_spline <- spline_obj$data_spline
    data_spline
  })
  names(spline_dfs) <- paste0("spline_", nc_spline)
  if (verbose) {
    print("splines generated")
  }
  ############################################################################
  
  
  #### Generate Random Forest features #####
  nc_rf <- c(10, 20, 50, 100, 200)
  rf.scenarios <- expand.grid(fr = c("rf_only", "rf_plus"), ncomp = nc_rf)
  ## First, run random forest
  
  form <- reformulate(covs, response = "Y")
  rfmod <- randomForest(form, data = data.c.train, ntree = 100)
  
  # Next, extract feat reps
  rf_obj <- extract_rf_features(rfmod, data, n_components = 200)
  data_rf_plus <- rf_obj$data_rf
  data_rf_only <- rf_obj$features %>% dplyr::mutate(Z = data$Z, Y = data$Y)
  ############################################################################
  
  raw.scenarios <- expand.grid(fr = "raw", ncomp = 0)
  ### Run the actual simulation ###
  scenarios <- rbind(raw.scenarios, rf.scenarios, spline.scenarios)
  
  out <- lapply(1:nrow(scenarios), function(i){
    fr <- scenarios[i,1]
    nc <- scenarios[i,2]
    feat_rep <- paste(fr, nc, sep = "_")
    if (fr == "raw") {
      dataset <- data
    } else if (fr == "rf_only") {
      dataset <- data_rf_only %>% 
        dplyr::select(Y, Z, all_of(paste0("PC", 1:nc)))
    } else if (fr == "rf_plus") {
      dataset <- data_rf_plus %>% 
        dplyr::select(Y, Z, all_of(paste0("PC", 1:nc)), all_of(covs))
    } else if (fr == "spline") {
      dataset <- spline_dfs[[paste0("spline_", nc)]]
    }
     
    # compute ATT, balance metrics, error metrics, and so forth for each estimator 
    bw_l2 <- balancingWeights(data = dataset, true_att = treat.true, feat_rep = feat_rep, type = "l2", verbose = FALSE)
    bw_inf <- balancingWeights(data = dataset, true_att = treat.true, feat_rep = feat_rep, type = "inf", verbose = FALSE)
    ipw <- logisticIPW(data = dataset, true_att = treat.true, feat_rep = feat_rep, verbose = FALSE)
    aug.l2 <- augmentedBalWeights(data = dataset, true_att = treat.true, feat_rep = feat_rep, out.mod = "ridge", ps.mod = "l2", verbose = FALSE)
    list(bw_l2 = bw_l2, bw_inf = bw_inf, ipw = ipw, aug.l2 = aug.l2)
  })
  out
}