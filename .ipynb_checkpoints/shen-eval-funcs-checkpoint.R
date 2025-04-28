eval_data <- function(dat, treat.true, verbose = FALSE) {
  
  data <- dat
  ################################################################################################
  
  data.c <- data %>% dplyr::filter(Z==0)
  data.t <- data %>% dplyr::filter(Z==1)
  
  set.seed(1031)
  train <- sample(nrow(data.c), round(0.05 * nrow(data.c)))
  length(train)
  data.c.train <- data.c[train, ]
  data.c.test <- data.c[-train, ]
  
  data <- rbind(data.t, data.c.test)
  # gets rid of the training sample used to select interactions & nonlinear covariates
  
  ##### FIT GAM (should probably throw into a different function) #####
  covs <- names(data.c.train)[!names(data.c.train) %in% c("Y", "Z")]
  covs_gam_mod <- covs#[grepl("^x_[0-9]+$", covs)]
  gam_fit <- fit_gam(data.c.train, covs_gam = covs_gam_mod, cutoff = 3)
  nonlinear_covs <- gam_fit$nonlinear_covs
  if (verbose) {
    print("gam fitted")
  }
  ############################################################################
  
  
  ### Generate splines ###
  spline_comps <- c(3, 6, 10, 15)
  spline.scenarios <- expand.grid(feat_rep = c("spline"), ncomp = spline_comps)
  
  spline_dfs <- lapply(spline_comps, function(comp) {
    spline_obj <- generate_splines(data, nonlinear_covs = nonlinear_covs, df = comp)
    data_spline <- spline_obj$data_spline
    data_spline
  })
  names(spline_dfs) <- paste0("spline_", spline_comps)
  
  if (verbose) {
    print("splines generated")
  }
  ############################################################################
  #### Generate Random Forest features #####
  ## First, run random forest
  form <- reformulate(covs, response = "Y")
  rfmod <- randomForest(form, data = data.c.train, ntree = 100)
  
  # Next, extract feat reps 
  # specify number of components
  rf_comps <- c(5, 10, 25, 50, 100, 200)
  rf.scenarios <- expand.grid(feat_rep = c("rf_only", "rf_plus"), ncomp = rf_comps)
  
  rf_obj <- extract_rf_features(rfmod, data, n_components = 200)
  data_rf_only <- rf_obj$features %>% 
    dplyr::mutate(Y = data$Y, Z = data$Z)
  data_rf_plus <- rf_obj$data_rf

  ############################################################################
  
  # do separate expand.grids() for spline, raw (0), bart+rf
  
  ### Run the actual simulation ###
  data_raw <- data
  raw.scenarios <- expand.grid(feat_rep = "raw", ncomp = 0)
  
  scenarios <- rbind(raw.scenarios, rf.scenarios, spline.scenarios)
  
  out <- map(
    1:nrow(scenarios),
    function(i) {
      fr <- scenarios[i, 1] 
      nc <- scenarios[i, 2]
      scenario <- paste(fr, nc, sep = "_")
      print(paste("Running scenario", i, "with ncomp", nc, "and fr", fr))
      if (fr != "regression") {
        if (fr == "raw") {
          input_data <- data_raw
          input_covs <- covs
        } else if (fr == "rf_plus") {
          input_covs <- c(covs, paste0("PC", 1:nc))
          input_data <- data_rf_plus %>% dplyr::select(Y, Z, all_of(input_covs))
        } else if (fr == "rf_only") {
          input_covs <- c(paste0("PC", 1:nc))
          input_data <- data_rf_only %>% dplyr::select(Y, Z, all_of(input_covs))
        } else if(fr == "bart_plus") {
          input_data <- all_bart[[paste0("bart_", nc)]]$data_bart_plus
          input_covs <- c(covs, paste0("BART_", 1:nc))
        } else if (fr == "bart_only") {
          input_data <- all_bart[[paste0("bart_", nc)]]$data_bart_only
          input_covs <- c(paste0("BART_", 1:nc))
        } else if (fr == "spline") {
          input_data <- spline_dfs[[paste0("spline_", nc)]]
        }
        if (verbose) print("STARTING L2 BW MODELS")
        bw_l2 <- balancingWeights(data = input_data, true_att = treat.true, feat_rep = scenario, type = "l2", verbose = verbose)
        bw_inf <- balancingWeights(data = input_data, true_att = treat.true, feat_rep = scenario, type = "inf", verbose = verbose)
        
        if (verbose) print("STARTING VANILLA IPW MODELS")
        ipw <- logisticIPW(data = input_data, true_att = treat.true, feat_rep = scenario, verbose = verbose)
        
        if (verbose) print("STARTING OLS MODELS")
        ols <- outcomeRegression(data = input_data, true_att = treat.true, feat_rep = scenario, type = "ols", verbose = verbose)
        # ridge <- outcomeRegression(data = dataset, true_att = treat.true, feat_rep = scenario, type = "ridge", verbose = verbose)
        # lasso <- outcomeRegression(data = dataset, true_att = treat.true, feat_rep = scenario, type = "lasso", verbose = verbose)
        
        if (verbose) print("STARTING AUGMENTED MODELS")
        aug.l2 <- augmentedBalWeights(data = input_data, true_att = treat.true, feat_rep = scenario, out.mod = "ridge", ps.mod = "l2", verbose = verbose)
        aug.inf <- augmentedBalWeights(data = input_data, true_att = treat.true, feat_rep = scenario, out.mod = "lasso", ps.mod = "inf", verbose = verbose)
        
        list(bw_l2 = bw_l2, bw_inf = bw_inf, ipw = ipw, ols = ols, aug.l2 = aug.l2, aug.inf = aug.inf)
      }
    }
  )
  #scenarios_paste <- paste(scenarios[, 1], scenarios[,2], sep = ".")
  #names(out) <- scenarios_paste
  out
}  
  
  
#   out <- lapply(1:length(feat_reps), function(i){
#     fr <- feat_reps[i]
#     dataset <- switch(fr, "spline3" = data_spline_3, "spline6" = data_spline_6, 
#                       "raw" = data, "rf10" = data_rf_10, "rf20" = data_rf_20)
#     
#     # compute ATT, balance metrics, error metrics, and so forth for each estimator
#     if (verbose) print("STARTING L2 BW MODELS")
#     bw_l2 <- balancingWeights(data = dataset, true_att = treat.true, feat_rep = fr, type = "l2", verbose = verbose)
#     bw_inf <- balancingWeights(data = dataset, true_att = treat.true, feat_rep = fr, type = "inf", verbose = verbose)
#     
#     if (verbose) print("STARTING VANILLA IPW MODELS")
#     ipw <- logisticIPW(data = dataset, true_att = treat.true, feat_rep = fr, verbose = verbose)
#     
#     if (verbose) print("STARTING OLS MODELS")
#     ols <- outcomeRegression(data = dataset, true_att = treat.true, feat_rep = fr, type = "ols", verbose = verbose)
#     # ridge <- outcomeRegression(data = dataset, true_att = treat.true, feat_rep = fr, type = "ridge", verbose = verbose)
#     # lasso <- outcomeRegression(data = dataset, true_att = treat.true, feat_rep = fr, type = "lasso", verbose = verbose)
#     
#     if (verbose) print("STARTING AUGMENTED MODELS")
#     aug.l2 <- augmentedBalWeights(data = dataset, true_att = treat.true, feat_rep = fr, out.mod = "ridge", ps.mod = "l2", verbose = verbose)
#     aug.inf <- augmentedBalWeights(data = dataset, true_att = treat.true, feat_rep = fr, out.mod = "lasso", ps.mod = "inf", verbose = verbose)
#     
#     list(bw_l2 = bw_l2, bw_inf = bw_inf, ipw = ipw, ols = ols, aug.l2 = aug.l2, aug.inf = aug.inf)
#   })
#   names(out) <- feat_reps
#   out
# }