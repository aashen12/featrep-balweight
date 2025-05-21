# hyperparameter estimation
estimate_regularization <- function(data, covs) {
  data.c <- data %>% dplyr::filter(treat==0)
  lambda.reg <- lm(reformulate(covs, response="scale(Y)"), data=data.c)
  l <- var(lambda.reg$resid)
  l
}

msm.out <- function(obj){
  SE <- coef(summary(obj))[2,2] 
  beta <- coef(obj)[2]
  lcl <- (beta - abs(qnorm(.025))*SE)
  ucl <- (beta + abs(qnorm(.025))*SE)
  return(c("beta" = beta, "lcl" = lcl, "ucl" = ucl, "SE" = SE))
}


create_sbw_constraints <- function(X, target, tol) {
  
  # balance constraint
  A1 <- t(X) / nrow(X)
  l1 <- target - tol
  u1 <- target + tol
  
  # sum to n constraint
  A2 <- rep(1, nrow(X))
  l2 <- nrow(X)
  u2 <- nrow(X)
  
  # positivity constraint
  A3 <- Matrix::Diagonal(nrow(X))
  l3 <- numeric(nrow(X))
  u3 <- rep(Inf, nrow(X))
  
  return(list(A = rbind(A1, A2, A3),
              l = c(l1, l2, l3),
              u = c(u1, u2, u3)))
}

#' SBW implmenetation with different solver
#' @param X Matrix of covariates (standardized!)
#' @param  Z Vector of treatment assignment
#' @param tol Tolerance for imbalance
#' @return List with weights for all units (treated weights are zero) and imbalance vector
sbw_osqp <- function(X, Z, tol, ...) {
  
  # compute target
  target <- colMeans(X[Z == 1, ])
  
  # P matrix
  n0 <- sum(Z == 0)
  P <- Matrix::sparseMatrix(1:n0, 1:n0, x = rep(1, n0))
  
  consts <- create_sbw_constraints(X[Z == 0,, drop = F], target, tol)
  
  pars <- do.call(osqp::osqpSettings, list(...))
  
  sol <- osqp::solve_osqp(P = P, A = consts$A, l = consts$l, u = consts$u)
  
  wts <- numeric(nrow(X))
  wts[Z == 0] <- sol$x
  
  imbal <- c(t(sol$x) %*% X[Z == 0, ] / sum(Z == 0)) - target
  
  return(list(wts = wts, imbalance = imbal))
  
}



balancingWeights <- function(data, true_att, feat_rep, type = "l2", simplex = TRUE, verbose = FALSE) {
  # l2 or inf
  #print(paste("True ATT:", round(true_att, 3)))
  if ("Z" %in% colnames(data)) {
    data <- data %>% 
      dplyr::rename(treat = Z) #%>% dplyr::select(-Y0, -Y1)
  }
  
  covs <- names(data)[!names(data) %in% c("Y", "treat")]
  print(covs)
  basis <- c(covs, "-1")
  
  contains_rf_only <- grepl("rf_only", feat_rep)
  if (contains_rf_only) {
    X <- model.matrix(reformulate(basis), data)
    # X <- scale(model.matrix(reformulate(basis), data))
  } else {
    X <- scale(model.matrix(reformulate(basis), data))
  }
  
  trt <- data$treat
  n <- nrow(data)
  
  if (type == "l2") {
    # balancing weights
    l <- estimate_regularization(data, covs = covs)
    if (simplex) {
      bal_weights <- multilevel_qp(X, trt, Z = rep(1,n), lambda = l, verbose= FALSE,
                                   exact_global = FALSE, scale_sample_size = FALSE)
    } else {
      bal_weights <- multilevel_qp(X, trt, Z = rep(1,n), lambda = l, verbose= FALSE,
                                   exact_global = FALSE, scale_sample_size = FALSE,
                                   lowlim = -Inf, uplim = Inf)
    }
    data$wts <- pmax(bal_weights$weights, 0)
  } else if (type == "inf") {
    bal_weights <- sbw_osqp(X, Z = trt, tol = 1.75)
    data$wts <- pmax(bal_weights$wts, 0)
    #print(bal_weights$imbal)
  }
  data$wts[data$treat == 1] <- 1
  if (verbose) {
    print("balancing weights estimated")
  }
  
  bal.data <- as.data.frame(X)
  bal.names <- names(bal.data)
  n_covs <- length(bal.names)
  bal.data$treat <- data$treat
  bal.data$wts <- data$wts
  
  # ESS
  data0 <- data %>% dplyr::filter(treat==0)
  n.0 <- (sum(data0$wts)^2) / (sum(data0$wts^2))
  data1 <- data %>% dplyr::filter(treat==1)
  n.1 <- (sum(data1$wts)^2) / (sum(data1$wts^2))
  
  if (verbose) {
    print("ESS calculated")
  }
  
  ## Balance diagnostics
  
  data.var <- bal.data %>% dplyr::group_by(treat) %>% 
    dplyr::summarize(across(all_of(bal.names), ~var(.x))) %>% as.data.frame()
  
  c.var <- as.numeric(data.var[1,])
  t.var <- as.numeric(data.var[2,])
  c.var <- c.var[-1]
  t.var <- t.var[-1]
  pooled.var <- sqrt((t.var + c.var)/2)
  
  ## Balance
  um.wt <- bal.data %>% dplyr::group_by(treat) %>% 
    dplyr::summarize(across(all_of(bal.names), ~mean(.x))) %>% as.data.frame()
  
  bal.st <- bal.data %>% dplyr::group_by(treat) %>% 
    dplyr::summarize(across(all_of(bal.names), ~ weighted.mean(.x, wts))) %>% as.data.frame()
  
  um.wt.tab <- matrix(NA, length(bal.names), 3)
  um.wt.tab[,1] <- unlist(um.wt[1,-1]) 
  um.wt.tab[,2] <- unlist(um.wt[2,-1])                        
  um.wt.tab[,3] <- (unlist(um.wt[2,-1]) - unlist(um.wt[1,-1]))/pooled.var
  
  bal.st.tab <- matrix(NA, length(bal.names), 3)
  bal.st.tab[,1] <- unlist(bal.st[1,-1]) 
  bal.st.tab[,2] <- unlist(bal.st[2,-1])                        
  bal.st.tab[,3] <- (unlist(bal.st[2,-1]) - unlist(bal.st[1,-1]))/pooled.var # SMD
  
  um.wt.bias <- um.wt.tab[,3]
  bal.bias <- bal.st.tab[,3] 
  
  ## Bias Reduction
  pbr.bal.wt <- (1 - (mean(abs(bal.bias))/mean(abs(um.wt.bias))))*100
  if (verbose) {
    print("PBR calculated")
  }
  
  # Outcomes for weighting
  bal.wt <- geeglm(Y ~ treat, data = data, std.err = 'san.se', 
                         weights = wts, id=1:nrow(data),
                         corstr="independence")
  att.bw <- msm.out(bal.wt)
  att.bw
  sr.att.bw <- sum((data$treat -(1 - data$treat)*data$wts)*data$Y)/sum(data$treat)
  sr.att.bw # singly robust
  
  with(data, {
    mean(Y[treat==1]) - sum(wts[treat==0] * Y[treat==0]) / sum(wts[treat==0])
  })
  
  sum(data$wts[data$treat == 0])
  sum(data$treat)
  
  bias.bw <- sr.att.bw - true_att
  bias.bw
  # ATT of ACIC-17 around 0.118
  
  coverage.bw <- (true_att >= att.bw["lcl.treat"]) & (true_att <= att.bw["ucl.treat"])
  names(coverage.bw) <- NULL
  
  if (verbose) {
    print("simulation finished")
  }
  
  if (bias.bw > 3) {
    warning("Bias is large for balancing weights. Caution advised.")
  }
  if (type == "l2") {
    est_final = ifelse(simplex, "bw.simp", "bw.nosimp")
  } else {
    est_final = "bw.inf"
  }
  data.frame(est = est_final, "feat_rep" = feat_rep, 
             "bias" = bias.bw, "cvg" = coverage.bw, 
             "pbr" = pbr.bal.wt, "ess" = n.0, "est.att" = sr.att.bw)
}


logisticIPW <- function(data, true_att, feat_rep, verbose = FALSE) {
  
  if ("Z" %in% colnames(data)) {
    data <- data %>% 
      dplyr::rename(treat = Z) #%>% dplyr::select(-Y0, -Y1)
  }
  
  covs <- names(data)[!names(data) %in% c("Y", "treat")]
  basis <- c(covs, "-1")
  X <- scale(model.matrix(reformulate(basis), data))
  trt <- data$treat
  n <- nrow(data)
  
  
  # traditional IPW weights
  psmod <- glm(reformulate(covs, response = "treat"),  family = binomial(), data = data)
  pscore <- psmod$fitted.values
  pscore <- pmax (0.05, pmin (0.95, pscore))
  ip <- ifelse(data$treat == 0, pscore / (1 - pscore), 1) # e(x) / (1 - e(x)) for Z=0               
  data$ip.wts <- ip
  if (verbose) {
    print("IPW weights estimated")
  }
  
  bal.data <- as.data.frame(X)
  bal.names <- names(bal.data)
  n_covs <- length(bal.names)
  bal.data$treat <- data$treat
  bal.data$ip.wts <- data$ip.wts
  
  # ESS
  data0 <- data %>% dplyr::filter(treat==0)
  data1 <- data %>% dplyr::filter(treat==1)
  ## IP Weights ESS
  n.0.ip <- (sum(data0$ip.wts)^2) / (sum(data0$ip.wts^2))
  n.1.ip <- (sum(data1$ip.wts)^2) / (sum(data1$ip.wts^2))
  
  if (verbose) {
    print("ESS calculated")
  }
  
  ## Balance diagnostics
  
  data.var <- bal.data %>% dplyr::group_by(treat) %>% 
    dplyr::summarize(across(all_of(bal.names), ~var(.x))) %>% as.data.frame()
  
  c.var <- as.numeric(data.var[1,])
  t.var <- as.numeric(data.var[2,])
  c.var <- c.var[-1]
  t.var <- t.var[-1]
  pooled.var <- sqrt((t.var + c.var)/2)
  
  ## Balance
  um.wt <- bal.data %>% dplyr::group_by(treat) %>% 
    dplyr::summarize(across(all_of(bal.names), ~mean(.x))) %>% as.data.frame()
  
  bal.ip <- bal.data %>% dplyr::group_by(treat) %>% 
    dplyr::summarize(across(all_of(bal.names), ~ weighted.mean(.x, ip.wts))) %>% as.data.frame()  
  
  um.wt.tab <- matrix(NA, length(bal.names), 3)
  um.wt.tab[,1] <- unlist(um.wt[1,-1]) 
  um.wt.tab[,2] <- unlist(um.wt[2,-1])                        
  um.wt.tab[,3] <- (unlist(um.wt[2,-1]) - unlist(um.wt[1,-1]))/pooled.var
  
  bal.ip.tab <- matrix(NA, length(bal.names), 3)
  bal.ip.tab[,1] <- unlist(bal.ip[1,-1]) 
  bal.ip.tab[,2] <- unlist(bal.ip[2,-1])                        
  bal.ip.tab[,3] <- (unlist(bal.ip[2,-1]) - unlist(bal.ip[1,-1]))/pooled.var
  
  um.wt.bias <- um.wt.tab[,3]
  ip.wt.bias <- bal.ip.tab[,3]
  
  ## Bias Reduction
  pbr.ip.wt <- (1 - (mean(abs(ip.wt.bias))/mean(abs(um.wt.bias))))*100
  pbr.ip.wt
  
  if (verbose) {
    print("PBR calculated")
  }
  
  # Outcomes for weighting
  ipw.wt <- bw <- geeglm(Y ~ treat, data = data, std.err = 'san.se', 
                         weights = ip.wts, id=1:nrow(data),
                         corstr="independence")
  att.ipw <- msm.out(ipw.wt)
  att.ipw
  sr.att.ipw <- sum((data$treat -(1 - data$treat)*data$ip.wts)*data$Y)/sum(data$treat)
  sr.att.ipw # singly robust
  
  bias.ipw <- sr.att.ipw - true_att
  coverage.ipw <- (true_att >= att.ipw["lcl.treat"]) & (true_att <= att.ipw["ucl.treat"])
  names(coverage.ipw) <- NULL
  
  if (verbose) {
    print("simulation finished")
  }
  
  data.frame(est = "ipw", "feat_rep" = feat_rep, 
             "bias" = bias.ipw, "cvg" = coverage.ipw,
             "pbr" = pbr.ip.wt, "ess" = n.0.ip, "est.att" = sr.att.ipw)
}


outcomeRegression <- function(data, true_att, feat_rep, type = "ols", verbose = FALSE) {
  if ("Z" %in% colnames(data)) {
    data <- data %>% 
      dplyr::rename(treat = Z) #%>% dplyr::select(-Y0, -Y1)
  }
  
  covs <- names(data)[!names(data) %in% c("Y", "treat")]
  covs <- c(covs, "-1")
  X <- scale(model.matrix(reformulate(covs), data))
  trt <- data$treat
  n <- nrow(data)
  
  ## Outcome regression
  if (type == "ols") {
    out.reg <- lm(reformulate(c("treat", covs[covs != "-1"]), response = "Y"), data)
    est.att <- coef(out.reg)["treat"]
    names(est.att) <- NULL
    # Compute robust variance-covariance matrix (HC0)
    vcov_hc <- vcovHC(out.reg, type = "HC3")
    # Extract robust standard error for 'treat'
    se_treat <- sqrt(diag(vcov_hc))["treat"]
    ci_lower <- est.att - abs(qnorm(.025)) * se_treat
    ci_upper <- est.att + abs(qnorm(.025)) * se_treat
    bias <- est.att - true_att
    # coverage
    coverage <- true_att >= ci_lower & true_att <= ci_upper
    names(coverage) <- NULL
  } else if (type == "lasso") {
    out.reg <- cv.glmnet(X[trt == 0, ], data$Y[trt == 0], alpha = 1, family = "gaussian", nfolds = 10)
    preds <- predict(out.reg, s="lambda.min", newx = X[trt == 1, ])
    est.att <- mean(data$Y[trt == 1]) - mean(preds)
    bias <- est.att - true_att
    coverage <- NA
  } else if (type == "ridge") {
    out.reg <- cv.glmnet(X[trt == 0, ], data$Y[trt == 0], alpha = 0, family = "gaussian", nfolds = 10)
    preds <- predict(out.reg, s="lambda.min", newx = X[trt == 1, ])
    est.att <- mean(data$Y[trt == 1]) - mean(preds)
    bias <- est.att - true_att
    coverage <- NA
  } else if (type == "rf") {
    input_data <- data
    outcome <- "Y"
    data.0 <- input_data %>% filter(treat == 0) %>% dplyr::select(-treat)
    data.1 <- input_data %>% filter(treat == 1) %>% dplyr::select(-treat)
    rfmod <- randomForest::randomForest(as.formula(paste(outcome, "~ .")), 
                                        data = data.0, ntree=500)
    rf_pred <- predict(rfmod, newdata = data.1)
    est.att <- mean(data.1[[outcome]]) - mean(rf_pred)
    bias <- est.att - true_att
    coverage <- NA
  }
  if (verbose) {
    print("simulation finished")
  }
  
  data.frame(est = type, "feat_rep" = feat_rep, 
             "bias" = bias, "cvg" = coverage, pbr = NA, ess = NA, "est.att" = est.att)
}


augmentedBalWeights <- function(data, true_att, feat_rep, out.mod = "ols", ps.mod = "l2", verbose = FALSE) {
  # out.mod: ols, lasso, ridge
  # ps.mod: sbw, ipw, l2
  # do not mix lasso with l2, or ridge with sbw
  if ("Z" %in% colnames(data)) {
    data <- data %>% 
      dplyr::rename(treat = Z) #%>% dplyr::select(-Y0, -Y1)
  }
  
  covs <- names(data)[!names(data) %in% c("Y", "treat")]
  covs <- c(covs, "-1")
  X <- scale(model.matrix(reformulate(covs), data))
  trt <- data$treat
  n <- nrow(data)
  
  
  ### Propensity Score Estimation ###
  
  if (ps.mod == "inf") {
    bal_weights <- sbw_osqp(X, Z = trt, tol = 0.5)
    data$wts <- pmax(bal_weights$wts, 0)
    data$wts[data$treat == 1] <- 1
    if (verbose) {
      print("L-inf weights estimated")
    }
    if (verbose) {
      print("sbw weights estimated")
    }
  } else if (ps.mod == "ipw") {
    ps.glm <- glm(reformulate(covs[covs != "-1"], response = "treat"), data, family = "binomial")
    pscore <- predict(ps.glm, newdata = data, type = "response")
    pscore <- pmax (0.05, pmin (0.95, pscore))
    data$wts <- ifelse(data$treat == 1, 1, pscore / (1 - pscore))
    if (verbose) {
      print("ipw weights estimated")
    }
  } else if (ps.mod == "l2") {
    # balancing weights
    l <- estimate_regularization(data, covs = covs)
    bal_weights <- multilevel_qp(X, trt, Z = rep(1,n), lambda = l, verbose= FALSE, 
                                 exact_global = FALSE, scale_sample_size = TRUE)
    data$wts <- pmax(bal_weights$weights, 0)
    data$wts[data$treat == 1] <- 1
    if (verbose) {
      print("l2 weights estimated")
    }
  }
  
  
  if (verbose) {
    print("balancing weights estimated")
  }
  
  ### AIPW ###
  #Outcome Model for Control
  if (out.mod == "ols") {
    eta0.glm <- lm(reformulate(covs[covs != "-1"], response = "Y"), subset = treat == 0, data = data)
    x <- model.matrix(reformulate(covs[covs != "-1"]), data)
    mu0 <- predict(eta0.glm, newdata = data.frame(x), type = "response")
    if (verbose) {
      print("DR outcome model estimated")
    }
  } else if (out.mod == "lasso") {
    eta0.glm <- cv.glmnet(X[trt == 0, ], data$Y[trt == 0], alpha = 1, family = "gaussian", nfolds = 10)
    preds <- predict(eta0.glm, s = "lambda.min", newx = X)
    mu0 <- preds
    if (verbose) {
      print("DR outcome model estimated")
    }
  } else if (out.mod == "ridge") {
    eta0.glm <- cv.glmnet(X[trt == 0, ], data$Y[trt == 0], alpha = 0, family = "gaussian", nfolds = 10)
    preds <- predict(eta0.glm, s = "lambda.min", newx = X)
    mu0 <- preds
    if (verbose) {
      print("DR outcome model estimated")
    }
  }
  
  dr.att <- sum((data$treat -(1 - data$treat)*data$wts)*(data$Y - mu0))/sum(data$treat)
  dr.att
  
  mod.dr <- lm(reformulate(c("treat", covs[covs != "-1"]), response="Y"), data=data, weights = wts)
  mod.dr.att <- msm.out(mod.dr)
  mod.dr.att
  if (verbose) print("Doubly robust estimates computed")
  
  
  bias.dr <- dr.att - true_att
  coverage.dr <- true_att >= mod.dr.att["lcl.treat"] & true_att <= mod.dr.att["ucl.treat"]
  names(coverage.dr) <- NULL
  
  if (verbose) {
    print("simulation finished")
  }
  
  data.frame(est = paste(out.mod, ps.mod, sep = "_"), "feat_rep" = feat_rep, 
             bias = bias.dr, cvg = coverage.dr, pbr = NA, ess = NA,
             "true.att" = true_att, est.att = dr.att)
}

