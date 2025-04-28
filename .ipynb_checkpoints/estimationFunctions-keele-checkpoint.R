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



balancingWeights <- function(data, true_att, feat_rep, type = "l2", verbose = FALSE) {
  # l2 or inf
  #print(paste("True ATT:", round(true_att, 3)))
  if ("Z" %in% colnames(data)) {
    data <- data %>% 
      dplyr::rename(treat = Z) #%>% dplyr::select(-Y0, -Y1)
  }
  
  covs <- names(data)[!names(data) %in% c("Y", "treat")]
  basis <- c(covs, "-1")
  X <- scale(model.matrix(reformulate(basis), data))
  trt <- data$treat
  n <- nrow(data)
  
  if (type == "l2") {
    # balancing weights
    l <- estimate_regularization(data, covs = covs)
    bal_weights <- multilevel_qp(X, trt, Z = rep(1,n), lambda = l, verbose= FALSE,
                                 exact_global = FALSE, scale_sample_size = FALSE)
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
  
  ## Balance Subset
  rownames(um.wt.tab) <- bal.names
  rownames(bal.st.tab) <- bal.names
  
  lg.un <-  um.wt.tab[which(abs(um.wt.tab[,3]) > 0.2),] # large unweighted (0.2)
  lg.wt <- bal.st.tab[which(abs(um.wt.tab[,3]) > 0.2),]
  
  ### Plots and Tables
  n_covs <- nrow(lg.un)
  # var_names <- covs[covs != "-1"] # use if plotting all covariates
  all(row.names(lg.un) == row.names(lg.wt))
  var_names <- rownames(lg.un)
  
  #### Balance Plot
  data.plot <- c(lg.un[,3], lg.wt[,3])
  data.plot <- as.data.frame(data.plot)
  names(data.plot) <- "std.dif"
  data.plot$contrast <- c(rep(1, n_covs), rep(2, n_covs))
  data.plot$contrast <- factor(data.plot$contrast, levels = c(1,2), labels = c("Unweighted", "Balancing Weights"))
  data.plot$covariate <- as.factor(var_names)

  
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
  
  data.frame(est = paste("bw", type, sep = "_"), "feat_rep" = feat_rep, 
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
  
  ## Balance Subset
  rownames(um.wt.tab) <- bal.names
  rownames(bal.ip.tab) <- bal.names
  
  lg.un <-  um.wt.tab[which(abs(um.wt.tab[,3]) > 0.2),] # large unweighted (0.2)
  lg.ip <- bal.ip.tab[which(abs(um.wt.tab[,3]) > .2),]
  
  ### Plots and Tables
  n_covs <- nrow(lg.un)
  # var_names <- covs[covs != "-1"] # use if plotting all covariates
  all(row.names(lg.un) == row.names(lg.ip))
  var_names <- rownames(lg.un)
  
  #### Balance Plot
  data.plot <- c(lg.un[,3], lg.ip[,3])
  data.plot <- as.data.frame(data.plot)
  names(data.plot) <- "std.dif"
  data.plot$contrast <- c(rep(1, n_covs), rep(2, n_covs))
  data.plot$contrast <- factor(data.plot$contrast, levels = c(1,2), labels = c("Unweighted", "IPW"))
  data.plot$covariate <- as.factor(var_names)
  
  
  #setwd("~/Dropbox/Bal-Weights/draft/educ/figures")
  
  
  #pdf("overall-balance.pdf", width=8, height=12, onefile=FALSE, paper="special")
  # p <- ggplot(data=data.plot, aes(x=,(std.dif), y=covariate, color=factor(contrast)))  + 
  #   geom_point(size=3) + 
  #   scale_shape_manual(name= "Contrast", values=c(1,12, 20)) + 
  #   xlab("Standardized Difference") + ylab("Covariates") +
  #   scale_y_discrete(limits = rev(levels(data.plot$covariate))) +
  #   geom_vline(xintercept= 0) +
  #   geom_vline(xintercept= 0.2, linetype = "dashed") +
  #   geom_vline(xintercept= -0.2, linetype = "dashed") +
  #   theme_bw()
  
  # Outcomes for weighting
  ipw.wt <- bw <- geeglm(Y ~ treat, data = data, std.err = 'san.se', 
                         weights = ip.wts, id=1:nrow(data),
                         corstr="independence")
  att.ipw <- msm.out(ipw.wt)
  att.ipw
  sr.att.ipw <- sum((data$treat -(1 - data$treat)*data$ip.wts)*data$Y)/sum(data$treat)
  sr.att.ipw # singly robust
  
  bias.ipw <- sr.att.ipw - true_att
  coverage.ipw <- true_att >= att.ipw["lcl.treat"] & true_att <= att.ipw["ucl.treat"]
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
    ps <- predict(ps.glm, newdata = data, type = "response")
    data$wts <- ifelse(data$treat == 1, 1, ps / (1 - ps))
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

