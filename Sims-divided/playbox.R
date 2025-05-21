make_data <- function(n) {
  # Generate ZZ variables from standard normal distribution
  ZZ <- MASS::mvrnorm(n, mu = rep(0, 10), Sigma = diag(10))
  colnames(ZZ) <- paste0("ZZ", 1:10)
  
  # Define X variables using non-linear transformations of ZZ
  X1 <- exp(ZZ[, 1] / 2)
  X2 <- ZZ[, 2] / (1 + exp(ZZ[, 1]))
  X3 <- (ZZ[, 1] * ZZ[, 3] / 25 + 0.6)^3
  X4 <- (ZZ[, 2] + ZZ[, 4] + 20)^2
  X5 <- ZZ[, 5]
  X6 <- ZZ[, 6]
  X7 <- ZZ[, 7]
  X8 <- ZZ[, 8]
  X9 <- ZZ[, 9]
  X10 <- ZZ[, 10]
  
  # Assemble data frame for covariates
  X <- data.frame(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)
  
  # Define eta_3 using Weierstrass function form
  eta3 <- function(ZZ, a = 2, b = 13, N = 20) {
    Z_tilde <- (ZZ[,2] + ZZ[,4] + ZZ[,6] + ZZ[,8] + ZZ[,10]) / 5
    sapply(Z_tilde, function(x) {
      sum(sapply(0:N, function(n) a^(-n) * cos(b^n * pi * x)))
    })
  }
  
  
  # Compute treatment assignment probabilities
  logits <- -ZZ[, 1] - 0.1 * ZZ[, 4] + eta3(ZZ)
  prob_Z <- exp(logits) / (1 + exp(logits))
  Z <- rbinom(n, 1, prob_Z)
  
  
  # Compute outcome variable Y
  # Assumes ZZ is an n x 10 matrix, and T is a binary vector of length n
  
  interaction_term <- 27.4 * ZZ[,1] + 13.7 * ZZ[,2] + 13.7 * ZZ[,3] + 13.7 * ZZ[,4]
  
  Y0 <- 200 - 0.5 * interaction_term + rnorm(n)
  Y1 <- Y0 + 10 + 1.5 * interaction_term
  Y <- Z*Y1 + (1 - Z)*Y0
  
  
  # Return the full dataset
  out.df <- data.frame(X, Z = Z, Y = Y)
  true.att <- mean(Y1[Z == 1]) - mean(Y0[Z == 1])
  
  return(list(out.df=out.df, true.att = true.att))
}

bdat.obj  = make_data(1000)
bdat <- bdat.obj$out.df

pilot.dat.obj <- make_data(1000)
pilot.dat <- pilot.dat.obj$out.df %>% dplyr::filter(Z == 0)

true.att.bdat <- bdat.obj$true.att

true.att.bdat

edat = eval_data(dat=bdat, pilot.dat=pilot.dat, treat.true = true.att.bdat, verbose = FALSE)

#edat$id = id
#edat
out <- lapply(1:length(edat), function(i) {
  resi <- edat[[i]]
  dplyr::bind_rows(resi)
})
o <- dplyr::bind_rows(out) %>% dplyr::mutate(id = id)


o %>% 
  data.frame() %>% 
  tibble() %>% 
  filter(est == "bw.simp") %>% 
  mutate(bias = abs(bias)) %>% 
  arrange(bias) %>% print(n=90)


o %>% 
  data.frame() %>% 
  tibble() %>% 
  filter(est == "ipw") %>% 
  mutate(bias = abs(bias)) %>% 
  arrange(bias) %>% print(n=90)


raw0_lines <- o %>%
  data.frame() %>%
  tibble() %>%
  mutate(bias = abs(bias)) %>%
  mutate(ncomp = as.numeric(str_extract(feat_rep, "\\d+")),
         feat_rep = str_remove(feat_rep, "_\\d+")) %>%
  filter(feat_rep == "raw") %>%
  mutate(method = recode(est, 
                         bw_l2 = "Bal. Weights - L2", 
                         bw_inf = "Bal. Weights - L-Inf",
                         lasso_inf = "Augmented Linf BW with Lasso",
                         ipw = "Trimmed IP Weights",
                         ols = "OLS Outcome Regression",
                         ols_ipw = "Augmented IP Weights with OLS",
                         ridge_l2 = "Augmented L2 BW with Ridge",
                         bw.simp = "Simplex Balancing Weights",
                         rf = "Random Forest Outcome",
                         bw.nosimp = "Off-Simplex Balancing Weights")) %>%
  select(method, raw0_bias = bias)

o_aug <- o %>%
  data.frame() %>%
  tibble() %>%
  mutate(bias = abs(bias)) %>%
  mutate(ncomp = as.numeric(str_extract(feat_rep, "\\d+")),
         feat_rep = str_remove(feat_rep, "_\\d+")) %>%
  mutate(method = recode(est, 
                         bw_l2 = "Bal. Weights - L2", 
                         bw_inf = "Bal. Weights - L-Inf",
                         lasso_inf = "Augmented Linf BW with Lasso",
                         ipw = "Trimmed IP Weights",
                         ols = "OLS Outcome Regression",
                         ols_ipw = "Augmented IP Weights with OLS",
                         ridge_l2 = "Augmented L2 BW with Ridge",
                         bw.simp = "Simplex Balancing Weights",
                         rf = "Random Forest Outcome",
                         bw.nosimp = "Off-Simplex Balancing Weights")) %>%
  left_join(raw0_lines, by = "method") %>%
  filter(!est %in% c("ipw", "bw.nosimp"))

o_aug %>%
  filter(!feat_rep %in% c("kbal_only")) %>% 
  ggplot(aes(x = ncomp, y = bias, color = factor(feat_rep))) +
  geom_point(position = position_dodge(width = 0.4), size = 4) +
  geom_line(aes(group = feat_rep), position = position_dodge(width = 0.4), linewidth = 1.5) +
  geom_hline(aes(yintercept = raw0_bias), linetype = "dashed", color = "black") +
  facet_wrap(~method, scales = "fixed") +
  theme_bw() +
  labs(x = "Number of PCs", y = "Absolute Bias", title = "", color = "Feat. Rep") +
  theme(text = element_text(size = 20))

