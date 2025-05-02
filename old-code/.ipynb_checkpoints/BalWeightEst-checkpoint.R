# Next step would be to fit the balancing weights and estimate treatment effects
library(foreign)
library(balancer)
library(tidyverse)
library(sandwich)
library(splines)

orig_data <- read_csv("data/acic_raw_data.csv")

ground_truth <- orig_data %>% dplyr::select(W, Y0, Y1)

true_att <- mean(ground_truth$Y1[ground_truth$W == 1]) - mean(ground_truth$Y0[ground_truth$W == 1])
true_att

true_ate <- mean(ground_truth$Y1) - mean(ground_truth$Y0)
true_ate

data <- orig_data %>% 
  dplyr::rename(treat = W) #%>% dplyr::select(-Y0, -Y1)

covs <- names(data)[!names(data) %in% c("Y", "treat", "alpha", "Y1", "Y0")]
covs <- c(covs, "-1")
X <- scale(model.matrix(reformulate(covs), data))
trt <- data$treat
n <- nrow(data)

# hyperparameter estimation
estimate_regularization <- function(data) {
  data.c <- data %>% dplyr::filter(treat==0)
  lambda.reg <- lm(reformulate(covs, response="Y"), data=data.c)
  l <- var(lambda.reg$resid)
  l
}

msm.out <- function(obj){
  SE <- sqrt(diag(vcovHC(obj, type="HC0")))[2] # robust standard errors
  beta <- coef(obj)[2]
  lcl <- (beta - abs(qnorm(.025))*SE)
  ucl <- (beta + abs(qnorm(.025))*SE)
  return(c("beta" = beta, "lcl" = lcl, "ucl" = ucl, "SE" = SE))
}

l <- estimate_regularization(data = data)

bal_weights <- multilevel_qp(X, trt, Z = rep(1,n), lambda = l, verbose= TRUE, 
                             exact_global = FALSE, scale_sample_size = TRUE)
data$wts <- pmax(bal_weights$weights, 0)
data$wts[data$treat == 1] <- 1
summary(data$wts)
mean(data$wts==0)

# data$wts <- bal_weights$weights

bal.data <- as.data.frame(X)
bal.names <- names(bal.data)
n_covs <- length(bal.names)
bal.data$treat <- data$treat
bal.data$wts <- data$wts
# bal.data$ip.wts <- data$ip.wts

# ESS
nrow(data)
table(data$treat)


data0 <- data %>% filter(treat==0)
n.0 <- (sum(data0$wts)^2) / (sum(data0$wts^2))
data1 <- data %>% filter(treat==1)
n.1 <- (sum(data1$wts)^2) / (sum(data1$wts^2))

n.0
n.1

## Balance diagnostics

data.var <- bal.data %>% group_by(treat) %>% 
  summarize(across(all_of(bal.names), ~var(.x))) %>% as.data.frame()

c.var <- as.numeric(data.var[1,])
t.var <- as.numeric(data.var[2,])
c.var <- c.var[-1]
t.var <- t.var[-1]
pooled.var <- sqrt((t.var + c.var)/2)

## Balance
um.wt <- bal.data %>% group_by(treat) %>% 
  summarize(across(all_of(bal.names), ~mean(.x))) %>% as.data.frame()

bal.st <- bal.data %>% group_by(treat) %>% 
  summarize(across(all_of(bal.names), ~ weighted.mean(.x, wts))) %>% as.data.frame()

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
# ip.wt.bias <- bal.ip.tab[,3]

## Bias Reduction
pbr.bal.wt <- (1 - (mean(abs(bal.bias))/mean(abs(um.wt.bias))))*100
# pbr.ip.wt <- (1 - (mean(abs(ip.wt.bias))/mean(abs(um.wt.bias))))*100

pbr.bal.wt
# pbr.ip.wt

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


#setwd("~/Dropbox/Bal-Weights/draft/educ/figures")


#pdf("overall-balance.pdf", width=8, height=12, onefile=FALSE, paper="special")
ggplot(data=data.plot, aes(x=,(std.dif), y=covariate, color=factor(contrast)))  + 
  geom_point(size=3) + 
  scale_shape_manual(name= "Contrast", values=c(1,12, 20)) + 
  xlab("Standardized Difference") + ylab("Covariates") +
  scale_y_discrete(limits = rev(levels(data.plot$covariate))) +
  geom_vline(xintercept= 0) +
  geom_vline(xintercept= 0.2, linetype = "dashed") +
  geom_vline(xintercept= -0.2, linetype = "dashed") +
  theme_bw()
#dev.off()  


## Outcome Estimates



# Outcomes
bal.wt <- lm(Y ~ treat, data=data, weights = wts)
att <- msm.out(bal.wt)
att
sr.att <- sum((data$treat -(1 - data$treat)*data$wts)*data$Y)/sum(data$treat)
sr.att # singly robust

bias <- sr.att - true_att
bias
# ATT of ACIC-17 around 0.118

coverage(bal.wt)

lm_ancova <- lm(reformulate(c(covs, "treat"), response = "Y"), data)
summary(lm_ancova)$coef[nrow(summary(lm_ancova)$coef),]

bias_ancova <- summary(lm_ancova)$coef[nrow(summary(lm_ancova)$coef),1] - true_att
bias_ancova

      