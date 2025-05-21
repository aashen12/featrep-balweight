# generate Luke's data from Simulation 1 from Bal-Wgt-Educ
library(MASS)
library(tidyverse)
options(list(dplyr.summarise.inform = FALSE))
library(balancer)
library(geepack)
library(GenericML)
library(WeightIt)
library(glmnet)
library(kbal)

rm(list=ls())

source("../utils-keele.R")
source("../fitSplines.R")
source("../randomForestFeatures.R")
source("../estimationFunctions-keele.R")
source("../shen-eval-funcs.R")

library(doFuture)
library(foreach)
library(doParallel)
library(parallel)
library(furrr)


make_data <- function(n,c, treat.true){
  
  mu=rep(0, 3) 
  var.cov.matrix=matrix(c(2, 1, -1, 1, 1, -0.5, -1, -0.5, 1), nrow=3, ncol=3)
  unif.min=-3
  unif.max=3
  chisq.df=1
  binom.prob=0.5
  C = mvrnorm(n=n, mu=mu, Sigma=var.cov.matrix)
  X1 <- C[,1]
  X2 <- C[,2]
  X3 <- C[,3]
  X4 = runif(n=n, min=unif.min, max=unif.max)
  X5 = rbinom(n=n, size=1, prob=binom.prob)
  
  z.star <- 1.5*X1 + 1.5*X2 - 2*X3 + X4 + 0.5*X5
  Z <- as.numeric((z.star / c + runif(n) - 0.5) > 0)
  Y0 <- (X1 + X2 + X3)^2 + mvrnorm(n=1, mu=rep(0, nrow(C)), Sigma=diag(nrow(C)))
  Y1 <- Y0 + treat.true
  Y <- Z*Y1 + (1 - Z)*Y0
  simDat <- data.frame(cbind(Z, Y, X1, X2, X3, X4, X5))
  return(simDat)
}

data <- make_data(n = 1000, c = 7.5, treat.true = 5) #%>% dplyr::rename(treat = Z)

data.c.train <- make_data(n = 1000, c = 7.5, treat.true = 5) %>% dplyr::filter(Z == 0) #%>% dplyr::rename(treat = Z)

bdat <- data
pilot.dat = data.c.train

table(data$treat)

covs <- rf_fit_covs <- paste0("X", 1:5)

edat = eval_data(dat=bdat, pilot.dat=pilot.dat, treat.true=5, verbose = FALSE)

#edat$id = id
#edat
out <- lapply(1:length(edat), function(i) {
  resi <- edat[[i]]
  dplyr::bind_rows(resi)
})
sim_df <- dplyr::bind_rows(out) %>% tibble()
sim_df











