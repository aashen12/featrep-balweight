rm(list=ls())
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


out_filename <- paste0("logs/sim-imai-", format(Sys.time(), "%b-%d-%X-%Y"), ".txt")

write("", out_filename, append = FALSE)   ##### ADDED (overwrite existing file)

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


numCores <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
plan(multisession, workers = numCores)

### Sim Run
sim_reps = 1000
set.seed(23967)


run_scenario = function() {
  
  log_message <- paste("Starting Kosuke Simulation at", Sys.time(), "\n")
  cat(log_message, file = out_filename, append = TRUE)
  
  
  # Run the Simulation              
  reps_qs0 = mclapply( 1:sim_reps, function( id ) {
    if (id < 20) cat(paste("Starting simulation", id, "at", Sys.time(), "\n"), file = out_filename, append = TRUE)
    if (id > 975) cat(paste("Starting simulation", id, "at", Sys.time(), "\n"), file = out_filename, append = TRUE) 
    if (id %% 20 == 0) cat(paste("Starting simulation", id, "at", Sys.time(), "\n"), file = out_filename, append = TRUE)
    
    bdat.obj  = make_data(1000)
    bdat <- bdat.obj$out.df
    
    pilot.dat.obj <- make_data(1000) # 1000 typically
    pilot.dat <- pilot.dat.obj$out.df %>% dplyr::filter(Z == 0)
    
    true.att.bdat <- bdat.obj$true.att
    
    edat = eval_data(dat=bdat, pilot.dat=pilot.dat, treat.true = true.att.bdat, verbose = FALSE)
    
    #edat$id = id
    #edat
    out <- lapply(1:length(edat), function(i) {
      resi <- edat[[i]]
      dplyr::bind_rows(resi)
    })
    dplyr::bind_rows(out) %>% dplyr::mutate(id = id)
  }, mc.set.seed = TRUE, mc.cores = numCores - 1) 
  #cat("Sim Done")
  dplyr::bind_rows(reps_qs0)
}


scenarios <- run_scenario()

save(scenarios, file="simulations-imai-temp.RData")

filename <- paste0("results/simulation-imai-", format(Sys.time(), "%b-%d-%Y"), ".RData")

save(scenarios, file=filename)
print("saved file")



