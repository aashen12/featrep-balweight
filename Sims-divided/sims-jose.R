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


out_filename <- paste0("logs/sim-jose-", overlap, "-", format(Sys.time(), "%b-%d-%X-%Y"), ".txt")

write("", out_filename, append = FALSE)   ##### ADDED (overwrite existing file)

make_data <- function(n, overlap = "strong") {
  # Generate ZZ variables from standard normal distribution
  
  sig.ep <- ifelse(overlap == "strong", 100, 30)
  #print(sig.ep)
  sig.123 <- diag(c(2,1,1))
  sig.123[1,2] <- 1; sig.123[1,3] <- -1; sig.123[2,3] <- -0.5;
  sig.123 <- Matrix::forceSymmetric(sig.123)
  beta_coef <- c(1,2,-2,-1,-0.5,1)
  
  X.123 <- as.matrix(mvrnorm(n, mu = rep(0,3), Sigma = sig.123))
  colnames(X.123) <- paste0("X", 1:3)
  X4 <- runif(n,-3,3)
  X5 <- rchisq(n,1)
  X6 <- rbinom(n,1,0.5)
  X <- cbind(X.123, X4, X5, X6)
  Z <- ifelse(X %*% matrix(beta_coef, ncol = 1) + rnorm(n,0,sig.ep) > 0, 1, 0)
  Y <- (X.123[,1] + X.123[,2] + X5)^2 + rnorm(n,0,1)
  
  out.df <- data.frame(X, Z = Z, Y = Y)
  
  return(list(out.df=out.df, true.att = 0))
}

d <- make_data(n, overlap = overlap)
d$out.df

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
    
    bdat.obj  = make_data(1000, overlap = overlap)
    bdat <- bdat.obj$out.df
    
    pilot.dat.obj <- make_data(1000, overlap = overlap) # 1000 typically
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

save(scenarios, file="simulations-jose-temp.RData")

filename <- paste0("results/simulation-jose-overlap-", overlap, "-", format(Sys.time(), "%b-%d-%Y"), ".RData")

save(scenarios, file=filename)
print("saved file")



