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

args <- commandArgs(trailingOnly = TRUE)

overlap_degree <- args[1]
if (is.na(overlap_degree)) {
  stop("`", args[1], "` is not a valid argument.")
}

sim_num <- as.integer(args[2])


out_filename <- paste0("sim", sim_num, "-keele-", overlap_degree, ".txt")

write("", out_filename, append = FALSE)

numCores <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
plan(multisession, workers = numCores)


if (sim_num == 1) {
  make_data <- function(n,c, treat.true){
    sigma <- matrix(0,5,5)
    diag(sigma) <- 1
    xmat <- mvrnorm(n,c(0,0,0,0,0), sigma)
    X1 <- xmat[,1]
    X2 <- xmat[,2]
    X3 <- xmat[,3]
    X4 <- xmat[,4]
    X5 <- xmat[,5]
    z.star <- 1.5*X1 + 1.5*X2 + .75*(X1*X2)
    Z <- as.numeric((z.star / c + runif(n) - 0.5) > 0)
    Y0 <- X2 + X3 + rnorm(n)
    Y1 <- Y0 + treat.true
    Y <- Z*Y1 + (1 - Z)*Y0
    simDat <- data.frame(cbind(Z, Y, X1, X2, X3, X4, X5))
    return(simDat)
  }
} else if (sim_num == 2) {
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
}



### Sim Run
sim_reps = 1000
set.seed(23967)

if (overlap_degree == "low") {
  scenarios = expand_grid(c = 1)
} else if (overlap_degree == "med") {
  scenarios = expand_grid(c = c(2.5, 5))
} else {
  scenarios = expand_grid(c = c(7.5, 10))
}

run_scenario = function( c ) {
  
  log_message <- paste("Starting Simulation:", c, "at", Sys.time(), "\n")
  cat(log_message, file = out_filename, append = TRUE)
  
  
  # Run the Simulation              
  reps_qs0 = map( 1:sim_reps, function( id ) {
    if (id <= 20) cat(paste("Starting simulation", id, "for overlap", c, "at", Sys.time(), "\n"), file = out_filename, append = TRUE)
    if (id > 975) cat(paste("Starting simulation", id, "for overlap", c, "at", Sys.time(), "\n"), file = out_filename, append = TRUE)
    if (id %% 20 == 0) cat(paste("Starting simulation", id, "for overlap", c, "at", Sys.time(), "\n"), file = out_filename, append = TRUE)
    bdat  = make_data( 1000, c=c, treat.true=5 )
    pilot.dat <- make_data(300, c=c, treat.true=5 ) %>% dplyr::filter(Z == 0)
    edat = eval_data(dat=bdat, pilot.dat=pilot.dat, treat.true=5, verbose = FALSE)
    #edat$id = id
    #edat
    out <- lapply(1:length(edat), function(i) {
      resi <- edat[[i]]
      dplyr::bind_rows(resi)
    })
    dplyr::bind_rows(out) %>% dplyr::mutate(id = id)
  }) 
  #cat("Sim Done")
  dplyr::bind_rows(reps_qs0)
}


# timing
t_fexact <- system.time({ 
  
  scenarios$res = future_pmap( scenarios, run_scenario )
  
})

save(scenarios, file="simulation-", sim_num, "-keele-temp.RData")

filename <- paste0("results/simulation-", sim_num, "-keele-", sim_reps, "-", overlap_degree, ".RData")

save(scenarios, file=filename)
print("saved file")


