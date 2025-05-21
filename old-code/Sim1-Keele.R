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

source("utils-keele.R")
source("fitSplines.R")
source("randomForestFeatures.R")
source("estimationFunctions-keele.R")
source("shen-eval-funcs.R")

library(doFuture)
library(foreach)
library(doParallel)
library(parallel)
library(furrr)
write("", "sim1-keele.txt", append = FALSE)   ##### ADDED (overwrite existing file)


numCores <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
plan(multisession, workers = numCores)

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

### Sim Run
sim_reps = 1000
set.seed(23967)

scenarios = expand_grid(c = c(1, 2.5, 5, 7.5, 10))

run_scenario = function( c ) {
  
  log_message <- paste("Starting Simulation:", c, "at", Sys.time(), "\n")
  cat(log_message, file = "sim1-keele.txt", append = TRUE)
  
  
  # Run the Simulation              
  reps_qs0 = map( 1:sim_reps, function( id ) {
    if (id %% 20 == 0) cat(paste("Starting simulation", id, "for overlap", c, "at", Sys.time(), "\n"), file = "sim1-keele.txt", append = TRUE)
    bdat  = make_data( 1000, c=c, treat.true=5 )
    pilot.dat <- make_data(200, c=c, treat.true=5 ) %>% dplyr::filter(Z == 0)
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

save(scenarios, file="simulation-1-keele-temp.RData")

filename <- paste0("simulation-1-keele-", sim_reps, ".RData")

save(scenarios, file=filename)
print("saved file")


# time required for computation, in minutes
time <- t_fexact[['elapsed']]/60
cat("Time in minutes: ", time, "\n")

