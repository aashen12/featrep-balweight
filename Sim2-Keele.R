# generate Luke's data from Simulation 1 from Bal-Wgt-Educ
library(MASS)
library(tidyverse)
options(list(dplyr.summarise.inform = FALSE))
library(balancer)
library(geepack)
library(GenericML)
library(WeightIt)
library(glmnet)
library(optweight)
library(kbal)


library(doFuture)
library(foreach)
library(doParallel)
library(parallel)
library(furrr)
numCores <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
plan(multisession, workers = numCores)


write("", file = "sim2-keele.txt")  # Clears the file before running

source("utils-keele.R")
source("fitSplines.R")
source("randomForestFeatures.R")
source("estimationFunctions-keele.R")
source("shen-eval-funcs.R")
print("External scripts loaded")

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

### Sim Run
sim_reps = 1000
set.seed(23967)

scenarios = expand_grid(c = c(1, 2.5, 5, 7.5, 10))

print("Starting simulations...")

run_scenario = function( c ) {
  
  log_message <- paste("Starting Simulation:", c, "at", Sys.time(), "\n")
  cat(log_message, file = "sim2-keele.txt", append = TRUE)
  
  # Run the Simulation              
  reps_qs0 = map( 1:sim_reps, function( id ) {
    if (id %% 20 == 0) cat(paste("Starting simulation", id, "for overlap", c, "at", Sys.time(), "\n"), file = "sim2-keele.txt", append = TRUE)
    bdat  = make_data( 1000, c = c, treat.true=5 )
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
  
  cat(paste0("Sim ", c, " done\n"), file = "sim2-keele.txt", append = TRUE)
  dplyr::bind_rows(reps_qs0)
}


# timing
t_fexact <- system.time({ 
  
  scenarios$res = future_pmap( scenarios, run_scenario )

})

save(scenarios, file="simulation-2-keele.RData")
print("saved file")

# time required for computation, in minutes
time <- t_fexact[['elapsed']]/60
cat("Time in minutes: ", time, "\n")

