# this script runs a simulation
# Next step would be to fit the balancing weights and estimate treatment effects
library(foreign)
library(balancer)
library(tidyverse)
library(sandwich)
library(splines)
library(mgcv)
library(aciccomp2017)
library(doParallel)
library(parallel)
library(foreach)
library(doFuture)

write("", file = "sim1-log.txt")  # Clears the file before running

print("Log file cleared")

source("generateACICData.R")
source("fitSplines.R")
source("randomForestFeatures.R")
source("estimationFunctions.R")

cat("Functions loaded \n", file = "sim1-log.txt", append = TRUE)


numCores <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))

plan(multisession, workers = numCores)

cat("Cores initialized \n", file = "sim1-log.txt", append = TRUE)

run_simulation <- function(dgp = 17, seed, verbose = FALSE) {
  if(seed > 250) {
    stop("Seed must be between 1-250.")
  }
  
  ### Generate data ###
  data <- generate_acic_data(dgp = dgp, seed = seed)
  orig_data <- data
  if (verbose) {
    print("data processed")
  }
  ################################################################################################
  
  
  ##### FIT GAM (should probably throw into a different function) #####
  covs <- names(data)[!names(data) %in% c("Y", "W", "alpha", "Y1", "Y0")]
  covs_gam_mod <- covs[grepl("^x_[0-9]+$", covs)]
  gam_fit <- fit_gam(data, covs_gam = covs_gam_mod, cutoff = 3)
  nonlinear_covs <- gam_fit$nonlinear_covs
  if (verbose) {
    print("gam fitted")
  }
  ############################################################################
  
  
  ### Generate splines ###
  spline_obj <- generate_splines(data, nonlinear_covs = nonlinear_covs, df = 3)
  data_spline <- spline_obj$data_spline
  num_nonlin_spline <- spline_obj$num_nonlin
  if (verbose) {
    print("splines generated")
  }
  ############################################################################
  
  
  #### Generate Random Forest features #####
  
  ## First, run random forest
  form <- reformulate(covs, response = "Y")
  rfmod <- randomForest(form, data = data, ntree = 100)
  
  # Next, extract feat reps
  rf_obj <- extract_rf_features(rfmod, data, n_components = 5)
  data_rf <- rf_obj$data_rf
  ############################################################################
  
  
  ### Run the actual simulation ###
  feat_reps <- c("spline", "raw", "rf")
  
  out <- lapply(1:length(feat_reps), function(i){
    fr <- feat_reps[i]
    dataset <- switch(fr, "spline" = data_spline, "raw" = data, "rf" = data_rf)
    
    # compute ATT, balance metrics, error metrics, and so forth for each estimator 
    bw <- balancingWeights(data = dataset, feat_rep = fr, dgp = dgp, verbose = verbose)
    ipw <- logisticIPW(data = dataset, feat_rep = fr, dgp = dgp, verbose = verbose)
    out.reg <- outcomeRegression(data = dataset, feat_rep = fr, dgp = dgp, verbose = verbose)
    aug.bw <- augmentedBalWeights(data = dataset, feat_rep = fr, dgp = dgp, verbose = verbose)
    list(bw = bw, ipw = ipw, out.reg = out.reg, aug.bw = aug.bw)
  })
  names(out) <- feat_reps
  out
}


dgps <- 1:32
num_reps <- 250

scenarios <- expand.grid(dgps, seq_len(num_reps))

cat("Starting sims... \n", file = "sim1-log.txt", append = TRUE)

t_exact <- system.time({
  res <- foreach(i = 1:nrow(scenarios), .options.future = list(seed = TRUE)) %dofuture% {
    if (i %% 100 == 0) {
      log_message <- paste(Sys.time(), "- Starting scenario", i, "out of", nrow(scenarios), "\n")
      # Write to log file
      cat(log_message, file = "sim1-log.txt", append = TRUE)
    }
    dgp <- scenarios[i, 1]
    seed <- scenarios[i, 2]
    sim <- run_simulation(dgp = dgp, seed = seed, verbose = FALSE)
    sim
  }
})

#res

t_exact

save(res, t_exact, file = "simulation1-results.Rdata")


