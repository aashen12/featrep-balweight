# Redirect output to a log file
library(purrr)
library(furrr)
library(doFuture)
library(foreach)
library(doParallel)
library(parallel)
sink("parallel_test.txt", append = TRUE)

tryCatch({
  # Load libraries
  
  
  # Define a slow function
  slow_function <- function(x) {
    Sys.sleep(10)  # simulate 5 seconds of work
    x^2
  }
  
  # --- Serial execution ---
  cat("Starting SERIAL pmap...\n")
  t1 <- Sys.time()
  
  serial_result <- pmap(list(1:6), function(x) slow_function(x))
  
  t2 <- Sys.time()
  cat("Serial pmap done.\n")
  cat("Serial duration (minutes): ", as.numeric(difftime(t2, t1, units = "mins")), "\n\n")
  
  # --- Parallel execution ---
  cat("Starting PARALLEL future_pmap...\n")
  
  numCores <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
  plan(multisession, workers = numCores)
  cat("numCores: ", numCores)
  
  t3 <- Sys.time()
  
  parallel_result <- future_pmap(list(1:6), function(x) slow_function(x))
  
  t4 <- Sys.time()
  cat("Parallel future_pmap done.\n")
  cat("Parallel duration (minutes): ", as.numeric(difftime(t4, t3, units = "mins")), "\n")
  sink()
  
}, finally = {
  sink()  # Close sink no matter what
})
