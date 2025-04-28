# Generate data
library(tidyverse)

# Set a reproducible seed
set.seed(123)

### 1. Generate True Covariates (Z)
# Function to generate a data.frame of four latent covariates (z1, z2, z3, z4)
generate_true_covariates <- function(n) {
  Z <- data.frame(
    z1 = rnorm(n),
    z2 = rnorm(n),
    z3 = rnorm(n),
    z4 = rnorm(n)
  )
  return(Z)
}

### 2. Transform True Covariates to Observed Covariates (X)
# Using the Kang & Schafer transformation:
#   x1 = exp(z1/2)
#   x2 = z2/(1 + exp(z1)) + 10
#   x3 = (z1 * z3/25 + 0.6)^3
#   x4 = (z2 + z4 + 20)^2
transform_covariates <- function(Z) {
  X <- data.frame(
    X1 = exp(Z$z1 / 2),
    X2 = Z$z2 / (1 + exp(Z$z1)) + 10,
    X3 = (Z$z1 * Z$z3 / 25 + 0.6)^3,
    X4 = (Z$z2 + Z$z4 + 20)^2
  )
  return(X)
}

### 3. Generate Treatment Assignment (W)
# The treatment is generated using a logistic model with linear predictor based on Z.
# For example, we use:
#   lp = -z1 + 0.5*z2 - 0.25*z3 - 0.1*z4
# and then compute the probability via the logistic function.
generate_treatment <- function(Z) {
  lp <- -Z$z1 + 0.5 * Z$z2 - 0.25 * Z$z3 - 0.1 * Z$z4
  p <- exp(lp) / (1 + exp(lp))
  W <- rbinom(n = length(p), size = 1, prob = p)
  return(W)
}

### 4. Generate Outcome (Y)
# The outcome is generated from a linear model of the form:
#   Y = 210 + 27.4*z1 + 13.7*z2 + 13.7*z3 + 13.7*z4 + error,
# where error is normally distributed.
generate_outcome <- function(Z, error.sd = 1) {
  Y <- 210 + 27.4 * Z$z1 + 13.7 * Z$z2 + 13.7 * Z$z3 + 13.7 * Z$z4 +
    rnorm(n = nrow(Z), mean = 0, sd = error.sd)
  return(Y)
}

### 5. Combine into a Modular Data Generation Function
# This function brings together all components. It returns a data.frame that
# includes the observed covariates X, treatment W, and outcome Y.
generate_dataset <- function(n, error.sd = 1) {
  Z <- generate_true_covariates(n)
  X <- transform_covariates(Z)
  W <- generate_treatment(Z)
  Y <- generate_outcome(Z, error.sd)
  
  # Combine into one data.frame
  data <- data.frame(X, W = W, Y = Y)
  return(data)
}

### Example Usage:
# Generate a dataset with 1000 observations
sim_data <- generate_dataset(n = 200, error.sd = 1)

# write to csv called "ks_data.csv"
write_csv(sim_data, "data/kang_data.csv")


