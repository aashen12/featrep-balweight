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

### 4. Generate Potential Outcomes (Y0 and Y1)
# Here, tau is the constant treatment effect.
generate_potential_outcomes <- function(Z, error.sd = 1, tau = -20) {
  # Baseline outcome under control
  Y0 <- 210 + 27.4 * Z$z1 + 13.7 * Z$z2 + 13.7 * Z$z3 + 13.7 * Z$z4 +
    rnorm(n = nrow(Z), mean = 0, sd = error.sd)
  # Outcome under treatment: add the treatment effect
  Y1 <- Y0 + tau
  list(Y0 = Y0, Y1 = Y1)
}

### 5. Combine into a Modular Data Generation Function
generate_dataset <- function(n, error.sd = 1, tau = -20) {
  Z <- generate_true_covariates(n)
  X <- transform_covariates(Z)
  W <- generate_treatment(Z)
  outcomes <- generate_potential_outcomes(Z, error.sd, tau)
  
  # Observed outcome: if treated, use Y1; if control, use Y0.
  Y <- ifelse(W == 1, outcomes$Y1, outcomes$Y0)
  
  # Return all components so we have the potential outcomes available
  data.frame(X, W = W, Y = Y, Y0 = outcomes$Y0, Y1 = outcomes$Y1)
}

### Example Usage:
sim_data <- generate_dataset(n = 1000, error.sd = 1, tau = 20)

write_csv(sim_data, "data/kang_data_GT.csv")

# Compute the true ATT using the stored potential outcomes
true_ATT <- mean(sim_data$Y1[sim_data$W == 1] - sim_data$Y0[sim_data$W == 1])
print(true_ATT)  # This should be close to -20 if tau is constant.