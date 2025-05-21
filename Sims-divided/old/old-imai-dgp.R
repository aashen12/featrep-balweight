make_data <- function(n, treat.true = -0.4) {
  # Helper function to simulate correlated variable
  sample_cor <- function(x, rho) {
    z <- rnorm(length(x))
    rho * scale(x) + sqrt(1 - rho^2) * scale(z)
  }

  # Data generation coefficients
  b <- c(0, 0.8, -0.25, 0.6, -0.4, -0.8, -0.5, 0.7)
  a <- c(-3.85, 0.3, -0.36, -0.73, -0.2, 0.71, -0.19, 0.26)
  g1 <- treat.true

  # Covariates
  w1 <- rnorm(n)
  w2 <- rnorm(n)
  w3 <- rnorm(n)
  w4 <- rnorm(n)
  w5 <- sample_cor(w1, 0.2)
  w6 <- sample_cor(w2, 0.9)
  w7 <- rnorm(n)
  w8 <- sample_cor(w3, 0.2)
  w9 <- sample_cor(w4, 0.9)
  w10 <- rnorm(n)

  # Dichotomize
  w1 <- as.numeric(w1 > mean(w1))
  w3 <- as.numeric(w3 > mean(w3))
  w5 <- as.numeric(w5 > mean(w5))
  w6 <- as.numeric(w6 > mean(w6))
  w8 <- as.numeric(w8 > mean(w8))
  w9 <- as.numeric(w9 > mean(w9))

  # Propensity score
  b0 <- b[1]; b1 <- b[2]; b2 <- b[3]; b3 <- b[4]; b4 <- b[5]
  b5 <- b[6]; b6 <- b[7]; b7 <- b[8]

  ps_linear <- b0 + b1*w1 + b2*w2 + b3*w3 + b4*w4 + b5*w5 + b6*w6 + b7*w7

  # Sim G
  z_trueps <- 1 / (1 + exp(-(ps_linear + b2*w2^2 + b4*w4^2 + b7*w7^2 +
                               b1*0.5*w1*w3 + b2*0.7*w2*w4 + b3*0.5*w3*w5 +
                               b4*0.7*w4*w6 + b5*0.5*w5*w7 + b1*0.5*w1*w6 +
                               b2*0.7*w2*w3 + b3*0.5*w3*w4 + b4*0.5*w4*w5 +
                               b5*0.5*w5*w6)))

  # Treatment assignment
  z <- as.numeric(z_trueps > runif(n))

  # Outcome
  y <- a[1] + a[2]*w1 + a[3]*w2 + a[4]*w3 + a[5]*w4 +
    a[6]*w8 + a[7]*w9 + a[8]*w10 + g1*z + sqrt(0.1)*rnorm(n)

  # Output data.frame
  df <- data.frame(
    X1 = w1, X2 = w2, X3 = w3, X4 = w4, X5 = w5, X6 = w6, X7 = w7,
    X8 = w8, X9 = w9, X10 = w10, Z = z, Y = y
  )

  return(df)
}
