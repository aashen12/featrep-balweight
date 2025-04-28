# This script will generate feature representations
library(tidyverse)
library(KRLS)
library(splines)

df_raw <- read_csv("data/acic_raw_data.csv")

df <- df_raw

# checking a few things
table(df$W)
mean(df$Y[df$W == 1]) - mean(df$Y[df$W == 0])

covs <- names(df)[!names(df) %in% c("Y", "W", "alpha", "Y1", "Y0")]

Y_on_X <- reformulate(covs, response = "Y")
Y_on_all <- reformulate(covs, response = "Y")

lm_Y_on_X <- lm(Y_on_X, data = df)
lm_Y_on_all <- lm(Y_on_all, data = df)
summary(lm_Y_on_X)
summary(lm_Y_on_all)

X <- model.matrix(lm_Y_on_X)[, -1]

# Gaussian kernel feature representation
X_gaussian <- gausskernel(X, sigma = 1) %>% data.frame() %>% tibble() # this is my phi(X)

df_gaussian <- data.frame(W = df$W, Y = df$Y, Y0 = df$Y0, Y1 = df$Y1, X_gaussian)
head(df_gaussian)
write_csv(df_gaussian, "data/gauss_kernel_acic.csv")

# Natural spline feature representation
X_ns <- ns(X, df = 5) %>% data.frame() %>% tibble() # this is my phi(X)

df_ns <- data.frame(W = df$W, Y = df$Y, Y0 = df$Y0, Y1 = df$Y1, X_ns)
head(df_ns)
write_csv(df_ns, "data/natural_spline_acic.csv")
