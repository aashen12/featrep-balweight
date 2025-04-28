# this script fits the splines, identifying nonlinearities via GAM
library(tidyverse)
library(foreign)
library(sandwich)
library(splines)
library(mgcv)

# df_raw <- read_csv("data/acic_raw_data.csv")
# 
# df <- df_raw

# checking a few things
# table(df$W)
# mean(df$Y[df$W == 1]) - mean(df$Y[df$W == 0])



fit_gam <- function(data, covs_gam, cutoff = 4) {
  # Count unique values for each predictor in df
  num_unique <- sapply(data[, covs_gam], function(x) length(unique(x)))
  
  # Remove predictors with too few unique values
  covs_gam <- covs_gam[num_unique >= 10]
  
  # Automatically generate the formula with s() applied to each covariate
  formula_str <- paste("Y", "~", paste(paste0("s(", covs_gam, ", fx = TRUE)"), collapse = " + "))
  
  # Convert to formula object
  gam_formula <- as.formula(formula_str)
  
  # Fit the GAM model
  non.lin <- gam(gam_formula, data = data)
  
  summ_gam <- summary(non.lin)
  summ_gam
  
  nonlinear_df <- summ_gam$s.table %>% data.frame() %>%  rownames_to_column("predictor") %>% 
    dplyr::arrange(desc(edf), p.value) %>% tibble()
  
  nonlinear_df
  
  nonlinear_covs <- gsub("^s\\((.*)\\)$", "\\1", nonlinear_df$predictor[nonlinear_df$edf > cutoff & nonlinear_df$p.value < 0.05])
  list(nonlinear_covs = nonlinear_covs, gam_mod = non.lin)
}


# extracts covs for splines

# Natural spline feature representation

generate_splines <- function(data, nonlinear_covs, df = 3) {
  if (length(nonlinear_covs) > 1) {
    nonlin_matrix <- tryCatch(
      {data %>% dplyr::select(all_of(nonlinear_covs))},
      error = function(e) {
        print("Nonlin bugged out for > 1 nonlin.")
        print(paste0("Nonlin covs: ", nonlinear_covs))
        return(data)
      }
    )
    spline_bases <- tryCatch({
      Map(function(col, name) {
        #print(paste("Processing column:", name))
        sp <- ns(col, df = df, intercept = FALSE)
        nos <- seq(1:ncol(sp))
        #print(nos)
        colnames(sp) <- paste0(name, "_spline", nos)
        sp
      }, nonlin_matrix, names(nonlin_matrix))
    }, error = function(e) {
      print("Generate splines bugged out for > 1 nonlin.")
      print(paste0("Nonlin covs: ", nonlinear_covs))
      return(data)
    })
    
    spline_data <- tryCatch({
      do.call(cbind, spline_bases) %>% data.frame() %>% tibble()
    }, error = function(e) {
      print("do call cbind bugged out for > 1 nonlin.")
      print(paste0("Nonlin covs: ", nonlinear_covs))
      return(data)
    })
    data_spline <- tryCatch({
      cbind(data, spline_data)
    }, error = function(e) {
      print("Append splines bugged out for > 1 nonlin.")
      print(paste0("Nonlin covs: ", nonlinear_covs))
      return(data)
    })
    return(list(data_spline = data_spline, num_nonlin = length(nonlinear_covs)))
  } else if (length(nonlinear_covs) == 1) {
    #print("1 non lin")
    nonlin_matrix <- data %>% dplyr::select(all_of(nonlinear_covs)) %>% dplyr::pull(!!nonlinear_covs)
    #return(nonlin_matrix)
    sp <- ns(nonlin_matrix, df = df, intercept = FALSE)
    #print("ns run")
    nos <- seq(1:ncol(sp))
    #print(nos)
    colnames(sp) <- paste0(nonlinear_covs, "_spline", nos)
    #print("colnames set")
    spline_bases <- sp
    spline_data <- spline_bases %>% data.frame() %>% tibble()
    #print("spline df created")
    data_spline <- cbind(data, spline_data)
    return(list(data_spline = data_spline, num_nonlin = length(nonlinear_covs)))
  } else {
    data_spline <- data
    return(list(data_spline = data_spline, num_nonlin = length(nonlinear_covs)))
  }
}
