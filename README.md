
# MCPModX

<!-- badges: start -->
[![R-CMD-check](https://github.com/cyang728/MCPModX/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cyang728/MCPModX/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of MCPModX package is to incorporate covariate adjustment to enhance the robustness and precision of dose-response analyses. Our approach leverages g-computation to address non-collapsibility and employs model-assisted estimation to ensure valid inferences, even in the presence of model misspecification.

## Installation

You can install the development version of MCPModX like so:

``` r
remotes::install_github("cyang728/MCPModX")
```

## Example

This is a basic example which shows you how to solve a common problem:

### 1. Power and Prediction in MCPMod with Key Covariates for Binary Endpoints.

``` r
library(MCPModX)

## Function for generating simulation data
gen_norm = function(n=1000, p1=10, rho=0.3){
  
  X = matrix(0,nrow=n,ncol=p1)
  X[,1] = rnorm(n)
  
  for(j in 2:p1){
    X[,j] = rho*X[,j-1]+sqrt(1-rho^2)*rnorm(n)
  }
  
  return(X)
}

# Function defining the effect of covariates
custom_fX = function(X){
  return(1*X[,1]-2.1*X[,2]+0.5*X[,3]+0.2*X[,4])
}

# Function defining the treatment effect based on dose
custom_fDose = function(dose){
  return(2 * dose)
}

# Function to generate trial data based on dose and covariates
gen_dat_trial = function(dose, p = 10, seed = 100, 
                         fX_func = custom_fX, 
                         rho = 0.3,
                         fDose_func = custom_fDose){
  
  set.seed(seed)
  n = length(dose)
  X = gen_norm(n = n, p1 = p, rho = rho)
  
  fX = fX_func(X)
  fDose = fDose_func(dose)
  
  fX_combined = fX + fDose
  
  probs = exp(fX_combined) / (1 + exp(fX_combined))
  probs[fX_combined > 100] = 1
  probs[fX_combined < -100] = 0
  y = rbinom(n, size = 1, prob = probs)
  dat = data.frame(X = X, dose = dose, y = y)
  
  colnames(dat)[1:p] = paste0("x", seq(p))
  
  return(dat)
}

# Generate data for the trial
n_dose = 30
dose = rep(c(0.00, 0.05, 0.20, 0.40, 0.70, 1.00), each = n_dose)
dat_cur = gen_dat_trial(dose, seed = 1)

# Apply MCPModX analysis to the generated data
results = MCPModX(data = dat_cur, 
                  covariates = c("x1", "x2", "x3", "x4"))

# Significance of the shape in the MCP Step
results$MCT_table

# Predicted log odds ratio using weighted AIC
results$weighted_estimand
```

### 2. Power and Prediction in MCPMod with Prognostic Score for Binary Endpoints. 

In this context, we assume that historical data can be leveraged to train a prognostic score using nonlinear machine learning methods, such as random forests, based on binary outcomes.

``` r
library(MCPModX)
library(randomForest)

## Function for generating simulation data
gen_norm = function(n=1000, p1=10, rho=0.3){
  
  X = matrix(0,nrow=n,ncol=p1)
  X[,1] = rnorm(n)
  
  for(j in 2:p1){
    X[,j] = rho*X[,j-1]+sqrt(1-rho^2)*rnorm(n)
  }
  
  return(X)
}

# Function defining the effect of covariates for historical data
hist_fX = function(X){
  return(1*X[,1]-2.1*X[,2]+0.5*X[,3]+0.2*X[,4])
}

# Function defining the effect of covariates for current trial data
custom_fX = function(X){
  return(1*X[,1]-2.1*X[,2]+0.5*X[,3]+0.2*X[,4])
}

# Function defining the treatment effect based on dose
custom_fDose = function(dose){
  return(2 * dose)
}

# Function to generate historical data based on dose and covariates
gen_dat_hist = function(n = n_hist, p_hist = 10, seed = 100, 
                        rf_col = 1:10, n_tree = 1000, 
                        rho = 0.3, 
                        fX_func = hist_fX){
  
  set.seed(seed)
  X = gen_norm(n = n, p1 = p_hist, rho = rho)
  fX = fX_func(X)
  probs = exp(fX) / (1 + exp(fX))
  probs[fX > 100] = 1
  probs[fX < -100] = 0
  y = rbinom(n, size = 1, prob = probs)
  dat = data.frame(X = X, y = y)
  
  colnames(dat)[1:p_hist] = paste0("x", seq(p_hist))
  
  rf = randomForest(dat[ ,rf_col], factor(dat[ ,(p_hist + 1)]), 
                    ntree = n_tree)
  
  out = list(dat = dat, rf = rf)
  return(out)
}

# Function to generate trial data based on dose and covariates
gen_dat_trial = function(dose, p = 10, seed = 100, 
                         fX_func = custom_fX, 
                         rho = 0.3,
                         fDose_func = custom_fDose){
  
  set.seed(seed)
  n = length(dose)
  X = gen_norm(n = n, p1 = p, rho = rho)
  
  fX = fX_func(X)
  fDose = fDose_func(dose)
  
  fX_combined = fX + fDose
  
  probs = exp(fX_combined) / (1 + exp(fX_combined))
  probs[fX_combined > 100] = 1
  probs[fX_combined < -100] = 0
  y = rbinom(n, size = 1, prob = probs)
  dat = data.frame(X = X, dose = dose, y = y)
  
  colnames(dat)[1:p] = paste0("x", seq(p))
  
  return(dat)
}

# Generate data for the historical data
n_hist = 1000
prog_covariates = 1:10
dat_hist = gen_dat_hist(n = n_hist, rf_col = prog_covariates, seed=1)
    
# Generate data for the current trial data
n_dose = 30
dose = rep(c(0.00, 0.05, 0.20, 0.40, 0.70, 1.00), each = n_dose)
dat_cur = gen_dat_trial(dose, seed = 1)
x_star= predict(dat_hist$rf, newdata = dat_cur, type = "prob")[,2]
x_star = pmax(pmin(x_star, 0.999), 0.001)
x_star = scale(log(x_star / (1 - x_star)))
dat_cur = data.frame(cbind(x_star, dat_cur))

# Apply MCPModX analysis to the generated data
results = MCPModX(data = dat_cur, 
                  covariates = c("x_star"))

# Significance of the shape in the MCP Step
results$MCT_table

# Predicted log odds ratio using weighted AIC
results$weighted_estimand
```


