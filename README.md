
# MCPModX

<!-- badges: start -->
[![R-CMD-check](https://github.com/cyang728/MCPModX/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cyang728/MCPModX/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of MCPModX is to ...

## Installation

You can install the development version of MCPModX like so:

``` r
remotes::install_github("cyang728/MCPModX")
```

## Example

This is a basic example which shows you how to solve a common problem:

### Power and Prediction in MCPMod with Key Covariates for Binary Endpoints.

``` r
library(MCPModX)

## function for generating simulaiton data
gen_norm = function(n=1000, p1=10, rho=0.3){
  
  X = matrix(0,nrow=n,ncol=p1)
  X[,1] = rnorm(n)
  
  for(j in 2:p1){
    X[,j] = rho*X[,j-1]+sqrt(1-rho^2)*rnorm(n)
  }
  
  return(X)
}

# covaraites effect 
custom_fX = function(X){
  return(1*X[,1]-2.1*X[,2]+0.5*X[,3]+0.2*X[,4])
}

# treatment effects
custom_fDose = function(dose){
  return(1 * dose)
}

dat_trial = function(dose, p = 10, seed = 100, 
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

# Generate data 
n_dose = 30
dose = rep(c(0.00, 0.05, 0.20, 0.40, 0.70, 1.00), each = n_dose)
dat_cur = dat_trial(dose, seed = 1)

results = MCPModX(data = dat_cur, 
                  covariates = c("x1", "x2", "x3", "x4"))

results$test_significant
results$weighted_logOR
```

