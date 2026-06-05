# MCPModX

**MCP-Mod with model-assisted covariate adjustment under missing outcomes.**

The goal of `MCPModX` is to incorporate covariate adjustment to enhance the
**robustness and precision** of dose-response analyses. The approach leverages
**g-computation** to address non-collapsibility of the binary log-odds estimand
and employs **model-assisted (augmented inverse-probability-weighted)
estimation** to ensure valid inference even when the working outcome model is
misspecified or when outcomes are **missing at random**. Inference uses the
**stacked M-estimation (sandwich) covariance**, which accounts for the
estimation of the observation- and outcome-model nuisance parameters
(Theorem 3.3 of the manuscript) and remains valid under the heteroscedasticity
intrinsic to binary outcomes.

## What it does

| Stage | Implementation |
|---|---|
| **Estimand** | marginal placebo-adjusted dose-response curve `delta_a = link(mu_a) - link(mu_0)`; for binary outcomes the marginal **log-odds** (collapsibility-correct), for continuous outcomes the mean difference |
| **Means** | IPW g-computation / AIPW: logistic observation model `p_hat`, IPW-weighted working outcome GLM, then standardization `mu_hat_a = mean_i m_hat_a(X_i)` |
| **Variance** | stacked-sandwich `V` (default; consistent whether or not either working model is correct) or the oracle plug-in (eq. 3.20, for comparison) |
| **MCP step** | optimal contrasts `c_m = Sigma^{-1} delta_m^0`, multiplicity-adjusted critical value from the joint normal (`mvtnorm::qmvnorm`) |
| **Mod step** | GLS fit of linear / Emax / quadratic shapes with generalised-AIC model averaging |

When there is no missingness the estimator reduces to complete-data
g-computation; when `covariates = NULL` it reduces to the unadjusted
(generalised MCP-Mod) analysis.

## Installation

You can install the development version of `MCPModX` like so:

```r
# install.packages("remotes")
remotes::install_github("cyang728/MCPModX")
```

or, from a local copy of the package directory:

```r
install.packages("mvtnorm")
R CMD INSTALL MCPModX        # or devtools::install("MCPModX")
```

The package depends only on `stats` and `mvtnorm`; `randomForest` and `MASS`
are used only by the examples below.

## Example 1 — Covariate adjustment for a binary endpoint

This is a basic example showing how MCP-ModX uses key prognostic covariates to
sharpen the dose-response test for a binary outcome.

```r
library(MCPModX)

## Function for generating correlated covariates
gen_norm = function(n = 1000, p1 = 10, rho = 0.3) {
  X = matrix(0, nrow = n, ncol = p1)
  X[, 1] = rnorm(n)
  for (j in 2:p1) X[, j] = rho * X[, j - 1] + sqrt(1 - rho^2) * rnorm(n)
  X
}

# Effect of the covariates on the outcome
custom_fX    = function(X)    1 * X[, 1] - 2.1 * X[, 2] + 0.5 * X[, 3] + 0.2 * X[, 4]
# Treatment effect as a function of dose
custom_fDose = function(dose) 2 * dose

# Generate a trial data set given a dose vector
gen_dat_trial = function(dose, p = 10, seed = 100, rho = 0.3,
                         fX_func = custom_fX, fDose_func = custom_fDose) {
  set.seed(seed)
  n  = length(dose)
  X  = gen_norm(n = n, p1 = p, rho = rho)
  lp = fX_func(X) + fDose_func(dose)
  y  = rbinom(n, size = 1, prob = plogis(lp))
  dat = data.frame(X = X, dose = dose, y = y)
  colnames(dat)[1:p] = paste0("x", seq(p))
  dat
}

# Generate the current trial: 6 doses, 30 patients per arm
n_dose  = 30
dose    = rep(c(0.00, 0.05, 0.20, 0.40, 0.70, 1.00), each = n_dose)
dat_cur = gen_dat_trial(dose, seed = 1)

# Candidate dose-response shapes (guesstimates on the standardized dose scale)
models = list(linear = NULL, emax = c(0.05, 0.20, 0.50), quadratic = -0.85)

# Apply the MCP-ModX analysis, adjusting for the prognostic covariates
fit = MCPModX(data       = dat_cur,
              outcome    = "y",
              dose       = "dose",
              models     = models,
              covariates = c("x1", "x2", "x3", "x4"),
              family     = binomial())

fit                         # one-line summary: significance, best shape, weights

# Significance of each candidate shape in the MCP step
fit$mcp$tstat               # optimal-contrast statistics (one per shape)
fit$mcp$crit                # multiplicity-adjusted critical value
fit$mcp$reject              # 1 if a dose-response signal is detected

# Predicted placebo-adjusted log-odds curve via gAIC model averaging
fit$mod$curve               # averaged contrasts on the active doses
fit$mod$weights             # gAIC model-averaging weights
```

## Example 2 — Prognostic score from historical data (random forest)

Here we assume **historical data** are available to train a prognostic score
with a nonlinear machine-learning method (a random forest) on the binary
outcome. The trained model is applied to the current trial to produce a single
strong covariate `x_star`, which is then passed to `MCPModX`.

```r
library(MCPModX)
library(randomForest)

gen_norm = function(n = 1000, p1 = 10, rho = 0.3) {
  X = matrix(0, nrow = n, ncol = p1)
  X[, 1] = rnorm(n)
  for (j in 2:p1) X[, j] = rho * X[, j - 1] + sqrt(1 - rho^2) * rnorm(n)
  X
}
hist_fX      = function(X)    1 * X[, 1] - 2.1 * X[, 2] + 0.5 * X[, 3] + 0.2 * X[, 4]
custom_fX    = function(X)    1 * X[, 1] - 2.1 * X[, 2] + 0.5 * X[, 3] + 0.2 * X[, 4]
custom_fDose = function(dose) 2 * dose

# Historical data -> train a random forest prognostic model
gen_dat_hist = function(n = 1000, p_hist = 10, seed = 100, rf_col = 1:10,
                        n_tree = 1000, rho = 0.3, fX_func = hist_fX) {
  set.seed(seed)
  X = gen_norm(n = n, p1 = p_hist, rho = rho)
  y = rbinom(n, size = 1, prob = plogis(fX_func(X)))
  dat = data.frame(X = X, y = y)
  colnames(dat)[1:p_hist] = paste0("x", seq(p_hist))
  rf = randomForest(dat[, rf_col], factor(dat[, p_hist + 1]), ntree = n_tree)
  list(dat = dat, rf = rf)
}

gen_dat_trial = function(dose, p = 10, seed = 100, rho = 0.3,
                         fX_func = custom_fX, fDose_func = custom_fDose) {
  set.seed(seed)
  n  = length(dose)
  X  = gen_norm(n = n, p1 = p, rho = rho)
  y  = rbinom(n, size = 1, prob = plogis(fX_func(X) + fDose_func(dose)))
  dat = data.frame(X = X, dose = dose, y = y)
  colnames(dat)[1:p] = paste0("x", seq(p))
  dat
}

# Train the prognostic model on historical data
dat_hist = gen_dat_hist(n = 1000, rf_col = 1:10, seed = 1)

# Current trial, then score it with the historical random forest
n_dose  = 30
dose    = rep(c(0.00, 0.05, 0.20, 0.40, 0.70, 1.00), each = n_dose)
dat_cur = gen_dat_trial(dose, seed = 1)

x_star  = predict(dat_hist$rf, newdata = dat_cur, type = "prob")[, 2]
x_star  = pmax(pmin(x_star, 0.999), 0.001)          # avoid 0/1
x_star  = scale(qlogis(x_star))                     # logit, then standardize
dat_cur = data.frame(x_star = x_star, dat_cur)

# A single, strong prognostic covariate is often all that is needed
models = list(linear = NULL, emax = c(0.05, 0.20, 0.50), quadratic = -0.85)
fit = MCPModX(data       = dat_cur,
              outcome    = "y",
              dose       = "dose",
              models     = models,
              covariates = "x_star",
              family     = binomial())

fit$mcp$tstat               # significance of each shape in the MCP step
fit$mod$curve               # predicted log-odds curve via weighted AIC
```

## Example 3 — Missing outcomes (inverse-probability weighting)

The distinguishing feature of `MCPModX` is valid inference when outcomes are
**missing at random**. Code missing responses as `NA` and supply
`obs_covariates`, the covariates of the logistic observation (propensity)
model; the analysis then uses IPW g-computation with the sandwich variance.

```r
library(MCPModX)

# ... reuse gen_norm / custom_fX / custom_fDose / gen_dat_trial from Example 1
dose    = rep(c(0.00, 0.05, 0.20, 0.40, 0.70, 1.00), each = 30)
dat_cur = gen_dat_trial(dose, seed = 1)

# Make outcomes missing at random: P(observed) depends on x1
set.seed(1)
p_obs = plogis(1.2 + 0.8 * dat_cur$x1)
dat_cur$y[rbinom(nrow(dat_cur), 1, p_obs) == 0] = NA
mean(is.na(dat_cur$y))      # ~ missing fraction

models = list(linear = NULL, emax = c(0.05, 0.20, 0.50), quadratic = -0.85)

# Complete-case generalized MCP-Mod (no adjustment) -- for comparison
fit_cc  = MCPModX(dat_cur, "y", "dose", models, family = binomial())

# MCP-ModX: outcome-model covariates + IPW observation model
fit_ipw = MCPModX(dat_cur, "y", "dose", models,
                  covariates     = c("x1", "x2", "x3", "x4"),  # outcome model
                  obs_covariates = "x1",                       # IPW propensity
                  family         = binomial(),
                  variance       = "sandwich")                 # Theorem 3.3

fit_ipw
fit_ipw$mcp$tstat           # IPW-adjusted optimal-contrast statistics
fit_ipw$delta               # placebo-adjusted contrasts (marginal log-odds)
fit_ipw$Sigma_delta         # their stacked-sandwich covariance
```

## Built-in example data

A small simulated trial with a binary outcome, three covariates, and ~21%
missing responses ships with the package:

```r
data(sim_dat)
str(sim_dat)
MCPModX(sim_dat, "y", "dose",
        models = list(linear = NULL, emax = c(0.05, 0.2), quadratic = -0.85),
        covariates = c("x1", "x2", "x3"), obs_covariates = c("x1", "x2"))
```

## The fitted object

`MCPModX()` returns an object of class `"MCPModX"`; the components most users
need are:

| Component | Meaning |
|---|---|
| `fit$mcp` | MCP step: `tstat` (per-shape optimal-contrast statistics), `crit` (multiplicity-adjusted critical value), `reject`, `pvalue`, `best`, `zmax` |
| `fit$mod` | Mod step: `curve` (gAIC-averaged placebo-adjusted dose-response), `weights`, `gaic` |
| `fit$delta` | placebo-adjusted contrasts (marginal log-odds for binary, mean difference for continuous) |
| `fit$Sigma_delta` | stacked-sandwich covariance of `fit$delta` |
| `fit$mu` | marginal arm means |

## Reproducing the manuscript simulations

See [`../MCPModX_simulation/`](../MCPModX_simulation/) for the full simulation
driver: data generation, all scenarios, the marginal-truth Monte-Carlo
integration, and the parallel grid runner.

## Testing

```r
devtools::test("MCPModX")    # requires the 'testthat' package
```
