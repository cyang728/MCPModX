[![R-CMD-check](https://github.com/cyang728/MCPModX/actions/workflows/r.yml/badge.svg)](https://github.com/cyang728/MCPModX/actions/workflows/r.yml)

# MCPModX

`MCPModX` is an R package for **Phase II dose-ranging trials**.

It extends the standard MCP-Mod framework to answer two practical questions that arise at the end of a dose-ranging study:

1. **Is there evidence of a dose-response signal?**
2. **What does the dose-response curve look like after accounting for baseline patient differences and missing outcomes?**

The package is designed for randomized trials with placebo and multiple active doses. It is especially useful when the endpoint is binary or non-normal, when important baseline covariates are available, or when some outcomes are missing.

---

## 🧭 Contents

* [🎯 Clinical target](#clinical-target)
* [❓ Why this package is needed](#why-this-package-is-needed)
* [⚙️ What MCPModX does](#what-mcpmodx-does)
* [✅ When to use MCPModX](#when-to-use-mcpmodx)
* [⚠️ When not to use MCPModX](#when-not-to-use-mcpmodx)
* [📦 Installation](#installation)
* [🚀 Quick start](#quick-start)
* [🔎 Interpreting the output](#interpreting-the-output)
* [🧪 Example: no missing outcomes](#example-no-missing-outcomes)
* [📊 Example: unadjusted analysis](#example-unadjusted-analysis)
* [🧬 Example: using a prognostic score](#example-using-a-prognostic-score)
* [🧾 Data format](#data-format)
* [📈 Choosing candidate dose-response shapes](#choosing-candidate-dose-response-shapes)
* [🧱 Main assumptions](#main-assumptions)
* [🧠 Technical summary for statisticians](#technical-summary-for-statisticians)
* [📚 Citation](#citation)

---

<a id="clinical-target"></a>

## 🎯 Clinical target

The target of `MCPModX` is the **population-average dose-response curve** in the trial population.

For each dose `a`, the package estimates:

```text
mu_a = average outcome if every patient in the trial population received dose a
```

This means:

* For a **binary endpoint**, `mu_a` is the response probability at dose `a`.
* For a **continuous endpoint**, `mu_a` is the mean outcome at dose `a`.

The treatment effect at an active dose is then compared with placebo:

```text
delta_a = link(mu_a) - link(mu_0)
```

where `mu_0` is the placebo response.

For a binary endpoint with a logit link, this is a **marginal log-odds contrast versus placebo**. Importantly, this is not just the coefficient of dose from a logistic regression. It is the dose-response curve after averaging over the baseline covariate distribution of the trial population.

This is the curve used for both:

* the **proof-of-concept test**, and
* the **dose-response modeling / dose-selection step**.

---

<a id="why-this-package-is-needed"></a>

## ❓ Why this package is needed

Standard generalized MCP-Mod is useful for dose-ranging trials, but two common features of real trials are not fully addressed.

### 1. Baseline covariates are often prognostic

Patients may differ in baseline disease severity, biomarkers, demographics, or other pretreatment variables. These variables can strongly predict the outcome.

Ignoring them can reduce power and make the estimated dose-response curve less precise.

`MCPModX` uses baseline covariates to improve precision while still targeting the population-average dose-response curve.

### 2. Outcomes may be missing

In Phase II studies, some outcomes may be missing because of dropout, non-adherence, delayed follow-up, or administrative reasons.

A complete-case analysis can lose information and may be biased when missingness is related to observed patient characteristics.

`MCPModX` can incorporate a model for the probability that an outcome is observed, so that the dose-response analysis remains valid under a missing-at-random assumption.

---

<a id="what-mcpmodx-does"></a>

## ⚙️ What MCPModX does

At a high level, the package performs three steps.

### Step 1. Estimate the adjusted response at each dose

The package estimates the average response at each dose after adjusting for baseline covariates.

For a binary endpoint, this means estimating the adjusted response probability for placebo and each active dose.

### Step 2. Test for a dose-response signal

The package performs the MCP step of MCP-Mod.

It compares the observed dose-response pattern with several prespecified candidate shapes, such as:

* linear,
* Emax,
* quadratic.

The result is a multiplicity-adjusted test of whether there is any evidence of a non-flat dose-response relationship.

### Step 3. Estimate the dose-response curve

If there is evidence of a dose-response signal, the package performs the Mod step.

It fits or averages candidate dose-response shapes to estimate the placebo-adjusted dose-response curve.

---

<a id="when-to-use-mcpmodx"></a>

## ✅ When to use MCPModX

`MCPModX` is appropriate when you have:

* a randomized Phase II dose-ranging trial,
* placebo and multiple active doses,
* a binary or continuous endpoint,
* baseline covariates measured before treatment,
* possible missing outcomes coded as `NA`,
* prespecified candidate dose-response shapes.

The package is most useful when some baseline covariates are expected to be prognostic, such as baseline disease severity, prior treatment history, risk score, biomarker status, or an externally trained prognostic score.

---

<a id="when-not-to-use-mcpmodx"></a>

## ⚠️ When not to use MCPModX

`MCPModX` is not intended to solve every dose-finding problem.

It is not a substitute for:

* a design for non-randomized observational data,
* a sensitivity analysis for missing-not-at-random outcomes,
* a full utility-based benefit-risk dose optimization method,
* an adaptive dose-escalation design,
* a method for selecting candidate dose-response shapes after looking at the data.

Candidate dose-response shapes should be chosen before the analysis based on clinical and pharmacologic knowledge.

---

<a id="installation"></a>

## 📦 Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("cyang728/MCPModX")
```

Then load the package:

```r
library(MCPModX)
```

If you see a GitHub authentication error such as `HTTP error 401: Bad credentials`, the problem is usually an expired GitHub token on your machine, not the package itself. For a public repository, you can try:

```r
remotes::install_github("cyang728/MCPModX", auth_token = NULL)
```

---

<a id="quick-start"></a>

## 🚀 Quick start

The package includes a small simulated dose-ranging dataset called `sim_dat`.

```r
library(MCPModX)

data(sim_dat)
str(sim_dat)
```

The dataset contains:

* `y`: binary outcome, with some missing values coded as `NA`,
* `dose`: assigned dose,
* `x1`, `x2`, `x3`: baseline covariates.

Define the candidate dose-response shapes:

```r
models <- list(
  linear    = NULL,
  emax      = c(0.05, 0.20),
  quadratic = -0.85
)
```

Run MCP-ModX:

```r
fit <- MCPModX(
  data           = sim_dat,
  outcome        = "y",
  dose           = "dose",
  models         = models,
  covariates     = c("x1", "x2", "x3"),
  obs_covariates = c("x1", "x2"),
  family         = binomial()
)

fit
```

In this example:

* `covariates` are used to improve precision of the dose-response estimates.
* `obs_covariates` are used to model the probability that the outcome is observed.
* `family = binomial()` indicates a binary endpoint.

---

<a id="interpreting-the-output"></a>

## 🔎 Interpreting the output

The fitted object contains the main quantities needed for dose-response decision-making.

```r
fit$mcp$reject
```

Whether the global dose-response test rejects the flat null. A value of `TRUE` means there is evidence of a dose-response signal.

```r
fit$mcp$tstat
fit$mcp$crit
```

The MCP test statistics for each candidate shape and the multiplicity-adjusted critical value.

```r
fit$mu
```

The adjusted marginal response estimate at each dose.

For a binary endpoint, these are adjusted response probabilities.

```r
fit$delta
```

The placebo-adjusted treatment effects.

For a binary endpoint with a logit link, these are marginal log-odds contrasts versus placebo.

```r
fit$mod$curve
```

The estimated placebo-adjusted dose-response curve from the Mod step.

```r
fit$mod$weights
```

The model-averaging weights assigned to the candidate dose-response shapes.

---

<a id="example-no-missing-outcomes"></a>

## 🧪 Example: no missing outcomes

When all outcomes are observed, you can omit `obs_covariates`.

```r
fit_complete <- MCPModX(
  data       = sim_dat,
  outcome    = "y",
  dose       = "dose",
  models     = models,
  covariates = c("x1", "x2", "x3"),
  family     = binomial()
)
```

In this setting, `MCPModX` uses covariate adjustment to improve precision.

---

<a id="example-unadjusted-analysis"></a>

## 📊 Example: unadjusted analysis

To run a generalized MCP-Mod analysis without covariate adjustment, omit both `covariates` and `obs_covariates`.

```r
fit_unadjusted <- MCPModX(
  data    = sim_dat,
  outcome = "y",
  dose    = "dose",
  models  = models,
  family  = binomial()
)
```

This can be useful as a comparison with the adjusted analysis.

---

<a id="example-using-a-prognostic-score"></a>

## 🧬 Example: using a prognostic score

A prognostic score can be used as a single covariate if it was trained before the current trial analysis.

For example, suppose `x_star` is a baseline risk score or predicted response score from historical or external data.

```r
fit_score <- MCPModX(
  data       = sim_dat,
  outcome    = "y",
  dose       = "dose",
  models     = models,
  covariates = "x_star",
  family     = binomial()
)
```

This is often useful when many baseline variables are available but the final dose-response analysis should remain simple and interpretable.

The prognostic score should be defined before analyzing the current trial outcome.

---

<a id="data-format"></a>

## 🧾 Data format

Your trial dataset should have one row per patient.

| Column type            | Description                                                            |
| ---------------------- | ---------------------------------------------------------------------- |
| Outcome                | The clinical endpoint, such as response status or change from baseline |
| Dose                   | The randomized dose assignment                                         |
| Baseline covariates    | Pretreatment variables used for precision adjustment                   |
| Observation covariates | Variables used to explain whether the outcome is observed              |

Missing outcomes should be coded as `NA`.

A typical dataset may look like this:

```r
head(sim_dat)
```

```text
  y  dose    x1    x2    x3
  1  0.00  ...
  0  0.05  ...
 NA  0.20  ...
```

Placebo should usually be coded as dose `0`. Active doses should be numeric and should be on the same scale used when specifying the candidate dose-response models.

---

<a id="choosing-candidate-dose-response-shapes"></a>

## 📈 Choosing candidate dose-response shapes

MCP-Mod requires candidate dose-response shapes to be specified before the analysis.

For example:

```r
models <- list(
  linear    = NULL,
  emax      = c(0.05, 0.20, 0.50),
  quadratic = -0.85
)
```

These shapes represent clinically plausible ways the treatment effect may change with dose.

The goal is not to know the true shape in advance. The goal is to include a reasonable set of possible shapes and then perform a multiplicity-adjusted test across them.

---

<a id="main-assumptions"></a>

## 🧱 Main assumptions

The analysis relies on the following assumptions.

### Randomization

Treatment assignment is randomized, so baseline covariates are not confounders of treatment assignment.

### Prespecified candidate shapes

The candidate dose-response shapes are chosen before looking at the outcome data.

### Missing at random

When outcomes are missing, the probability of observing the outcome may depend on observed variables such as dose and baseline covariates.

However, after conditioning on these observed variables, missingness should not depend on the unobserved outcome itself.

### Positivity

Every patient should have a nonzero probability of having the outcome observed within the levels of the variables used in the observation model.

---

<a id="technical-summary-for-statisticians"></a>

## 🧠 Technical summary for statisticians

`MCPModX` estimates the vector of marginal mean responses across doses and then supplies the corresponding covariance matrix to the MCP and Mod steps.

For each dose `a`, the target is:

```text
mu_a = E[Y(a)]
```

and the placebo-adjusted effect is:

```text
delta_a = g(mu_a) - g(mu_0)
```

where `g` is the link function.

The estimator uses standardization over the empirical baseline covariate distribution. When outcomes are missing, the outcome model is fitted with inverse-probability weights from an observation model. The resulting estimator targets the marginal dose-response curve and can be interpreted as an augmented inverse-probability-weighted standardization estimator.

The covariance matrix accounts for estimation of both the outcome model and the observation model. This covariance is then used in:

1. the MCP step, for optimal contrast testing and multiplicity adjustment;
2. the Mod step, for generalized least-squares fitting and model averaging.

When there is no missingness, the method reduces to a complete-data covariate-adjusted MCP-Mod analysis.

When no covariates are supplied, the method reduces to an unadjusted generalized MCP-Mod analysis.

---

<a id="citation"></a>

## 📚 Citation

If you use `MCPModX`, please cite the MCP-ModX manuscript.

Citation details will be added when the manuscript or preprint is publicly available.
