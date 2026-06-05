## ---------------------------------------------------------------------------
## Core estimation: marginal arm means by IPW g-computation / AIPW.
##
##   observation model : R ~ obs_covariates           (logistic)  -> p_hat
##   working outcome   : Y ~ dose(factor) + covariates (GLM, family),
##                       fitted on observed rows with IPW weights 1/p_hat
##   marginal means    : mu_hat_a = mean_i m_hat_a(X_i)
##                       with the AIPW augmentation built into the influence
##                       function used for the variance.
##
## When there is no missingness (obs_covariates = NULL or R == 1 throughout),
## p_hat == 1 and the estimator reduces to complete-data g-computation.
## When covariates = NULL the estimator reduces to the unadjusted (generalised
## MCP-Mod) arm means.
## ---------------------------------------------------------------------------

EPS_PROB <- 1e-6

#' Estimate marginal arm means (internal core)
#'
#' @param data data frame.
#' @param outcome name of the outcome column; missing values coded NA.
#' @param dose name of the (numeric) dose column.
#' @param covariates character vector of outcome-model covariates, or NULL.
#' @param obs_covariates character vector of observation-model covariates, or
#'   NULL for no missingness handling (complete-case / complete-data).
#' @param family a \code{\link[stats]{family}} object (binomial or gaussian).
#' @return a list with the fitted pieces needed for inference.
#' @keywords internal
mcpmodx_estimate <- function(data, outcome, dose, covariates = NULL,
                             obs_covariates = NULL, family = stats::binomial()) {
  n <- nrow(data)
  y_full <- data[[outcome]]
  dvec <- data[[dose]]
  doses <- sort(unique(dvec))
  J1 <- length(doses)
  pi_a <- as.numeric(table(factor(dvec, levels = doses))) / n

  ## observation indicator R (1 = outcome observed)
  R <- as.integer(!is.na(y_full))
  has_miss <- any(R == 0L) && !is.null(obs_covariates)

  dat <- data
  dat$.dosef <- factor(dvec, levels = doses)
  dat$.R <- R
  dat$.y <- y_full

  ## propensity / IPW weights
  if (!has_miss) {
    p_hat <- rep(1, n)
  } else {
    fo <- stats::as.formula(paste(".R ~", paste(obs_covariates, collapse = "+")))
    p_hat <- stats::fitted(stats::glm(fo, data = dat, family = stats::binomial()))
    p_hat <- pmin(pmax(p_hat, EPS_PROB), 1 - EPS_PROB)
  }
  w <- R / p_hat

  ## working outcome model on observed rows, IPW-weighted
  obs <- which(R == 1L)
  rhs <- if (is.null(covariates)) ".dosef"
         else paste(".dosef +", paste(covariates, collapse = "+"))
  fit <- suppressWarnings(stats::glm(
    stats::as.formula(paste(".y ~", rhs)), data = dat[obs, ],
    family = family, weights = 1 / p_hat[obs]))

  ## g-computation predictions m_hat_a(X_i) for all i and the marginal means
  M <- matrix(0, n, J1)
  for (a in seq_len(J1)) {
    dtmp <- dat; dtmp$.dosef <- factor(doses[a], levels = doses)
    M[, a] <- stats::predict(fit, newdata = dtmp, type = "response")
  }
  mu_hat <- colMeans(M)
  if (family$family == "binomial")
    mu_hat <- pmin(pmax(mu_hat, EPS_PROB), 1 - EPS_PROB)

  list(doses = doses, J1 = J1, pi_a = pi_a, mu_hat = mu_hat,
       M = M, fit = fit, p_hat = p_hat, w = w, R = R, y_full = y_full,
       dat = dat, covariates = covariates, obs_covariates = obs_covariates,
       family = family, has_miss = has_miss, n = n)
}
