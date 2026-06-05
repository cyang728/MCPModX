#' MCP-ModX: covariate-adjusted MCP-Mod under missing outcomes
#'
#' Fits the model-assisted MCP-Mod analysis of a dose-finding trial.  Marginal
#' arm means are estimated by inverse-probability-weighted g-computation / AIPW
#' (handling outcomes missing at random), inference uses the stacked-sandwich
#' covariance (Theorem 3.3), the MCP step performs the optimal-contrast test
#' with a multiplicity-adjusted critical value, and the Mod step model-averages
#' the candidate dose-response shapes by generalised AIC.
#'
#' @param data a data frame.
#' @param outcome name of the outcome column; missing values coded \code{NA}.
#' @param dose name of the numeric dose column (placebo is the smallest value).
#' @param models named list of candidate shapes and guesstimates passed to
#'   \code{\link{build_models}}, e.g.
#'   \code{list(linear = NULL, emax = c(0.05, 0.2, 0.5), quadratic = -0.85)}.
#' @param covariates character vector of outcome-model adjustment covariates, or
#'   \code{NULL} for the unadjusted (generalised MCP-Mod) analysis.
#' @param obs_covariates character vector of observation-model covariates used
#'   for the IPW propensity, or \code{NULL} for complete-case / complete-data.
#' @param family a \code{\link[stats]{family}}: \code{binomial()} (default,
#'   marginal log-odds estimand) or \code{gaussian()} (mean-difference estimand).
#' @param alpha family-wise type-I error for the MCP step.
#' @param variance \code{"sandwich"} (default, Theorem 3.3) or \code{"oracle"}
#'   (the eq. 3.20 plug-in, for comparison).
#' @param mod logical; run the Mod step (default \code{TRUE}).
#'
#' @return an object of class \code{"MCPModX"}: a list with \code{mu} (marginal
#'   means), \code{delta} and \code{Sigma_delta} (placebo-adjusted contrasts and
#'   their covariance), \code{mcp} (test result), \code{mod} (model-averaging
#'   result), and the candidate-shape matrix \code{models}.
#' @export
#' @examples
#' data(sim_dat)
#'
#' ## unadjusted, complete-case (generalised MCP-Mod)
#' MCPModX(sim_dat, outcome = "y", dose = "dose",
#'         models = list(linear = NULL, emax = c(0.05, 0.2), quadratic = -0.85))
#'
#' ## covariate-adjusted with inverse-probability weighting for the missing
#' ## outcomes (the sandwich variance is used by default)
#' fit <- MCPModX(sim_dat, outcome = "y", dose = "dose",
#'                models = list(linear = NULL, emax = c(0.05, 0.2), quadratic = -0.85),
#'                covariates = c("x1", "x2", "x3"),
#'                obs_covariates = c("x1", "x2"))
#' fit
#' fit$delta            # placebo-adjusted contrasts
#' fit$mod$curve        # model-averaged dose-response curve
MCPModX <- function(data, outcome, dose, models,
                    covariates = NULL, obs_covariates = NULL,
                    family = stats::binomial(), alpha = 0.05,
                    variance = c("sandwich", "oracle"), mod = TRUE) {
  variance <- match.arg(variance)
  if (is.character(family)) family <- get(family, mode = "function")()
  if (is.function(family)) family <- family()

  est <- mcpmodx_estimate(data, outcome, dose, covariates, obs_covariates, family)

  V <- switch(variance,
    sandwich = tryCatch(mcpmodx_sandwich(est), error = function(e) NULL),
    oracle   = mcpmodx_oracle(est))
  if (is.null(V)) V <- mcpmodx_oracle(est)           # robust fallback
  dl <- .delta_from_means(est, V)

  D0 <- build_models(models, est$doses)
  mcp <- mcp_test(dl$delta, dl$Sigma, D0, alpha = alpha)
  modres <- if (mod) mod_average(dl$delta, dl$Sigma, est$doses) else NULL

  structure(list(
    mu = setNames(est$mu_hat, paste0("d", est$doses)),
    delta = setNames(dl$delta, paste0("d", est$doses[-1])),
    Sigma_delta = dl$Sigma, V_mean = V,
    doses = est$doses, models = D0, mcp = mcp, mod = modres,
    family = family$family, variance = variance, alpha = alpha,
    covariates = covariates, obs_covariates = obs_covariates,
    n_missing = sum(est$R == 0L), n = est$n), class = "MCPModX")
}

#' @export
print.MCPModX <- function(x, ...) {
  cat("MCP-ModX fit\n")
  cat(sprintf("  family: %s  |  estimand: %s\n", x$family,
              if (x$family == "binomial") "marginal log-odds" else "marginal mean diff."))
  cat(sprintf("  n = %d  (%d missing outcomes)  |  variance: %s\n",
              x$n, x$n_missing, x$variance))
  cat(sprintf("  covariates: %s  |  obs-model: %s\n",
              if (is.null(x$covariates)) "(none)" else paste(x$covariates, collapse = ", "),
              if (is.null(x$obs_covariates)) "(none)" else paste(x$obs_covariates, collapse = ", ")))
  cat("\n  Placebo-adjusted contrasts (delta):\n")
  print(round(x$delta, 4))
  cat(sprintf("\n  MCP step (alpha = %.3f): %s",
              x$alpha, if (isTRUE(x$mcp$reject == 1)) "SIGNIFICANT" else "not significant"))
  if (isTRUE(x$mcp$reject == 1))
    cat(sprintf("  (best shape: %s)", colnames(x$models)[x$mcp$best]))
  cat(sprintf("\n    z_max = %.3f, crit = %.3f, p = %.4f\n",
              x$mcp$zmax, x$mcp$crit, x$mcp$pvalue))
  if (!is.null(x$mod)) {
    cat("\n  Mod step gAIC weights:\n")
    print(round(x$mod$weights, 3))
  }
  invisible(x)
}
