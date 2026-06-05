## ---------------------------------------------------------------------------
## MCP step: optimal-contrast test with a multiplicity-adjusted critical value.
##
## For each candidate shape delta_m^0 the optimal contrast is
##   c_m = Sigma^{-1} delta_m^0,
## and the test statistic z_m = c_m^T delta_hat / sqrt(c_m^T Sigma c_m).  Under
## the global null the (z_1,...,z_K) are jointly normal with correlation
## cov2cor(C^T Sigma C); the critical value is the equicoordinate upper-alpha
## quantile of that distribution (mvtnorm::qmvnorm), reducing the family-wise
## error to alpha.  A Bonferroni value is used as a fallback.
## ---------------------------------------------------------------------------

#' Optimal contrasts for a set of candidate shapes
#' @keywords internal
optimal_contrasts <- function(D0, Sigma) {
  Sin <- solve(Sigma)
  C <- Sin %*% D0
  apply(C, 2, function(c) c / sqrt(sum(c^2)))        # unit-norm (cosmetic)
}

#' MCP optimal-contrast test
#' @param delta estimated placebo-adjusted dose-response contrasts (length J).
#' @param Sigma covariance of \code{delta}.
#' @param D0 standardised candidate-shape matrix (active doses x K).
#' @param alpha family-wise type-I error.
#' @return list(reject, tstat, crit, pvalue, best, zmax).
#' @keywords internal
mcp_test <- function(delta, Sigma, D0, alpha = 0.05) {
  Sin <- tryCatch(solve(Sigma), error = function(e) NULL)
  bad <- is.null(Sin) || any(!is.finite(delta)) || any(!is.finite(Sigma)) ||
         any(diag(Sigma) <= 0)
  if (bad)
    return(list(reject = NA, tstat = NA, crit = NA, pvalue = NA,
                best = NA, zmax = NA))
  C <- Sin %*% D0
  num <- as.numeric(crossprod(C, delta))
  den <- sqrt(pmax(diag(t(C) %*% Sigma %*% C), 0))
  z <- num / den
  if (any(!is.finite(z)))
    return(list(reject = NA, tstat = z, crit = NA, pvalue = NA,
                best = NA, zmax = NA))
  Gam <- cov2cor(t(C) %*% Sigma %*% C)
  crit <- tryCatch(
    mvtnorm::qmvnorm(1 - alpha, tail = "lower.tail", corr = Gam)$quantile,
    error = function(e) stats::qnorm(1 - alpha / ncol(D0)))
  pval <- tryCatch(
    1 - mvtnorm::pmvnorm(lower = -Inf, upper = rep(max(z), ncol(D0)), corr = Gam)[1],
    error = function(e) min(1, ncol(D0) * stats::pnorm(max(z), lower.tail = FALSE)))
  list(reject = as.integer(max(z) > crit),
       tstat = z, crit = crit, pvalue = pval,
       best = if (max(z) > crit) which.max(z) else 0L,
       zmax = max(z))
}
