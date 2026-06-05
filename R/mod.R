## ---------------------------------------------------------------------------
## Mod step: GLS fit of each candidate dose-response shape to (delta_hat, Sigma),
## generalised-AIC model averaging, and the averaged placebo-adjusted curve.
##
##   GLS objective : (delta - f(d;theta))^T Sigma^{-1} (delta - f(d;theta))
##   gAIC          : GLS_min + 2 * (#nonlinear params)
##   weights       : exp(-(gAIC - min)/2), normalised.
## Linear and quadratic are linear in their parameters (closed form); Emax is
## profiled over log(ED50).
## ---------------------------------------------------------------------------

#' Mod-step GLS fit + generalised-AIC model averaging
#' @param delta estimated contrasts (length J).
#' @param Sigma covariance of \code{delta}.
#' @param doses doses including placebo.
#' @param shapes character vector of shapes to average over.
#' @return list(curve = averaged contrasts on active doses, weights, gaic, fits).
#' @keywords internal
mod_average <- function(delta, Sigma, doses,
                        shapes = c("linear", "emax", "quadratic")) {
  doses <- sort(unique(doses)); da <- doses[-1]
  Sin <- tryCatch(solve(Sigma), error = function(e) NULL)
  if (is.null(Sin) || any(!is.finite(delta)))
    return(list(curve = rep(NA_real_, length(da)), weights = NA, gaic = NA,
                fits = NULL))
  gls <- function(pred) { r <- delta - pred; as.numeric(t(r) %*% Sin %*% r) }
  fits <- list()
  if ("linear" %in% shapes) {
    th1 <- as.numeric(crossprod(da, Sin %*% delta) / (t(da) %*% Sin %*% da))
    fits$linear <- list(pred = th1 * da, k = 1)
  }
  if ("quadratic" %in% shapes) {
    G <- cbind(da, da^2)
    thQ <- solve(t(G) %*% Sin %*% G, t(G) %*% Sin %*% delta)
    fits$quadratic <- list(pred = as.numeric(G %*% thQ), k = 2)
  }
  if ("emax" %in% shapes) {
    obj <- function(le) { e <- exp(le); b <- da / (da + e)
      a1 <- as.numeric(crossprod(b, Sin %*% delta) / (t(b) %*% Sin %*% b)); gls(a1 * b) }
    opt <- stats::optimize(obj, lower = log(1e-3), upper = log(max(da) * 5))
    e <- exp(opt$minimum); b <- da / (da + e)
    a1 <- as.numeric(crossprod(b, Sin %*% delta) / (t(b) %*% Sin %*% b))
    fits$emax <- list(pred = a1 * b, k = 2)
  }
  gaic <- sapply(fits, function(f) gls(f$pred) + 2 * f$k)
  wts <- exp(-(gaic - min(gaic)) / 2); wts <- wts / sum(wts)
  curve <- Reduce("+", Map(function(f, w) w * f$pred, fits, wts))
  list(curve = curve, weights = wts, gaic = gaic, fits = fits)
}
