## ---------------------------------------------------------------------------
## Covariance of the marginal arm means.
##
## Stacked M-estimation (sandwich) covariance V (Theorem 3.3 / SI).  Stacking the
## observation-, outcome- and mean-moment estimating functions
##   Psi_alpha = U (R - p),
##   Psi_eta   = (R/p) Z (Y - m),                        m = g^{-1}(zeta),
##   Psi_mu,a  = (R/p)(I(A=a)/pi_a)(Y - m_a) + m_a - mu_a,
## the influence function of mu_hat_a is
##   phi^mu = Psi_mu - A_mu_eta A_eta_eta^{-1} Psi_eta
##                   - (A_mu_alpha - A_mu_eta A_eta_eta^{-1} A_eta_alpha)
##                       A_alpha_alpha^{-1} Psi_alpha,
##   V = (1/n) sum_i phi^mu_i phi^mu_i^T.
## The GLM working weight m(1-m) is replaced by family$mu.eta(eta) so the same
## code serves binomial (logit) and gaussian (identity) outcomes.
##
## The oracle plug-in v^orc (eq. 3.20) drops the eta/alpha corrections; it is
## consistent for V^orc >= V and is exposed for comparison only.
## ---------------------------------------------------------------------------

#' Stacked-sandwich covariance of the marginal arm means
#' @param est output of \code{mcpmodx_estimate}.
#' @return a (J+1)x(J+1) matrix V (the covariance of sqrt(n)(mu_hat - mu)).
#' @keywords internal
mcpmodx_sandwich <- function(est) {
  with(est, {
    fam <- family
    dat <- est$dat; n <- est$n; J1 <- est$J1; doses <- est$doses
    R <- est$R; w <- est$w; p_hat <- est$p_hat
    y <- est$y_full; y[is.na(y)] <- 0   # missing rows carry weight w = 0 (0*NA -> 0)
    tt <- stats::delete.response(stats::terms(fit))
    Zobs   <- stats::model.matrix(tt, dat, xlev = fit$xlevels)
    eta_o  <- stats::predict(fit, newdata = dat, type = "link")
    m_obs  <- fam$linkinv(eta_o)
    g_obs  <- fam$mu.eta(eta_o)                       # dm/deta at observed dose
    qeta   <- ncol(Zobs)
    ## design and dm/deta with dose set to each level
    Za <- vector("list", J1); ga <- vector("list", J1)
    for (a in seq_len(J1)) {
      d <- dat; d$.dosef <- factor(doses[a], levels = doses)
      Za[[a]] <- stats::model.matrix(tt, d, xlev = fit$xlevels)
      ga[[a]] <- fam$mu.eta(stats::predict(fit, newdata = d, type = "link"))
    }
    ## outcome-score block
    Psi_eta <- (w * (y - m_obs)) * Zobs
    Psi_eta[!is.finite(Psi_eta)] <- 0
    Aee <- crossprod(Zobs, (w * g_obs) * Zobs) / n
    Aee_inv <- solve(Aee)
    ## AIPW moment and cross-derivative A_mu_eta
    Psi_mu <- matrix(0, n, J1); Ame <- matrix(0, J1, qeta)
    for (a in seq_len(J1)) {
      ind <- as.numeric(dat$.dosef == doses[a])
      Psi_mu[, a] <- w * ind / pi_a[a] * (y - M[, a]) + M[, a] - mu_hat[a]
      Ame[a, ] <- colSums((-(1 - w * ind / pi_a[a]) * ga[[a]]) * Za[[a]]) / n
    }
    Psi_mu[!is.finite(Psi_mu)] <- 0
    phi <- Psi_mu - Psi_eta %*% (Aee_inv %*% t(Ame))
    ## observation-model (alpha) block, only under missingness
    if (has_miss) {
      U <- stats::model.matrix(
        stats::as.formula(paste("~", paste(obs_covariates, collapse = "+"))), dat)
      Psi_a <- (R - p_hat) * U
      Aaa <- crossprod(U, (p_hat * (1 - p_hat)) * U) / n
      Aea <- crossprod((w * (1 - p_hat) * (y - m_obs)) * Zobs, U) / n
      Ama <- matrix(0, J1, ncol(U))
      for (a in seq_len(J1)) {
        ind <- as.numeric(dat$.dosef == doses[a])
        Ama[a, ] <- colSums((w * (1 - p_hat) * ind / pi_a[a] * (y - M[, a])) * U) / n
      }
      Bma <- Ama - Ame %*% (Aee_inv %*% Aea)
      phi <- phi - Psi_a %*% (solve(Aaa) %*% t(Bma))
    }
    crossprod(phi) / n
  })
}

#' Oracle plug-in covariance v^orc (eq. 3.20), for comparison
#' @keywords internal
mcpmodx_oracle <- function(est) {
  with(est, {
    y <- y_full; y[is.na(y)] <- 0   # missing rows carry weight w = 0 (0*NA -> 0)
    Phi <- matrix(0, n, J1)
    for (a in seq_len(J1)) {
      ind <- as.numeric(dat$.dosef == doses[a])
      Phi[, a] <- w * ind / pi_a[a] * (y - M[, a]) + M[, a]
    }
    Phi[!is.finite(Phi)] <- 0
    crossprod(scale(Phi, center = mu_hat, scale = FALSE)) / n
  })
}

## Map the mean-scale covariance V to the placebo-adjusted link-scale contrasts
## delta_a = link(mu_a) - link(mu_0), and return Sigma_delta = D V D^T / n.
.delta_from_means <- function(est, V) {
  fam <- est$family; mu <- est$mu_hat; J1 <- est$J1
  dlink <- 1 / fam$mu.eta(fam$linkfun(mu))            # d link / d mu
  D <- matrix(0, J1 - 1, J1)
  for (a in 2:J1) { D[a - 1, a] <- dlink[a]; D[a - 1, 1] <- -dlink[1] }
  delta <- fam$linkfun(mu)[-1] - fam$linkfun(mu)[1]
  list(delta = delta, Sigma = (D %*% V %*% t(D)) / est$n,
       doses = est$doses)
}
