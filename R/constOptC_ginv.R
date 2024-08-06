#' Calculate Constrained Optimal Contrast Using Generalized Inverse
#'
#' This function calculates the optimal contrast under the constraint that the control and
#' the active treatment groups have different signs in the contrast.
#'
#' @param mu A numeric vector of means.
#' @param Sinv A matrix representing the generalized inverse of the covariance matrix. Default is NULL.
#' @param placAdj Logical, whether to adjust for placebo (default is FALSE).
#' @param direction Character, either "increasing" or "decreasing" indicating the direction of the effect.
#'
#' @return A numeric vector representing the constrained optimal contrast.
#' @export
#'
constOptC_ginv = function(mu, Sinv = NULL, placAdj = FALSE, direction){
  ## calculate optimal contrasts under the additional constraint that
  ## the control and the active treatment groups have a different sign
  ## in the contrast
  S = ginv(Sinv) # ugly fix, we should use S as argument
  if(!placAdj){
    k = length(mu)
    CC = cbind(-1,diag(k-1))
    SPa = CC%*%S%*%t(CC)
    muPa = as.numeric(CC%*%mu)
  } else {
    k = length(mu)+1
    SPa = S
    muPa = mu
  }
  ## determine direction of effect
  unContr = ginv(SPa)%*%muPa # unconstrained optimal contrast
  mult = ifelse(direction == "increasing", 1, -1) # 1 increasing, -1 decreasing
  ## prepare call of quadprog::ginv.QP
  D = SPa
  d = rep(0,k-1)
  tA = rbind(muPa,
             mult*diag(k-1))
  A = t(tA)
  bvec = c(1,rep(0,k-1))
  contr = quadprog::ginv.QP(D, d, A, bvec, meq=1)$solution
  contr[abs(contr) < 1e-10] = 0
  if(!placAdj)
    contr = c(-sum(contr), contr)
  contr/sqrt(sum(contr^2))
}
