#' Calculate Optimal Contrast Using Generalized Inverse
#'
#' This function calculates the optimal contrast for a given mean vector `mu`
#' and the generalized inverse of the covariance matrix `Sinv`.
#'
#' @param mu A numeric vector of means.
#' @param Sinv A matrix representing the generalized inverse of the covariance matrix. Default is NULL.
#' @param placAdj Logical, whether to adjust for placebo (default is FALSE).
#'
#' @return A numeric vector representing the optimal contrast.
#' @export
#'
optC_ginv = function(mu, Sinv = NULL, placAdj = FALSE){
  ## calculate optimal contrast for given mu and Sinv (Sinv = proportional to inv covariance matrix)
  if(!placAdj){
    aux = rowSums(Sinv)  # Sinv %*% 1
    mn = sum(mu * aux)/sum(aux) # formula is: S^(-1)(mu-mu*S^(-1)*1/(1*S^(-1)1)1)
    val = Sinv %*% (mu - mn)
    ## now center so that sum is 0
    ## and standardize to have norm 1
    val = val - sum(val)
  } else { # placAdj = TRUE
    val = Sinv %*% mu
  }
  val/sqrt(sum(val^2))
}


