#' Calculate Variance Estimate for Different Measures
#'
#' This function calculates the variance estimate for different measures (log Odds Ratio, Risk Difference, or log Risk Ratio)
#' based on the provided `mu` vector and the variance-covariance matrix `vr`.
#'
#' @param mu A numeric vector of probabilities (e.g., `mu_hat` values) for each dose level.
#' @param vr A variance-covariance matrix.
#' @param estimate A character string specifying the type of measure.
#'        Options are "logOR" (log Odds Ratio), "RD" (Risk Difference), and "logRR" (log Risk Ratio).
#'        Default is "logOR".
#'
#' @return A matrix representing the variance estimate for the chosen measure.
#' @export
#'
#' @examples
#' # Example usage:
#' muhat <- c(0.2, 0.5, 0.8)
#' vrhat <- diag(3)  # Example variance-covariance matrix
#' calculate_variance_estimate(mu = muhat, vr = vrhat, estimate = "logOR")
calculate_variance_estimate <- function(mu = muhat, vr = vrhat,
                                        estimate = c("logOR", "RD", "logRR")) {

  # Ensure that the estimate is one of the allowed values
  estimate <- match.arg(estimate)

  # Use switch to select the correct partial function based on the estimate
  partial_f <- switch(estimate,
                      logOR = calculate_partial_logOR(mu = mu),
                      RD = calculate_partial_RD(mu = mu),
                      logRR = calculate_partial_logRR(mu = mu))

  # Calculate the variance estimate
  S <- partial_f %*% vr %*% t(partial_f)

  return(S)
}
