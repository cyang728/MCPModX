#' Calculate Partial Derivatives of the Log Odds Ratio
#'
#' This function computes the partial derivatives of the log odds ratio (log OR) for a given vector of probabilities (`mu`).
#'
#' @param mu A numeric vector of probabilities (e.g., `mu_hat` values) for each dose level.
#'
#' @return A matrix representing the partial derivatives of the log odds ratio.
#' @export
#'
#' @examples
#' # Example usage:
#' # Assuming `mu_hat` is a vector of probabilities
#' # mu_hat = c(0.01, 0.2, 0.4, 0.7)
#' # partial_derivatives <- calculate_partial_logOR(mu_hat)
calculate_partial_logOR <- function(mu) {

  # Ensure the input is a numeric vector
  if (!is.numeric(mu)) {
    stop("Input 'mu' must be a numeric vector.")
  }

  ndoses <- length(mu)

  # First row of the derivative matrix
  firstrow <- matrix(rep(-1 / (mu[1] * (1 - mu[1])), ndoses - 1))

  # Initialize the contrast vector matrix
  contr_vec <- matrix(0, ndoses - 1, ndoses - 1)

  # Fill the diagonal with partial derivatives for each dose level
  for (i in 2:ndoses) {
    diag(contr_vec)[i-1] <- 1 / (mu[i] * (1 - mu[i]))
  }

  # Combine the first row with the contrast vector matrix
  contr_vec <- cbind(firstrow, contr_vec)

  return(contr_vec)
}
