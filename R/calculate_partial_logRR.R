#' Calculate Partial Derivatives of Log Risk Ratio (logRR)
#'
#' This function computes the partial derivatives of the Log Risk Ratio (logRR)
#' for a given vector of probabilities (`mu`).
#'
#' @param mu A numeric vector of probabilities for each dose level.
#'
#' @return A matrix representing the partial derivatives of the Log Risk Ratio (logRR).
#' @export
#'
#' @examples
#' # Example usage:
#' mu <- c(0.2, 0.4, 0.6)
#' calculate_partial_logRR(mu)
calculate_partial_logRR <- function(mu) {

  # Ensure the input is a numeric vector
  if (!is.numeric(mu)) {
    stop("Input 'mu' must be a numeric vector.")
  }

  ndoses <- length(mu)

  # First row of the derivative matrix
  firstrow <- matrix(rep(-1 / mu[1], ndoses - 1), ndoses - 1)

  # Initialize the contrast vector matrix
  contr_vec <- matrix(0, nrow = ndoses - 1, ncol = ndoses - 1)

  # Fill the diagonal with partial derivatives for each dose level
  diag(contr_vec) <- 1 / mu[-1]

  # Combine the first row with the contrast vector matrix
  contr_vec <- cbind(firstrow, contr_vec)

  return(contr_vec)
}
