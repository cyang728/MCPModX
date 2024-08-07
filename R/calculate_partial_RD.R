#' Calculate Partial Derivatives of Risk Difference (RD)
#'
#' This function computes the partial derivatives of the Risk Difference (RD)
#' for a given vector of probabilities (`mu`).
#'
#' @param mu A numeric vector of probabilities for each dose level.
#'
#' @return A matrix representing the partial derivatives of the Risk Difference (RD).
#' @export
#'
#' @examples
#' # Example usage:
#' mu <- c(0.2, 0.4, 0.6)
#' calculate_partial_RD(mu)
calculate_partial_RD <- function(mu) {

  # Ensure the input is a numeric vector
  if (!is.numeric(mu)) {
    stop("Input 'mu' must be a numeric vector.")
  }

  ndoses <- length(mu)

  # First row of the derivative matrix
  firstrow <- matrix(rep(-1, ndoses - 1), nrow = 1)

  # Initialize the contrast vector matrix
  contr_vec <- matrix(0, nrow = ndoses - 1, ncol = ndoses - 1)

  # Fill the diagonal with partial derivatives for each dose level
  diag(contr_vec) <- 1

  # Combine the first row with the contrast vector matrix
  contr_vec <- cbind(firstrow, contr_vec)

  return(contr_vec)
}
