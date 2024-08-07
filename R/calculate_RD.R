#' Calculate Log Risk Difference
#'
#' This function computes the risk difference (RD) between two probabilities.
#'
#' @param p1 A numeric value representing the probability of the event occurring in group 1.
#' @param p0 A numeric value representing the probability of the event occurring in group 0.
#'
#' @return A numeric value representing the risk difference between `p1` and `p0`.
#' @export
#'
#' @examples
#' # Example usage:
#' # Assuming p1 and p0 are the probabilities for two groups
#' # RD <- calculate_RD(0.7, 0.3)
calculate_RD <- function(p1, p0) {

  # Ensure the inputs are numeric and within the range (0, 1)
  if (!is.numeric(p1) || !is.numeric(p0) || any(c(p1, p0) <= 0 | c(p1, p0) >= 1)) {
    stop("Both 'p1' and 'p0' must be numeric values between 0 and 1 (exclusive).")
  }

  # Calculate the log odds ratio
  RD <- p1 - p0

  return(RD)
}
