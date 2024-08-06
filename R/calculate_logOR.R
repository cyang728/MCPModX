#' Calculate Log Odds Ratio
#'
#' This function computes the log odds ratio (log OR) between two probabilities.
#'
#' @param p1 A numeric value representing the probability of the event occurring in group 1.
#' @param p0 A numeric value representing the probability of the event occurring in group 0.
#'
#' @return A numeric value representing the log odds ratio between `p1` and `p0`.
#' @export
#'
#' @examples
#' # Example usage:
#' # Assuming p1 and p0 are the probabilities for two groups
#' # logOR <- calculate_logOR(0.7, 0.3)
calculate_logOR <- function(p1, p0) {

  # Ensure the inputs are numeric and within the range (0, 1)
  if (!is.numeric(p1) || !is.numeric(p0) || any(c(p1, p0) <= 0 | c(p1, p0) >= 1)) {
    stop("Both 'p1' and 'p0' must be numeric values between 0 and 1 (exclusive).")
  }

  # Calculate the log odds ratio
  logOR <- log((p1 / (1 - p1)) / (p0 / (1 - p0)))

  return(logOR)
}
