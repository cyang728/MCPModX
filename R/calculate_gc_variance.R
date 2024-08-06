#' Calculate Variance of Model-Assisted Estimates for Binary Outcomes
#'
#' This function computes the variance of model-assisted estimates for binary outcomes.
#' The input data must contain columns `dose` and `y` (the binary outcome).
#'
#' @param data A data frame containing the variables `dose` and `y`.
#' @param glm_obj A fitted GLM object, typically created by the `glm` function.
#'                Default is `anovaMod`.
#'
#' @return A variance-covariance matrix.
#' @export
#'
calculate_gc_variance <- function(data, glm_obj) {

  # Ensure the necessary columns are present
  if (!all(c("dose", "y") %in% colnames(data))) {
    stop("Data must contain 'dose' and 'y' columns.")
  }

  # Unique dose levels
  dose_all <- unique(data$dose)

  # Initialize the variance-covariance matrix
  S <- matrix(0, nrow = length(dose_all), ncol = length(dose_all))
  y_resp <- as.numeric(data$y)

  # Loop over each unique dose level
  for (j in seq_along(dose_all)) {

    idx_j <- which(data$dose == dose_all[j])
    pi_mean <- mean(data$dose == dose_all[j])

    # Predict for the j-th dose level
    dat_tmp <- data
    dat_tmp$dose <- dose_all[j]
    mu_t <- predict(glm_obj, dat_tmp, type = "response")

    # Calculate components of the variance
    S_rt <- var(y_resp[idx_j] - mu_t[idx_j])
    Qytt <- cov(y_resp[idx_j], mu_t[idx_j])
    S_mt <- var(mu_t)
    v_tt <- S_rt/pi_mean + 2 * Qytt - S_mt
    S[j, j] <- v_tt

    # Loop over pairs of dose levels
    for (k in 1:j) {
      if (k != j) {
        idx_k <- which(data$dose == dose_all[k])
        dat_tmp$dose <- dose_all[k]
        mu_s <- predict(glm_obj, dat_tmp, type = "response")
        Qyts <- cov(y_resp[idx_j], mu_s[idx_j])
        Qyst <- cov(y_resp[idx_k], mu_t[idx_k])
        Qmts <- cov(mu_t, mu_s)
        S[j, k] <- Qyts + Qyst - Qmts
      }
    }
  }

  # Symmetrize the matrix and scale by the number of observations
  S <- S + t(S)
  diag(S) <- diag(S)/2
  S <- S / nrow(data)

  return(S)
}
