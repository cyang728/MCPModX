#' Calculate Mean Response Estimates for Each Dose Level
#'
#' This function computes the mean response estimates (`mu_hat`) for each dose level using a fitted GLM object.
#' The input data must contain a `dose` column.
#'
#' @param data A data frame containing the `dose` column and any other necessary predictors used in the GLM.
#' @param glm_obj A fitted GLM object, typically created by the `glm` function.
#'                Default is `anovaMod`.
#'
#' @return A numeric vector of mean response estimates (`mu_hat`) for each dose level.
#' @export
#'
calculate_gc_muhat <- function(data, glm_obj) {

  # Ensure the necessary column is present
  if (!"dose" %in% colnames(data)) {
    stop("Data must contain a 'dose' column.")
  }

  # Sort and get unique dose levels
  dose_all <- sort(unique(data$dose))

  # Initialize vector to store mean response estimates
  mu_hat <- numeric(length(dose_all))

  # Loop over each dose level
  for (j in seq_along(dose_all)) {
    dat_tmp <- data
    dat_tmp$dose <- dose_all[j]

    # Predict the mean response for the j-th dose level
    mu_t <- predict(glm_obj, dat_tmp, type = "response")
    mu_hat[j] <- mean(mu_t)
  }

  return(mu_hat)
}
