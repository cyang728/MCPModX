#' MCPModX: Modified MCPMod for Binary Outcomes with Covariates
#'
#' This function performs MCPMod analysis with model-assisted variance for generalized outcomes such as binary and count.
#'
#' @param data A data frame containing the variables `dose`, `y`, and any covariates specified.
#' @param family A family for generalized linear model with different outcomes
#' @param covariates A character vector of covariate names to include in the model. Default is "x_star".
#' @param gauss_models A list of models for the dose-response relationship. Default includes linear, emax, and quadratic models.
#' @param estimand A character string specifying the type of measure to estimate.
#'        Options are "logOR" (log Odds Ratio), "RD" (Risk Difference), and "logRR" (log Risk Ratio).
#'        Default is "logOR".
#' @param alpha Significance level for the test. Default is 0.05.
#' @param bnds A numeric vector specifying the bounds for the optimization. Default is c(0, 0.2).
#'
#' @return A list containing the significant tests and the generalized weighted AIC fitted values.
#' @export
#'
#' @examples
#' # Example usage:
#' # Assuming `dat` is your dataset with `dose`, `y`, and covariates
#' # result <- MCPModX(data = dat, covariates = "x_star")
MCPModX <- function(data,
                    covariates = "x_star",
                    family = binomial(link = "logit"),
                    gauss_models = Mods(linear = NULL, emax = c(0.05, 0.20, 0.50), quadratic = -0.85,
                                        doses = unique(data$dose)),
                    estimand = c("logOR", "RD", "logRR"),
                    alpha = 0.05,
                    bnds = c(0, 2)) {

  # Ensure data contains 'dose' and 'y'
  if (!all(c("dose", "y") %in% colnames(data))) {
    stop("Data must contain 'dose' and 'y' columns.")
  }

  # Extract unique doses
  doses <- unique(sort(data$dose))

  # Check that placebo (dose = 0) is included
  if (min(doses) != 0) {
    stop("Data must include a placebo group with dose set to 0.")
  }

  # Ensure that the estimand is one of the allowed values
  estimand <- match.arg(estimand)

  # Construct formula based on whether covariates are used
  if(covariates == 1 || covariates == "1") covariates = NULL

  if (!is.null(covariates)) {
    formula_str <- paste("y ~ factor(dose) +", paste(covariates, collapse = "+"))
  } else {
    formula_str <- "y ~ factor(dose)"
  }
  formula <- as.formula(formula_str)

  # Fit the generalized linear model
  anovaMod <- glm(formula = formula, data = data, family = family)

  # Estimate the response means
  muhat <- calculate_gc_muhat(data = data, glm_obj = anovaMod)

  # Compute log odds ratios
  thetahat <- sapply(2:length(doses), function(i) {
    switch(estimand,
           logOR = calculate_logOR(muhat[i], muhat[1]),
           logRR = calculate_logRR(muhat[i], muhat[1]),
           RD = calculate_RD(muhat[i], muhat[1]))
  })

  # Calculate partial derivatives and variance
  vrhat <- calculate_gc_variance(data = data, glm_obj = anovaMod)
  S <- calculate_variance_estimate(mu = muhat, vr = vrhat, estimate = estimand)

  # Adjust variance-covariance matrix if no covariates are used
  if (is.null(covariates)) {
    S <- vcov(anovaMod)
    S <- S[-1, -1]
  }

  # Dose levels excluding placebo for power calculations
  doses_power <- doses[-1]

  # Check if the determinant of S is near zero and adjust accordingly
  if (det(S) <= 1e-100) {
    stop("determinant of adjusted covariance is singular.")
  } else {
    contMat <- optContr(gauss_models, doses_power, S = S, placAdj = TRUE)$contMat
  }
  rownames(contMat) <- doses_power

  # Compute test statistics
  ct <- as.vector(thetahat %*% contMat)
  covMat <- t(contMat) %*% S %*% contMat
  corMat <- cov2cor(covMat)
  den <- sqrt(diag(covMat))
  tStat <- ct / den
  corMat[is.nan(corMat)] <- 0

  # Calculate critical value using multivariate t-distribution
  ctrl <- do.call("mvtnorm.control", mvtnorm.control())
  crtl_val <- qmvt(1 - alpha,
                   tail = "lower.tail",
                   corr = corMat,
                   df = anovaMod$df.residual,
                   algorithm = ctrl,
                   interval = ctrl$interval)

  # Determine significant tests
  MCT_table = matrix(NA, nrow = length(tStat), ncol = 2)
  rownames(MCT_table) = names(tStat)
  colnames(MCT_table) = c("t-Stat", "p-value")
  MCT_table[,1] = tStat
  MCT_table[,2] = MCTpval(contMat,
                          corMat,
                          df = anovaMod$df.residual,
                          tStat,
                          alternative = "one.sided",
                          control = mvtnorm.control())

  # Fit models and calculate generalized weighted AIC
  fitting_models <- lapply(seq_along(gauss_models), function(i_model) {
    fitMod(doses[-1],
           thetahat,
           S = S,
           model = names(gauss_models)[i_model],
           type = "general",
           placAdj = TRUE,
           bnds = bnds)
  })

  aics = sapply(fitting_models, gAIC)
  min_aic = min(aics)
  aics = aics - min_aic

  fitted_values = sapply(fitting_models, function(model) predict(model, se.fit = FALSE, doseSeq = doses[-1], predType = "effect-curve"))
  gfitted_aic = rowSums(sapply(1:length(fitting_models), function(i) exp(-aics[i] / 2) * fitted_values[, i])) / sum(exp(-aics / 2))

  return(list(MCT_table = MCT_table,
              mu_hat = muhat,
              weighted_estimand = gfitted_aic,
              estimand = estimand))
}
