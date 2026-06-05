#' MCPModX: MCP-Mod with model-assisted covariate adjustment under missing outcomes
#'
#' Extends the Multiple Comparison Procedure and Modelling (MCP-Mod) framework
#' for dose-finding with model-assisted covariate adjustment and handling of
#' outcomes that are missing at random, via inverse-probability-weighted
#' g-computation / AIPW with stacked-sandwich (M-estimation) inference.
#'
#' @keywords internal
#' @importFrom stats as.formula binomial cov2cor delete.response fitted gaussian
#'   glm model.matrix optimize pnorm predict qnorm setNames terms
#' @importFrom mvtnorm qmvnorm pmvnorm
"_PACKAGE"
