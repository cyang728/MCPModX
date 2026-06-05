## ---------------------------------------------------------------------------
## Candidate dose-response shapes for the MCP and Mod steps.
## A "model set" is a named list whose names are shapes and whose values are the
## guesstimate parameter(s) of that shape, in the DoseFinding convention:
##   linear    = NULL                 f(d) = d
##   emax      = ED50 (>0), length>=1 f(d) = d / (ED50 + d)
##   sigEmax   = c(ED50, h)           f(d) = d^h / (ED50^h + d^h)
##   quadratic = delta (<0)           f(d) = d + delta d^2
##   exponential = tau (>0)           f(d) = exp(d/tau) - 1
## Each candidate is reduced to a *placebo-adjusted, unit-norm* shape on the
## active doses, which is all the optimal-contrast test and the GLS fit need.
## ---------------------------------------------------------------------------

#' Standardised candidate shapes on the active doses
#'
#' @param models named list of candidate shapes and guesstimates (see Details).
#' @param doses numeric vector of doses including placebo (the minimum).
#' @return a matrix with one column per candidate, rows = active doses, each
#'   column the placebo-adjusted shape scaled to unit Euclidean norm.
#' @export
build_models <- function(models, doses) {
  doses <- sort(unique(doses))
  d0 <- doses[1]; da <- doses[-1]
  shape <- function(name, par) {
    f <- switch(name,
      linear      = function(x) x,
      emax        = function(x) x / (par[1] + x),
      sigemax     = ,
      sigEmax     = function(x) x^par[2] / (par[1]^par[2] + x^par[2]),
      quadratic   = function(x) x + par[1] * x^2,
      exponential = function(x) exp(x / par[1]) - 1,
      stop(sprintf("unknown model shape '%s'", name)))
    f(da) - f(d0)
  }
  cols <- list(); nm <- character(0)
  for (i in seq_along(models)) {
    name <- names(models)[i]; par <- models[[i]]
    if (name %in% c("emax", "exponential") && length(par) > 1) {
      for (p in par) { cols[[length(cols) + 1]] <- shape(name, p)
        nm <- c(nm, sprintf("%s_%g", name, p)) }
    } else {
      cols[[length(cols) + 1]] <- shape(name, par); nm <- c(nm, name)
    }
  }
  M <- do.call(cbind, cols); colnames(M) <- nm
  apply(M, 2, function(v) v / sqrt(sum(v^2)))    # scale is irrelevant for contrasts
}

## Evaluate a fitted shape (for the Mod step) on the active doses.
.shape_pred <- function(name, par, da, d0) {
  f <- switch(name,
    linear    = function(x) x,
    emax      = function(x) x / (par[1] + x),
    quadratic = function(x) x + par[1] * x^2,
    stop("unsupported Mod shape"))
  f(da) - f(d0)
}
