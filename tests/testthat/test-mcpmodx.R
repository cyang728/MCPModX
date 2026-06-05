## End-to-end sanity tests for the MCP-ModX estimator.

make_data <- function(n = 120, doses = c(0, 0.05, 0.2, 0.4, 0.7, 1),
                      beta = 0.8, miss = 0, seed = 1) {
  set.seed(seed)
  N <- n * length(doses)
  dose <- rep(doses, each = n)
  x1 <- rnorm(N); x2 <- rbinom(N, 1, 0.5); x3 <- runif(N)
  lp <- -1.2 * x1 + 0.6 * x2 - x3 + beta * dose
  y <- rbinom(N, 1, plogis(lp))
  d <- data.frame(dose, x1, x2, x3, y)
  if (miss > 0) {                                     # MAR on x1, x2
    pr <- plogis(2 - 0.7 * x1 + 0.3 * x2)             # ~15-20% missing
    R <- rbinom(N, 1, pr); d$y[R == 0] <- NA
  }
  d
}

MODS <- list(linear = NULL, emax = c(0.05, 0.2, 0.5), quadratic = -0.85)

test_that("complete-data adjusted fit runs and is well-formed", {
  d <- make_data(beta = 0.8)
  fit <- MCPModX(d, "y", "dose", models = MODS, covariates = c("x1", "x2", "x3"))
  expect_s3_class(fit, "MCPModX")
  expect_length(fit$delta, 5)
  expect_true(all(is.finite(fit$delta)))
  expect_true(isSymmetric(unname(fit$Sigma_delta), tol = 1e-8))
  expect_true(all(eigen(fit$Sigma_delta, only.values = TRUE)$values > -1e-8))
  expect_true(fit$mcp$reject %in% c(0L, 1L))
})

test_that("a real dose effect is detected with high probability", {
  rej <- mean(replicate(40, {
    d <- make_data(beta = 1.2, seed = sample.int(1e6, 1))
    MCPModX(d, "y", "dose", models = MODS, covariates = c("x1", "x2", "x3"))$mcp$reject
  }))
  expect_gt(rej, 0.7)
})

test_that("flat truth controls type-I error near nominal", {
  rej <- mean(replicate(60, {
    d <- make_data(beta = 0, seed = sample.int(1e6, 1))
    MCPModX(d, "y", "dose", models = MODS, covariates = c("x1", "x2", "x3"))$mcp$reject
  }))
  expect_lt(rej, 0.18)                                # generous bound for 60 reps
})

test_that("missing outcomes are handled via IPW (sandwich runs)", {
  d <- make_data(beta = 0.8, miss = 0.2)
  fit <- MCPModX(d, "y", "dose", models = MODS,
                 covariates = c("x1", "x2", "x3"), obs_covariates = c("x1", "x2"))
  expect_true(fit$n_missing > 0)
  expect_true(all(is.finite(fit$delta)))
  expect_true(all(eigen(fit$Sigma_delta, only.values = TRUE)$values > -1e-8))
})

test_that("unadjusted analysis works (covariates = NULL)", {
  d <- make_data(beta = 0.8)
  fit <- MCPModX(d, "y", "dose", models = MODS, covariates = NULL)
  expect_length(fit$delta, 5)
})

test_that("gaussian (continuous) outcome is supported", {
  set.seed(3)
  doses <- c(0, 0.25, 0.5, 0.75, 1); n <- 100
  dose <- rep(doses, each = n); x1 <- rnorm(n * length(doses))
  y <- 1 + 2 * x1 + 1.5 * dose + rnorm(n * length(doses))
  d <- data.frame(dose, x1, y)
  fit <- MCPModX(d, "y", "dose", models = list(linear = NULL, emax = 0.2),
                 covariates = "x1", family = gaussian())
  expect_equal(unname(fit$family), "gaussian")
  expect_true(all(is.finite(fit$delta)))
})
