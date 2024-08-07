#test_that("multiplication works", {
#  expect_equal(2 * 2, 4)
#})

test_that("MCPModX runs without errors for standard input", {
  # Simulate a dataset
  set.seed(123)
  dat <- data.frame(
    dose = rep(c(0, 0.05, 0.2, 0.5), each = 25),
    y = rbinom(100, 1, prob = 0.5),
    x_star = rnorm(100)
  )

  # Run the function
  result <- MCPModX(data = dat)

  # Check that the result is a list
  expect_type(result, "list")

  # Check that the list contains expected elements
  expect_true("test_significant" %in% names(result))
  expect_true("weighted_estimand" %in% names(result))
})

test_that("MCPModX throws an error if 'dose' or 'y' columns are missing", {
  # Simulate a dataset without 'dose'
  dat_no_dose <- data.frame(
    y = rbinom(100, 1, prob = 0.5),
    x_star = rnorm(100)
  )

  # Expect an error
  expect_error(MCPModX(data = dat_no_dose), "Data must contain 'dose' and 'y' columns.")

  # Simulate a dataset without 'y'
  dat_no_y <- data.frame(
    dose = rep(c(0, 0.05, 0.2, 0.5), each = 25),
    x_star = rnorm(100)
  )

  # Expect an error
  expect_error(MCPModX(data = dat_no_y), "Data must contain 'dose' and 'y' columns.")
})

test_that("MCPModX throws an error if placebo group is missing", {
  # Simulate a dataset without a placebo group
  dat_no_placebo <- data.frame(
    dose = rep(c(0.05, 0.2, 0.5), each = 25),
    y = rbinom(75, 1, prob = 0.5),
    x_star = rnorm(75)
  )

  # Expect an error
  expect_error(MCPModX(data = dat_no_placebo), "Data must include a placebo group with dose set to 0.")
})

test_that("MCPModX handles edge cases for the bnds parameter", {
  # Simulate a dataset
  dat <- data.frame(
    dose = rep(c(0, 0.05, 0.2, 0.5), each = 25),
    y = rbinom(100, 1, prob = 0.5),
    x_star = rnorm(100)
  )

  # Test with different bounds
  result <- MCPModX(data = dat, bnds = c(0, 1))
  expect_type(result, "list")

  result <- MCPModX(data = dat, bnds = c(1, 2))
  expect_type(result, "list")
})

test_that("MCPModX handles different estimand values correctly", {
  # Simulate a dataset
  dat <- data.frame(
    dose = rep(c(0, 0.05, 0.2, 0.5), each = 25),
    y = rbinom(100, 1, prob = 0.5),
    x_star = rnorm(100)
  )

  # Test for logOR
  result_logOR <- MCPModX(data = dat, estimand = "logOR")
  expect_type(result_logOR, "list")

  # Test for RD
  result_RD <- MCPModX(data = dat, estimand = "RD")
  expect_type(result_RD, "list")

  # Test for logRR
  result_logRR <- MCPModX(data = dat, estimand = "logRR")
  expect_type(result_logRR, "list")
})
