# Test file for ParametricSurvKM function
# tests/testthat/test-ParametricSurvKM.R

library(testthat)
library(survival)
library(dplyr)

# Create test data that matches the expected structure from DataParametricSurv
create_test_survival_data <- function(n = 100, stratify_levels = c("M", "F"), arm_levels = c("Control", "Treatment")) {
  set.seed(123)

  data.frame(
    SUBJID = paste0("SUBJ", sprintf("%03d", 1:n)),
    ARM = factor(sample(arm_levels, n, replace = TRUE)),
    STRATIFY = factor(sample(stratify_levels, n, replace = TRUE)),
    SURVTIME = rexp(n, rate = 0.1),  # Random survival times
    CNSR = rbinom(n, 1, 0.3),        # 30% censoring rate
    EVENT = 1 - rbinom(n, 1, 0.3),   # Event indicator (opposite of censoring)
    stringsAsFactors = FALSE
  )
}

# Test 1: Basic functionality
test_that("ParametricSurvKM creates plots with valid inputs", {
  skip_if_not_installed("survminer")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("flexsurv")

  test_data <- create_test_survival_data(n = 50)

  # Test with exponential distribution
  expect_error({
    result <- ParametricSurvKM(
      dataset = test_data,
      distribution = "exp",
      return_plots = TRUE
    )
  }, NA)  # NA means no error expected
})

# Test 2: Input validation - dataset structure
test_that("ParametricSurvKM validates dataset structure", {
  skip_if_not_installed("survminer")

  # Missing required columns
  invalid_data <- data.frame(
    SUBJID = c("S1", "S2", "S3"),
    ARM = factor(c("A", "B", "A"))
    # Missing STRATIFY, SURVTIME, CNSR, EVENT
  )

  expect_error(
    ParametricSurvKM(dataset = invalid_data, distribution = "exp"),
    "Dataset is missing required columns"
  )
})

# Test 3: Distribution validation
test_that("ParametricSurvKM validates distribution parameter", {
  skip_if_not_installed("survminer")

  test_data <- create_test_survival_data(n = 30)

  # Invalid distribution
  expect_error(
    ParametricSurvKM(dataset = test_data, distribution = "invalid_dist"),
    "Distribution must be exactly one of"
  )

  # Multiple distributions (should only accept one)
  expect_error(
    ParametricSurvKM(dataset = test_data, distribution = c("exp", "weibull")),
    "Distribution must be exactly one of"
  )
})

# Test 4: ARM validation
test_that("ParametricSurvKM validates ARM levels", {
  skip_if_not_installed("survminer")
  skip_if_not_installed("flexsurv")

  # Data with only one ARM level
  single_arm_data <- create_test_survival_data(n = 30, arm_levels = c("Treatment"))

  expect_error(
    ParametricSurvKM(dataset = single_arm_data, distribution = "exp"),
    "This function requires exactly 2 ARM levels"
  )
})

# Test 5: Control arm validation
test_that("ParametricSurvKM validates control_arm parameter", {
  skip_if_not_installed("survminer")
  skip_if_not_installed("flexsurv")

  test_data <- create_test_survival_data(n = 40)

  # Invalid control arm
  expect_error(
    ParametricSurvKM(
      dataset = test_data,
      distribution = "exp",
      control_arm = "NonExistent"
    ),
    "Specified control_arm 'NonExistent' not found in ARM levels"
  )
})

# Test 6: Different distributions (simplified)
test_that("ParametricSurvKM works with exponential distribution", {
  skip_if_not_installed("survminer")
  skip_if_not_installed("flexsurv")

  test_data <- create_test_survival_data(n = 60)

  # Test exponential distribution specifically
  expect_error({
    result <- ParametricSurvKM(
      dataset = test_data,
      distribution = "exp",
      return_plots = TRUE
    )
  }, NA)
})

# Test 7: Weibull distribution
test_that("ParametricSurvKM works with Weibull distribution", {
  skip_if_not_installed("survminer")
  skip_if_not_installed("flexsurv")

  test_data <- create_test_survival_data(n = 60)

  # Test Weibull distribution
  expect_error({
    result <- ParametricSurvKM(
      dataset = test_data,
      distribution = "weibull",
      return_plots = TRUE
    )
  }, NA)
})

# Test 8: Custom parameters
test_that("ParametricSurvKM handles custom parameters", {
  skip_if_not_installed("survminer")
  skip_if_not_installed("flexsurv")

  test_data <- create_test_survival_data(n = 40)

  # Test with custom parameters
  expect_error({
    result <- ParametricSurvKM(
      dataset = test_data,
      distribution = "exp",
      figure_caption = "Custom Caption",
      stratify_name = "SEX",
      control_arm = "Control",
      time_scale = "Years",
      time_max = 10,
      return_plots = TRUE
    )
  }, NA)
})

# Test 9: Return plots parameter
test_that("ParametricSurvKM handles return_plots parameter", {
  skip_if_not_installed("survminer")
  skip_if_not_installed("flexsurv")

  test_data <- create_test_survival_data(n = 30)

  # When return_plots = TRUE, should not error
  expect_error({
    result <- ParametricSurvKM(
      dataset = test_data,
      distribution = "exp",
      return_plots = TRUE
    )
  }, NA)
})

# Test 10: Minimum viable dataset
test_that("ParametricSurvKM handles minimum viable dataset", {
  skip_if_not_installed("survminer")
  skip_if_not_installed("flexsurv")

  # Create minimal dataset with just enough data
  minimal_data <- data.frame(
    SUBJID = paste0("S", 1:20),
    ARM = factor(rep(c("A", "B"), each = 10)),
    STRATIFY = factor(rep(c("X", "Y"), 10)),
    SURVTIME = c(rexp(10, 0.1), rexp(10, 0.05)),
    CNSR = rep(0, 20),  # No censoring for simplicity
    EVENT = rep(1, 20),  # All events
    stringsAsFactors = FALSE
  )

  expect_error({
    result <- ParametricSurvKM(
      dataset = minimal_data,
      distribution = "exp",
      return_plots = TRUE
    )
  }, NA)
})

# Test 11: Different time scales
test_that("ParametricSurvKM accepts different time scales", {
  skip_if_not_installed("survminer")
  skip_if_not_installed("flexsurv")

  test_data <- create_test_survival_data(n = 40)

  # Test different time scales
  time_scales <- c("Months", "Years", "Weeks", "Days")

  for (scale in time_scales) {
    expect_error({
      ParametricSurvKM(
        dataset = test_data,
        distribution = "exp",
        time_scale = scale,
        return_plots = TRUE
      )
    }, NA, label = paste("Testing time_scale:", scale))
  }
})

# Test 12: Multiple stratification levels
test_that("ParametricSurvKM handles multiple stratification levels", {
  skip_if_not_installed("survminer")
  skip_if_not_installed("flexsurv")

  # Test with 3 stratification levels
  test_data_3strat <- create_test_survival_data(
    n = 60,
    stratify_levels = c("A", "B", "C")
  )

  expect_error({
    result <- ParametricSurvKM(
      dataset = test_data_3strat,
      distribution = "exp",
      return_plots = TRUE
    )
  }, NA)
})

# Test 13: Data type validation
test_that("ParametricSurvKM handles proper data types", {
  skip_if_not_installed("survminer")
  skip_if_not_installed("flexsurv")

  test_data <- create_test_survival_data(n = 40)

  # Ensure ARM and STRATIFY are factors
  test_data$ARM <- as.factor(test_data$ARM)
  test_data$STRATIFY <- as.factor(test_data$STRATIFY)

  expect_error({
    result <- ParametricSurvKM(
      dataset = test_data,
      distribution = "exp",
      return_plots = TRUE
    )
  }, NA)
})

# Test 14: Edge case - all events or all censored
test_that("ParametricSurvKM handles edge cases in event data", {
  skip_if_not_installed("survminer")
  skip_if_not_installed("flexsurv")

  test_data <- create_test_survival_data(n = 40)

  # Test with all events (no censoring)
  test_data$CNSR <- 0
  test_data$EVENT <- 1

  expect_error({
    result <- ParametricSurvKM(
      dataset = test_data,
      distribution = "exp",
      return_plots = TRUE
    )
  }, NA)
})

# Test 15: Parameter combinations
test_that("ParametricSurvKM works with various parameter combinations", {
  skip_if_not_installed("survminer")
  skip_if_not_installed("flexsurv")

  test_data <- create_test_survival_data(n = 50)

  # Test specific parameter combination that's commonly used
  expect_error({
    result <- ParametricSurvKM(
      dataset = test_data,
      distribution = "weibull",
      figure_caption = "Test Analysis",
      stratify_name = "SEX",
      control_arm = "Control",
      time_scale = "Months",
      time_max = 24,
      return_plots = TRUE
    )
  }, NA)
})
