# tests/testthat/test-FitSurvMods.R

test_that("FitSurvMods validates input dataset", {
  # Test with missing required columns
  bad_dataset <- data.frame(SUBJID = 1:10, ARM = rep("A", 10))

  expect_error(
    FitSurvMods(bad_dataset),
    "Dataset is missing required columns"
  )
})

test_that("FitSurvMods validates distributions", {
  # Create minimal valid dataset structure
  valid_dataset <- data.frame(
    SUBJID = 1:10,
    ARM = factor(rep(c("A", "B"), 5)),
    STRATIFY = factor(rep("Overall", 10)),
    SURVTIME = runif(10, 1, 100),
    CNSR = rbinom(10, 1, 0.3),
    EVENT = 1 - rbinom(10, 1, 0.3)
  )

  expect_error(
    FitSurvMods(valid_dataset, distributions = c("invalid_dist")),
    "Invalid distributions specified"
  )

  expect_error(
    FitSurvMods(valid_dataset, distributions = c("exp", "invalid")),
    "Invalid distributions specified"
  )
})

test_that("FitSurvMods returns correct output format", {
  skip("Requires more comprehensive test data setup")
})
