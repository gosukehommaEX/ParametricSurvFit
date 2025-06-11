# tests/testthat/test-ExtractParams.R

test_that("ExtractParams validates distribution parameter", {
  expect_error(
    ExtractParams(NULL, "invalid_dist"),
    "Distribution must be one of"
  )

  expect_error(
    ExtractParams(NULL, "weibul"),  # typo in distribution name
    "Distribution must be one of"
  )
})

test_that("ExtractParams returns correct structure", {
  skip("Requires fitted model objects")
  # This test would verify the structure of returned lists
})
