# tests/testthat/test-DataParametricSurv.R

test_that("DataParametricSurv validates parameters correctly", {
  # Test invalid population
  expect_error(
    DataParametricSurv(
      adslph_path = "dummy.sas7bdat",
      adtte_path = "dummy.sas7bdat",
      population = "INVALID"
    ),
    "Population must be one of"
  )
})

test_that("DataParametricSurv handles stratification correctly", {
  # Test with NULL stratification
  # This would need actual test data to run properly
  skip("Requires test data files")
})
