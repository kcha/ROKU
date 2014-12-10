context("Transform vector by Tukey biweight")

test_that("Vector subtracted by Tukey biweight", {
  x <- c(1, 5, 10)
  tu <- affy::tukey.biweight(x)
  expected <- x - tu
  target <- preprocess_with_tukey(x)
  expect_equal(target, expected)
})

test_that("Vector with NA subtracted by Tukey biweight", {
  x <- c(1, 5, NA, 10)
  tu <- affy::tukey.biweight(x[!is.na(x)])
  expected <- x - tu
  target <- preprocess_with_tukey(x)
  expect_equal(target, expected)
})