context("Build all outlier combinations")

z <- c(1, 5, 10, 20)

test_that("All possible models are returned for n = 4", {
  target <- build_outlier_combinations(z)
  expect_equal(length(target), 8)
})

test_that("All possible models are return for n = 3", {
  target <- build_outlier_combinations(z[1:3])
  expected <- list(c(1, 0, 0),
                   c(1, 1, 0),
                   c(0, 0, 1),
                   c(1, 0, 1),
                   c(0, 0, 0))
  expect_equal(length(target), length(expected))
})

