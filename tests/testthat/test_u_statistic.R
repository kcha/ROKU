context("U statistic calculation")

test_that("U statistic calculated for n = 3", {
  x <- c(1, 2, 10)
  target <- u_statistic(x, c(0, 0, 1))
  expected <- -0.2027353
})