u_statistic <- function(x, outliers, epsilon=1e-04) {
  stopifnot(length(outliers) == length(x))
  non.outliers <- outliers == 0
  sigma <- sd(x[non.outliers], na.rm=T)
  n <- sum(non.outliers, na.rm=T)
  s <- sum(outliers, na.rm=T)
  U <- n*log(sigma + epsilon) + sqrt(2)*s*(log(factorial(n))/n)
  return(U)
}