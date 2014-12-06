entropy_shannon <- function(x, unit = c("log", "log2", "log10")) {
  epsilon <- 1e-08
  unit <- match.arg(unit)
  p <- x / sum(x)
  p <- p + epsilon
  H <- -sum(p * log(p))
  if (unit == "log2")
    H <- H/log(2)
  if (unit == "log10")
    H <- H/log(10)
  return(H)
}