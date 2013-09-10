# 2013-09-10
#
# Implementation of "ROKU: a novel method for identification of tissue-specific genes"
# from Kadota et al., BMC Bioinformatics, 2006
#
# General steps:
# 1) Process each vector: x -> x' using Tukey Biweight
# 2) Calculate entropy: H(x')
# 3) Assign tissues detected as outliers using AIC
library(entropy)
library(affy)
library(parallel)


preprocess_with_tukey <- function(x) {
  # Step 1 of ROKU
  tu <- tukey.biweight(x)
  new.x <- abs(x - tu)
  return(new.x)
}

find_tissue_outliers <- function(x) {
  # Step 3 of ROKU - finding tissue outliers
  
  # Get Z-scores
  z <- scale(x)
  
  # Sort Z-scores
  z.sorted <- sort(z)
  
  # For each combination of U:
  #   Find best combunation of U
  U <- rep(NA, length(z.sorted))
  outlier.models <- build_outlier_combinations(length(z.sorted))
  for (i in 1:length(U)) {
    U[i] <- U.statistic(z.sorted, outlier.models[[i]])
  }
  print(outlier.models)
  print(U)
  return(outlier.models[[which.min(U)]])
}

build_outlier_combinations <- function(n) {
  obs <- list()
  mid <- floor(n/2)
  for (i in 1:mid) {
    obs[[i]] <- c(rep(1, i), rep(0, n-i))
  }
  
  for (i in (mid+1):n) {
    obs[[i]] <- c(rep(0, i-1), rep(1, n-i+1))
  }
  return(obs)
}

U.statistic <- function(x, outliers) {
  stopifnot(length(outliers) == length(x))
  non.outliers <- outliers == 0
  sigma <- sd(x[non.outliers])
  n <- sum(non.outliers)
  s <- sum(outliers)
  U <- n*log(sigma) + sqrt(2)*s*(log(factorial(n))/n)
  return(U)
}