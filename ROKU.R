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

find_tissue_outliers <- function(x, epsilon=1e-04) {
  # Step 3 of ROKU - finding tissue outliers
  
  # Get Z-scores
  z <- scale(x)
#   z <- (x - median(x))/(mad(x) + epsilon)
  
  # Sort Z-scores
#   st <- order(z)
  if (attr(z, "scaled:scale") == 0) {
    return(rep(0, length(z)))
  }
  
  # For each combination of U:
  #   Find best combunation of U
  outlier.models <- build_outlier_combinations(z)
  U <- rep(NA, length(outlier.models))
  for (i in 1:length(U)) {
    U[i] <- U.statistic(z, outlier.models[[i]])
  }
  print(outlier.models)
  print(U)
  best.model <- outlier.models[[which.min(U)]]
#   outliers <- rep(0,length(z))
#   outliers[st[best.model == 1]] <- 1
  
  return(best.model)
}

build_outlier_combinations <- function(z, c=0.5) {
  obs <- list()
  
  ix <- which(abs(z) > c)
  
  for (i in 1:length(ix)) {
    comb <- combn(ix, i)
    for (j in 1:ncol(comb)) {
      values <- rep(0, length(z))
      values[comb[,j]] <- 1
      obs[[length(obs)+1]] <- values
    }
  }
  return(obs)
}

U.statistic <- function(x, outliers, epsilon=1e-04) {
  stopifnot(length(outliers) == length(x))
  non.outliers <- outliers == 0
  sigma <- sd(x[non.outliers])
  n <- sum(non.outliers)
  s <- sum(outliers)
  U <- n*log(sigma + epsilon) + sqrt(2)*s*(log(factorial(n))/n)
  return(U)
}

ROKU <- function(m, cores=8) {
  if (nrow(m) < 10) cores <- 2
  
  # pre-process
  x.prime <- mclapply(1:nrow(m), function(i) preprocess_with_tukey(m[i,]), 
                      mc.cores=cores)
  x.prime <- do.call("rbind", x.prime)
  dimnames(x.prime) <- list(rownames(m), colnames(m))
  
  # calculate Shannon entropy
  H <- mclapply(1:nrow(x.prime), function(i) entropy(x.prime[i,], unit="log2"),
                mc.cores=cores)
  H <- unlist(H)
  names(H) <- rownames(m)
  
  # find outliers
  O <- mclapply(1:nrow(m), function(i) find_tissue_outliers(m[i,]), mc.cores=cores)
  O <- do.call("rbind", O)
  dimnames(O) <- list(rownames(m), colnames(m))
  
  return(list(entropy=H, tissue.specific.outliers=O))
}


