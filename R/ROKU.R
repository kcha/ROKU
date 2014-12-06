# Copyright (c) 2013, Kevin Ha
# University of Toronto
#
# Implementation of "ROKU: a novel method for identification of tissue-specific genes"
# from Kadota et al., BMC Bioinformatics, 2006
# http://www.biomedcentral.com/1471-2105/7/294
#
# General steps:
# 1) Process each vector: x -> x' using Tukey Biweight
# 2) Calculate entropy: H(x')
# 3) Assign tissues detected as outliers using AIC
library(entropy)
library(affy)
library(parallel)


ROKU <- function(m, cores=8) {
  if (nrow(m) < 10) cores <- 2
  
  # pre-process
  m <- as.matrix(m)
  x.prime <- mclapply(1:nrow(m), function(i) preprocess_with_tukey(m[i,]), 
                      mc.cores=cores)
  
  # calculate Shannon entropy
  H <- mclapply(x.prime, function(i) entropy(abs(i[!is.na(i)]), unit="log2"),
                mc.cores=cores)
  H <- unlist(H)
  names(H) <- rownames(m)
  H.pct <- H / log2(rowSums(!is.na(m)))
 
  # find outliers
  O <- mclapply(x.prime, function(i) find_tissue_outliers(abs(i)), mc.cores=cores)
  M <- do.call("rbind", lapply(O, "[[", 1))
  M2 <- do.call("rbind", lapply(O, "[[", 2))
  dimnames(M) <- list(rownames(m), colnames(m))
  
  x.prime <- do.call("rbind", x.prime)
  dimnames(x.prime) <- list(rownames(m), colnames(m))

  v <- M*x.prime
  v <- data.frame(v, Outlier.Detection.Method=M2)
  
  return(list(entropy=H, entropy.pct=H.pct, outliers=v))
}

ROKU.2 <- function(m, cores) {
#   m <- as.matrix(m)
  
  R <- mclapply(1:nrow(m), function(i) ROKU.run(m[i,]), mc.cores=cores)
  R <- do.call("rbind", R)
#   dimnames(R) <- list(rownames(m), 
#                       c("Entropy", "Entropy.Normalized", colnames(m)))
  
  return(R)
}

ROKU.run <- function(x) {
#   x <- t(x)
  x.prime <- preprocess_with_tukey(x)
  H <- entropy(x.prime[!is.na(x.prime)], unit="log2")
  O <- find_tissue_outliers(t(x))
  values <- cbind(Entropy=H, Entropy.Normalized=H/log2(sum(!is.na(x))), O*x.prime)
  return(values)
}

ROKU.post_filter.index <- function(df, minDiff=15) {
  t <- which(rowSums(df > minDiff, na.rm=T) >= 1)
  return(t)
}

ROKU.post_filter <- function(r, ...) {
  t <- ROKU.post_filter.index(r$outliers, ...)
  post <- list(entropy=r$entropy[t],
               entropy.pct=r$entropy.pct[t],
               outliers=r$outliers[t,])
  return(post)
}

ROKU.as.df <- function(r, meta=NULL) {
  df <- data.frame(Entropy=r$entropy,
                   Entropy.Normalized=r$entropy.pct,
                   r$outliers)
  if (!is.null(meta)) {
    if (nrow(meta)!=nrow(r$outliers)) {
      stop("Number of rows don't match")
    }
    df <- data.frame(meta, df)  
  }
  return(df)
}
