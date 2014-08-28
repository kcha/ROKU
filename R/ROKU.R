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

preprocess_with_tukey <- function(x) {
  # Step 1 of ROKU
  tu <- tukey.biweight(x[!is.na(x)])
  new.x <- x - tu
  return(new.x)
}

find_tissue_outliers <- function(x, epsilon=1e-04) {
  # Step 3 of ROKU - finding tissue outliers
  best.model <- list(Model=rep(0, length(x)), Outlier.Detection.Method=NA)
  best.model[[1]][is.na(x)] <- NA
    
  sd <- sd(x, na.rm=T)
  
#   if (attr(z, "scaled:scale") == 0) return(best.model)
  if (all(is.na(x)))
    return(best.model)
  if (is.na(sd)) 
    return(best.model)
  if (sd < epsilon)
    return(best.model)
  
  # Get Z-scores
  z <- scale(x)
#   z <- (x - median(x, na.rm=T))/(mad(x, na.rm=T) + epsilon)
  ord <- order(z)
  
  # For each combination of U:
  #   Find best combunation of U
  outlier.models <- build_outlier_combinations(z[ord])
  if (is.null(outlier.models)) return(best.model)
  
  U <- rep(NA, length(outlier.models))
  for (i in 1:length(U)) {
    U[i] <- U.statistic(z[ord], outlier.models[[i]])
  }

#   print(outlier.models)
#   print(U)
#   browser()
  
  min.U <- which.min(U)
  if (length(min.U) != 0) {
    s <- outlier.models[[min.U]]
    
    # 2 step approach:
    #   1) If an outlier is found, assign outlier labels to the original vector
    #   2) If not, proceed to doing k-means clustering (k=2)
    if (sum(s, na.rm=T) != 0) {
      best.model[[1]][ord[s == 1]] <- 1
      best.model[[2]] <- "AIC"
    } else {
      km <- kmeans(z[!is.na(z)], 2)
      if (diff(km$size) > 0.1*length(x)) {
        smallest.cluster <- which.min(km$size)
        best.model[[1]][!is.na(z)][km$cluster == smallest.cluster] <- 1
        best.model[[2]] <- "KM"
      }
    }
    
  } 
  return(best.model)
}

build_outlier_combinations <- function(z) {
  obs <- list()
  n <- length(z)
  mid <- floor(n/2)
  for (i in 1:mid) {
    obs[[i]] <- c(rep(1, i), rep(0, n-i))
  }
  for (i in (mid+1):n) {
    obs[[i]] <- c(rep(0, i-1), rep(1, n-i+1))
  }
  for (i in 1:(mid-1)) {
    for (j in (mid+2):n) {
      obs[[length(obs)+1]] <- obs[[i]] + obs[[j]] 
    }
  }
  obs[[length(obs)+1]] <- rep(0, n)
  
  for (i in 1:length(obs)) {
    obs[[i]][is.na(z)] <- NA
  }
  
#   ix <- which(abs(z) > c)
#   
#   if (length(ix) == 0) {
#     return(NULL)
#   }
#   
#   for (i in 1:length(ix)) {
#     comb <- combn(ix, i)
#     for (j in 1:ncol(comb)) {
#       values <- rep(0, length(z))
#       values[comb[,j]] <- 1
#       obs[[length(obs)+1]] <- values
#     }
#   }
  return(unique(obs))
}

U.statistic <- function(x, outliers, epsilon=1e-04) {
  stopifnot(length(outliers) == length(x))
  non.outliers <- outliers == 0
  sigma <- sd(x[non.outliers], na.rm=T)
  n <- sum(non.outliers, na.rm=T)
  s <- sum(outliers, na.rm=T)
  U <- n*log(sigma + epsilon) + sqrt(2)*s*(log(factorial(n))/n)
  return(U)
}

entropy.shannon <- function(x, unit = c("log", "log2", "log10")) {
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
