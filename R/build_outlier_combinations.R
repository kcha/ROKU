
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