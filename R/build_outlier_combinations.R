#'  Build outlier combinations
#'  
#'  Build all possible combinations of over- and under-expressed outliers.
#'  
#'  @details
#'  Generate all possible combinations of outliers, starting from the top and
#'  bottom of the sorted input.
#'  
#'  @param z a vector of sorted z-scores
#'  @seealso \code{\link{find_tissue_outliers}}
#'  @return a list of binary vectors. Each vector represents an outlier model 
#'  where outliers are denoted by 1, otherwise 0.
#'  @export
build_outlier_combinations <- function(z) {
  obs <- list()
  n <- length(z)
  mid <- ceiling(n/2)
  for (i in 1:mid) {
    obs[[i]] <- c(rep(1, i), rep(0, n-i))
  }
  for (i in (mid+1):n) {
    obs[[i]] <- c(rep(0, i-1), rep(1, n-i+1))
  }
  for (i in 1:(mid)) {
    for (j in (mid+1):n) {
      if (all(obs[[i]] + obs[[j]] == 1)) {
        next;
      } 
      obs[[length(obs)+1]] <- obs[[i]] + obs[[j]] 
    }
  }
  obs[[length(obs)+1]] <- rep(0, n)
  
  for (i in 1:length(obs)) {
    obs[[i]][is.na(z)] <- NA
  }
  
  return(unique(obs))
}