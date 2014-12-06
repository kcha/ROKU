#' ROKU post filtering
#' 
#' @export
filter_outliers <- function(r, ...) {
  ROKU.post_filter.index <- function(df, minDiff=15) {
    t <- which(rowSums(df > minDiff, na.rm=T) >= 1)
    return(t)
  }
  
  t <- ROKU.post_filter.index(r$outliers, ...)
  post <- list(entropy=r$entropy[t],
               entropy.pct=r$entropy.pct[t],
               outliers=r$outliers[t,])
  return(post)
}