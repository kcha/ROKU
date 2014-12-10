#' ROKU post filtering
#' 
#' After running \code{\link{ROKU}}, filter out genes with outliers that have
#' Shannon Entropy greater than \code{min_diff}.
#' 
#' @param r List outputted by \code{\link{ROKU}}
#' @param min_diff Minimum entropy to filter outliers
#' @return List like \code{r}, but filtered.
#' @export
#' @examples
#' \dontrun{
#' rk <- ROKU(psidata)
#' candidates <- filter_outliers(rk, min_diff = 15)
#' }
filter_outliers <- function(r, min_diff=15) {
  ROKU.post_filter.index <- function(df, min_diff) {
    t <- which(rowSums(abs(df) > min_diff, na.rm=T) >= 1)
    return(t)
  }
  
  t <- ROKU.post_filter.index(r$Outliers, min_diff)
  post <- list(Entropy=r$Entropy[t],
               Entropy.Normalized=r$Entropy.Normalized[t],
               Outlier.Detection.Method=r$Outlier.Detection.Method[t],
               Outliers=r$Outliers[t,])
  return(post)
}