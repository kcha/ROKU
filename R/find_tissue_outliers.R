#' Find over- or under-expressed tissue-specific genes
#' 
#' The third step in ROKU algorithm. Identify tissue-specific outliers.
#' 
#' @details
#' See the \href{http://www.biomedcentral.com/1471-2105/7/294#sec4}{Methods} from 
#' the manuscript for more details.
#' 
#' Identifiers over- or under-expressed tissue-specific genes using
#' the U statistic from 
#' \href{http://www.ncbi.nlm.nih.gov/pubmed/12499447}{Kadota et al. (2003)}. See
#' \code{\link{u_statistic}} for more details.
#'
#' The combination of outliers that where \emph{U} is minimized is returned.
#' 
#' If no outliers are found, then k-means clustering is performed to cluster
#' the expression values.
#' 
#' @param x a vector of expression values
#' @param epsilon a very small number for offsetting zero values
#' @return a binary vector representing the ouliers (1) and non-outliers (0)
#' @references
#' Kadota et al (2006). ROKU: a novel method for identification of tissue-specific genes. BMC Genomics, 7:294. 
#' \url{http://www.biomedcentral.com/1471-2105/7/294}
#' 
#' Kadota et al (2003). Detection of genes with tissue-specific expression patterns using Akaike's information criterion procedure. Physiol Genomics, 12(3):251-9. \url{http://www.ncbi.nlm.nih.gov/pubmed/12499447}
#' @seealso
#' \code{\link{u_statistic}} for calculation of \emph{U}
#' \code{\link{build_outlier_combinations}} for building all possible outlier models
#' @export
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
    U[i] <- u_statistic(z[ord], outlier.models[[i]])
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