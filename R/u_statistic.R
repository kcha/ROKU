#' U statistic
#' 
#' Calculate U statistic based on method described by 
#' \href{http://www.ncbi.nlm.nih.gov/pubmed/12499447}{Kadota et al. (2003)}
#' 
#' @details
#' Formula:
#' 
#' \deqn{U = n\log{\sigma} + \sqrt{2} \times s \times \frac{\log{n!}}{n}}
#' 
#' where \emph{n} and \emph{s} denote the numbers of non-outlier and outlier candidates, 
#' and \eqn{\sigma} denotes the standard deviation of the observations of the 
#' \emph{n} non-outlier candidates.
#' 
#' @param x Vector of expression values
#' @param outliers Binary vector describing potential outlier model
#' @param epsilon Offset for log-transform
#' @references
#' Kadota et al (2006). ROKU: a novel method for identification of tissue-specific genes. BMC Genomics, 7:294. 
#' \url{http://www.biomedcentral.com/1471-2105/7/294}
#' 
#' Kadota et al (2003). Detection of genes with tissue-specific expression patterns using Akaike's information criterion procedure. Physiol Genomics, 12(3):251-9. \url{http://www.ncbi.nlm.nih.gov/pubmed/12499447}
#' @export
#' @return U statistic (numeric)
#' @seealso \code{\link{find_tissue_outliers}}
u_statistic <- function(x, outliers, epsilon=1e-04) {
  stopifnot(length(outliers) == length(x))
  non.outliers <- outliers == 0
  sigma <- sd(x[non.outliers], na.rm=T)
  n <- sum(non.outliers, na.rm=T)
  s <- sum(outliers, na.rm=T)
  U <- n*log(sigma + epsilon) + sqrt(2)*s*(log(factorial(n))/n)
  return(U)
}