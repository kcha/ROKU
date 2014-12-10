#' ROKU
#' 
#' A method for identification of tissue-specific genes
#' 
#' @details
#' See the \href{http://www.biomedcentral.com/1471-2105/7/294#sec4}{Methods} from the manuscript for more details.
#' 
#' General steps:
#' \enumerate{
#'  \item Process each vector: x -> x' using Tukey biweight
#'  \item Calculate entropy: H(x')
#'  \item Assign tissues detected as outliers using AIC
#'  }
#' @param m a \emph{n} by \emph{m} data frame where \emph{n} is the number of genes
#' and \emph{m} is the number of samples/tissues
#' @param cores number of cores to use for parallel processing. Default is 1. 
#' @return List:
#' \itemize{
#'  \item{Entropy}{Vector of Shannon entropy for each gene} 
#'  \item{Entropy.Normalized}{Vector of Shannon entropy normalized by all entropies}
#'  \item{Outlier.Detection.Method}{Method used to detect outliers: AIC, or KM}
#'  \item{Outliers}{Data frame with dimensions similar to \code{m}. Non-zero values
#'  indicate predicted outliers.}
#' }
#' @importFrom parallel mclapply
#' @import entropy
#' @export
#' @references
#' Kadota et al (2006). ROKU: a novel method for identification of tissue-specific genes. BMC Genomics, 7:294. 
#' \url{http://www.biomedcentral.com/1471-2105/7/294}
#' @examples
#' rk <- ROKU(psidata)
ROKU <- function(m, cores=1) {
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
  O <- mclapply(1:nrow(m), function(i) find_tissue_outliers(m[i,]), mc.cores=cores)
  M <- do.call("rbind", lapply(O, "[[", 1))
  M2 <- sapply(O, "[[", 2)
  names(M2) <- rownames(m)
  dimnames(M) <- list(rownames(m), colnames(m))
  
  x.prime <- do.call("rbind", x.prime)
  dimnames(x.prime) <- list(rownames(m), colnames(m))

  v <- M*x.prime
  v <- data.frame(v)
  
  return(list(Entropy=H, Entropy.Normalized=H.pct, 
              Outlier.Detection.Method=M2,
              Outliers=v))
}

# Alternative method for running ROKU in parallel
ROKU.2 <- function(m, cores) {
#   m <- as.matrix(m)
  ROKU.run <- function(x) {
    #   x <- t(x)
    x.prime <- preprocess_with_tukey(x)
    H <- entropy(abs(x.prime[!is.na(x.prime)]), unit="log2")
    O <- find_tissue_outliers(x)
    values <- cbind(Entropy=H, Entropy.Normalized=H/log2(sum(!is.na(x))), 
                    Outlier.Detection.Method=O$Outlier.Detection.Method,
                    t(O$Model*x.prime))
    return(values)
  }  
  
  R <- mclapply(1:nrow(m), function(i) ROKU.run(m[i,]), mc.cores=cores)
  R <- do.call("rbind", R)
#   dimnames(R) <- list(rownames(m), 
#                       c("Entropy", "Entropy.Normalized", colnames(m)))
  
  return(as.data.frame(R))
}



#' Convert ROKU results as data frame
#' 
#' Take ROKU list and convert to data frame
#' 
#' @param r a list returned by \code{\link{ROKU}}
#' @param meta a optional data frame of additional metadata (column-wise) that will
#' be appended in front of the resulting data frame 
#' @return a data frame
#' @export
#' @seealso \code{\link{ROKU}}
#' @examples
#' rk <- ROKU(psidata)
#' df <- ROKU_as_df
#' head(df)
ROKU_as_df <- function(r, meta=NULL) {
  df <- data.frame(Entropy=r$Entropy,
                   Entropy.Normalized=r$Entropy.Normalized,
                   Outlier.Detection.Method=r$Outlier.Detection.Method,
                   r$Outliers)
  if (!is.null(meta)) {
    if (nrow(meta)!=nrow(r$outliers)) {
      stop("Number of rows don't match")
    }
    df <- data.frame(meta, df)  
  }
  return(df)
}
