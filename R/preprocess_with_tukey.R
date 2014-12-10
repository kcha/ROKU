#' Normalize by Tukey biweight
#' 
#' The first step in ROKU algorithm. Normalize the expression data using Tukey's
#' biweight. Values with \code{NA} are excluded.
#' 
#' @param x A vector of length \emph{m} where \emph{m} is the number of samples
#' or tissues
#' @return A transformed vector of length \emph{m}
#' @importFrom affy tukey.biweight
#' @references
#' Kadota et al (2006). ROKU: a novel method for identification of tissue-specific genes. BMC Genomics, 7:294. 
#' \url{http://www.biomedcentral.com/1471-2105/7/294}
#' 
#' @export
preprocess_with_tukey <- function(x) {
  # Step 1 of ROKU
  tu <- tukey.biweight(x[!is.na(x)])
  new.x <- x - tu
  return(new.x)
}