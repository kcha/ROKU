preprocess_with_tukey <- function(x) {
  # Step 1 of ROKU
  tu <- tukey.biweight(x[!is.na(x)])
  new.x <- x - tu
  return(new.x)
}