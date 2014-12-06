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