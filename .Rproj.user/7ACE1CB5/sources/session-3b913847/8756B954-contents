get.des <- function(dd, na.option = TRUE) {
  desc_stats <- data.frame(
    mean = apply(dd, 2, mean, na.rm = TRUE),
    # mean
    sd = apply(dd, 2, sd, na.rm = TRUE),
    # Standard deviation
    min = apply(dd, 2, min, na.rm = TRUE),
    # minimum
    max = apply(dd, 2, max, na.rm = TRUE),
    # Maximum
    skew = apply(dd, 2, e1071::skewness, na.rm = TRUE)
    #Q1=apply(dd, 2, quantile,probs=0.25),
    #Med = apply(dd, 2, median), # median
    #Q3=apply(dd, 2, quantile,probs=0.75)
  )
  
}



dataTransfomation = function(rawDataValues) {
  library(MASS)
  
  #rawDataValues=as.data.frame(rawDataValues)
  rawDataValuesNormlized = rawDataValues
  
  desc_stats <- get.des(rawDataValues)
  desc_stats$feature <- rownames(desc_stats)
  ll <- length(desc_stats$feature)
  rownames(desc_stats) <- NULL
  no.var <- dim(desc_stats)[1]
  #desc_stats<-desc_stats[,c("icol","feature","mean","sd","min","max","skew")]
  desc_stats$skew2 <- rep(NA, no.var)
  desc_stats$lambda1 <- rep(NA, no.var) #shape
  desc_stats$lambda2 <- rep(NA, no.var) #location
  desc_stats$note <- rep(NA, no.var)
  
  
  for (i in 1:nrow(desc_stats)) {
    #print(desc_stats[i,])
    min.value <- desc_stats$min[i]
    sd.value <- desc_stats$sd[i]
    skew.value <- abs(desc_stats$skew[i])
    
    
    if (sd.value == 0 | abs(skew.value) <= 1 | is.na(skew.value)) {
      desc_stats$note[i] <- 'Raw data'
    } else {
      if (min.value > 0) {
        #log(y) transformation
        rawDataValuesNormlized[, i] <- log(rawDataValues[, i])
        desc_stats$note[i] <- "Log trans."
        s2 <- e1071::skewness(rawDataValuesNormlized[, i], na.rm = TRUE)
        desc_stats$skew2[i] <- s2
        
        if (abs(s2) > 1) {
          # Boxcox transformation was used instead of log trans. if log trans. do not work well
          value <- rawDataValues[, i]
          out.boxcox <- boxcox(value ~ 1, plotit = FALSE)
          index.lambda <- which(out.boxcox$y == max(out.boxcox$y))
          lambda <- out.boxcox$x[index.lambda]
          desc_stats$lambda1[i] <- lambda
          rawDataValuesNormlized[, i] <-
            (((rawDataValues[, i]) ^ lambda) - 1) / lambda
          desc_stats$note[i] <- "BoxCox (1)"
        }
      }  else {
        lamb.TF <- TRUE
        tt <-
          try(geoR::boxcoxfit(rawDataValues[, i],
                        lambda = 0,
                        lambda2 = lamb.TF),
              silent = TRUE)
        #where lambda=0  as initial ==>logtransformation
        
        if (class(tt) != "try-error") {
          lambda <- tt$lambda[[1]] #this is lambda
          lambda2 <- tt$lambda[[2]]
          desc_stats$lambda1[i] <- lambda
          desc_stats$lambda2[i] <- lambda2
          
          if (lambda != 0) {
            rawDataValuesNormlized[, i] <-
              (((rawDataValues[, i] + lambda2) ^ lambda) - 1) / lambda
            desc_stats$note[i] <- "Boxcox trans.(2)"
          } else {
            rawDataValuesNormlized[, i] <-
              log(rawDataValues[, i] + lambda2)
            desc_stats$note[i] <- "Boxcox trans.(2) with log"
          }
          
          
          
        }   else
          ##### fail in boxcoxfit transformation
        {
          desc_stats$note[i] <-
            'Raw due to fail in Boxcox trans.(2)' #skew>1 , but fail to transform
        }
        
      }     ###### end of transformation
      
      
    } #end for the first level of if
    
    
  } #end for for
  
  return(list(rawDataValuesNormlized, desc_stats))
  
}
