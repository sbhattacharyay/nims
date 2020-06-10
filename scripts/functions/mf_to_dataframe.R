mf_to_dataframe <- function(motion_features,n,verbose = TRUE) {
  totallyMissingSet <- data.frame(matrix(ncol = 4, nrow = 0))
  names(totallyMissingSet) <- c("ptIdx","ftIdx","srIdx","mfIdx")
  
  # Here we also compile a complete set of all the recorded motion features across the set
  completeFeatureSet <- data.frame(matrix(ncol = 10, nrow = 0))
  names(completeFeatureSet) <- c("Bed","LA","LE","LW","RA","RE","RW","featureType","ptIdx","timeCount")
  
  for (mfIdx in 1:length(motion_features)) {
    curr_matrix <- ts(t(motion_features[[mfIdx]]),frequency = 720)
    colnames(curr_matrix) <- c("Bed","LA","LE","LW","RA","RE","RW")
    curr_matrix[is.nan(curr_matrix)] <- NA
    curr_matrix[is.infinite(curr_matrix)] <- NA
    if(mfIdx %% n == 0){
      ptIdx <- n
    } else {
      ptIdx <- mfIdx %% n
    }
    ftIdx <- ceiling(mfIdx/n)
    for (srIdx in 1:ncol(curr_matrix)) {
      if (all(curr_matrix[,srIdx] == 0,na.rm = TRUE)){
        curr_matrix[,srIdx] <- rep(NA, length(curr_matrix[,srIdx]))
      }
      if (all(is.na(curr_matrix[,srIdx]))) {
        totallyMissingSet<-rbind(totallyMissingSet,data.frame(ptIdx,ftIdx,srIdx,mfIdx))
      }
    }
    tempDF <- as.data.frame(curr_matrix)
    tempDF$featureType <- featureLabels[,ftIdx]
    tempDF$ptIdx <- ptIdx
    tempDF$timeCount <- 1:nrow(curr_matrix)
    completeFeatureSet <- rbind(completeFeatureSet,tempDF)
    if (verbose == TRUE){
      print(paste('set no',mfIdx,'complete'))
    }
  }
  return(list(completeFeatureSet,totallyMissingSet))
}