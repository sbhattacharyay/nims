test_prediction_model <- function(modelOutput,predictor_matrix,trueLabels){
  
  ordProbs <- as.data.frame(matrix(nrow = nrow(predictor_matrix), ncol = length(modelOutput)))
  
  for (labelIdx in 1:length(modelOutput)){
    curr_label_pred <- predict(modelOutput[[labelIdx]], predictor_matrix, type = "prob")
    ordProbs[,labelIdx] <- curr_label_pred$greater
  }
  
  classProbs <- as.data.frame(matrix(nrow = nrow(predictor_matrix), ncol = length(modelOutput)+1))
  
  classProbs[,1] <- 1-ordProbs[,1]
  classProbs[,ncol(classProbs)] <- ordProbs[,ncol(ordProbs)]  
  
  for (cpIdx in 2:(ncol(classProbs)-1)){
    classProbs[,cpIdx] <- ordProbs[,(cpIdx-1)] - ordProbs[,cpIdx]
  }
  
  predLabels <- as.factor(max.col(classProbs, ties.method = "first"))
  predLabels <- plyr::revalue(predLabels,c("1"="decrease","2"="no.change","3"="increase"))
  return(data.frame(trueLabels,predLabels))
}