train_caret_models<-function(train_CV,Y,currTrainIdx,k,classifier_choice){
  
  set.seed(2020)
  
  cvIndex <- createFolds(Y[currTrainIdx], k, returnTrain = T)
  tc <- trainControl(index = cvIndex,
                     method = 'cv', 
                     number = k,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)
 
  tempList <- lapply(classifier_choice,function(x) train(
    make.names(Y) ~., 
    data = cbind(as.data.frame(train_CV),Y[currTrainIdx])%>% rename(Y = `Y[currTrainIdx]`), 
    method = x,
    metric = "ROC",
    trControl = tc,
    preProcess = c("center","scale")  
  ))
  
  names(tempList)<-classifier_choice
  return(tempList)
} 