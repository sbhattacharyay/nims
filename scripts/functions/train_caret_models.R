train_caret_models<-function(train_CV,Y,currTrainIdx,k,classifier_choice){
  
  # This function tunes the hyperparameters for the ML models specified in `classifier_choice` and subsequently trains the tuned models on the training data
  
  set.seed(2020)
  
  cvIndex <- createFolds(Y[currTrainIdx], k, returnTrain = T)
  tc <- trainControl(index = cvIndex,
                     method = 'cv', 
                     number = k,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary,
                     verboseIter = FALSE,
                     savePredictions = "all",
                     returnResamp = "all")
 
  tempList <- lapply(classifier_choice,function(x) train(
    make.names(Y) ~., 
    data = cbind(as.data.frame(train_CV),Y[currTrainIdx])%>% rename(Y = `Y[currTrainIdx]`), 
    method = x,
    metric = "ROC",
    trControl = tc,
    preProcess = c("YeoJohnson","center","scale")  
  ))
  
  names(tempList)<-classifier_choice
  return(tempList)
} 