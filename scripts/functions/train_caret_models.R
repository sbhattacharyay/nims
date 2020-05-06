train_caret_models<-function(train_CV,Y,currTrainIdx,k,classifier_choice){
  
  tempList <- lapply(classifier_choice,function(x) train(
    make.names(Y) ~., 
    data = cbind(as.data.frame(train_CV),Y[currTrainIdx])%>% rename(Y = `Y[currTrainIdx]`), 
    method = x,
    metric = "ROC",
    trControl = trainControl("cv", number = k,classProbs=TRUE,summaryFunction = twoClassSummary),
    preProcess = c("center","scale")  
  ))
  
  names(tempList)<-classifier_choice
  return(tempList)
} 