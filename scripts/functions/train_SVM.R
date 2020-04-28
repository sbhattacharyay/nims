train_SVM<-function(train_CV,Y,currTrainIdx,k){
  
  bp_svm <- train(
    make.names(Y) ~., 
    data = cbind(as.data.frame(train_CV$band_power_train_cv),Y[currTrainIdx])%>% rename(Y = `Y[currTrainIdx]`), 
    method = "svmRadial",
    metric = "ROC",
    trControl = trainControl("cv", number = k,classProbs=TRUE,summaryFunction = twoClassSummary),
    preProcess = c("center","scale")  
    )
  
  fe_svm <- train(
    make.names(Y) ~., 
    data = cbind(as.data.frame(train_CV$freq_entropy_train_cv),Y[currTrainIdx])%>% rename(Y = `Y[currTrainIdx]`), 
    method = "svmRadial",
    metric = "ROC",
    trControl = trainControl("cv", number = k,classProbs=TRUE,summaryFunction = twoClassSummary),
    preProcess = c("center","scale")  
  )
  
  fp_svm <- train(
    make.names(Y) ~., 
    data = cbind(as.data.frame(train_CV$freq_pairs_train_cv),Y[currTrainIdx])%>% rename(Y = `Y[currTrainIdx]`), 
    method = "svmRadial",
    metric = "ROC",
    trControl = trainControl("cv", number = k,classProbs=TRUE,summaryFunction = twoClassSummary),
    preProcess = c("center","scale")  
  )
  
  mf_svm <- train(
    make.names(Y) ~., 
    data = cbind(as.data.frame(train_CV$med_freq_train_cv),Y[currTrainIdx])%>% rename(Y = `Y[currTrainIdx]`), 
    method = "svmRadial",
    metric = "ROC",
    trControl = trainControl("cv", number = k,classProbs=TRUE,summaryFunction = twoClassSummary),
    preProcess = c("center","scale")  
  )
  
  sm_svm <- train(
    make.names(Y) ~., 
    data = cbind(as.data.frame(train_CV$sma_train_cv),Y[currTrainIdx])%>% rename(Y = `Y[currTrainIdx]`), 
    method = "svmRadial",
    metric = "ROC",
    trControl = trainControl("cv", number = k,classProbs=TRUE,summaryFunction = twoClassSummary),
    preProcess = c("center","scale")  
  )
  
  wv_svm <- train(
    make.names(Y) ~., 
    data = cbind(as.data.frame(train_CV$wavelets_train_cv),Y[currTrainIdx])%>% rename(Y = `Y[currTrainIdx]`), 
    method = "svmRadial",
    metric = "ROC",
    trControl = trainControl("cv", number = k,classProbs=TRUE,summaryFunction = twoClassSummary),
    preProcess = c("center","scale")  
  )

  tempList<-list(bp_svm,fe_svm,fp_svm,mf_svm,sm_svm,wv_svm)
  names(tempList)<-c("bp_svm","fe_svm","fp_svm","mf_svm","sm_svm","wv_svm")
  return(tempList)
} 