train_SVM<-function(train_CV,Y,currTrainIdx,k){
  
  tempList <- lapply(train_CV,function(x) train(
    make.names(Y) ~., 
    data = cbind(as.data.frame(x),Y[currTrainIdx])%>% rename(Y = `Y[currTrainIdx]`), 
    method = "svmRadial",
    metric = "ROC",
    trControl = trainControl("cv", number = k,classProbs=TRUE,summaryFunction = twoClassSummary),
    preProcess = c("center","scale")  
  ))
  names(tempList)<-c(
    "bp_clin_svm",
    "fe_clin_svm",
    "fp_clin_svm",
    "mf_clin_svm",
    "sm_clin_svm",
    "wv_clin_svm",
    "clin_only_svm",
    "bp_only_svm",
    "fe_only_svm",
    "fp_only_svm",
    "mf_only_svm",
    "sm_only_svm",
    "wv_only_svm"
  )
  return(tempList)
} 