train_kNN<-function(train_CV,Y,currTrainIdx,k){
  tempList <- lapply(train_CV,function(x) train(
    make.names(Y) ~., 
    data = cbind(as.data.frame(x),Y[currTrainIdx])%>% rename(Y = `Y[currTrainIdx]`), 
    method = "knn",
    metric = "ROC",
    trControl = trainControl("cv", number = k,classProbs=TRUE,summaryFunction = twoClassSummary),
    preProcess = c("center","scale")  
  ))
  names(tempList)<-c(
    "bp_clin_knn",
    "fe_clin_knn",
    "fp_clin_knn",
    "mf_clin_knn",
    "sm_clin_knn",
    "wv_clin_knn",
    "clin_only_knn",
    "bp_only_knn",
    "fe_only_knn",
    "fp_only_knn",
    "mf_only_knn",
    "sm_only_knn",
    "wv_only_knn"
  )
  return(tempList)
} 