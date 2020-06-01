predict_caret_models<-function(X_models,CV_set,Y,IdxSet){
  
  tempList <- lapply(X_models, function(x) x%>%predict(CV_set,type = "prob")%>%dplyr::select(2))
  tempMat<-list.cbind(tempList)
  tempDF<-cbind(Y[IdxSet],as.data.frame(tempMat))
  names(tempDF)<-c('ground_truth',names(X_models))
  
  return(tempDF)
} 