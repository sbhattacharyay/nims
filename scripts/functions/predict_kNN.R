predict_kNN<-function(X_kNN,CV_set,Y,IdxSet){
  
  tempList <- mapply(function(x,y) x%>%predict(y,type = "prob")%>%dplyr::select(2),X_kNN,CV_set)
  tempMat<-list.cbind(tempList)
  
  tempDF<-cbind(Y[IdxSet],as.data.frame(tempMat))%>%rename(ground_truth=1)
  return(tempDF)
} 