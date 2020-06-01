predict_GLM<-function(X_GLM,CV_set,Y,IdxSet){
  
  tempMat <- mapply(function(x,y) predict(x,newx = data.matrix(y),type = "response",s="lambda.min"),X_GLM,CV_set)
  tempDF<-cbind(Y[IdxSet],as.data.frame(tempMat))%>%rename(ground_truth=1)
  return(tempDF)
} 