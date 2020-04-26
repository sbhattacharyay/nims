predict_GLM<-function(X_GLM,CV_set,Y,IdxSet,train_logical){
  
  if (train_logical == TRUE){
    bp_pred<-predict(X_GLM$bp_glm,newx = CV_set$band_power_train_cv,type = "response",s="lambda.min")
    fe_pred<-predict(X_GLM$fe_glm,newx = CV_set$freq_entropy_train_cv,type = "response",s="lambda.min")
    fp_pred<-predict(X_GLM$fp_glm,newx = CV_set$freq_pairs_train_cv,type = "response",s="lambda.min")
    mf_pred<-predict(X_GLM$mf_glm,newx = CV_set$med_freq_train_cv,type = "response",s="lambda.min")
    sm_pred<-predict(X_GLM$sm_glm,newx = CV_set$sma_train_cv,type = "response",s="lambda.min")
    wv_pred<-predict(X_GLM$wv_glm,newx = CV_set$wavelets_train_cv,type = "response",s="lambda.min")
  } else {
    bp_pred<-predict(X_GLM$bp_glm,newx = CV_set$band_power_test_cv,type = "response",s="lambda.min")
    fe_pred<-predict(X_GLM$fe_glm,newx = CV_set$freq_entropy_test_cv,type = "response",s="lambda.min")
    fp_pred<-predict(X_GLM$fp_glm,newx = CV_set$freq_pairs_test_cv,type = "response",s="lambda.min")
    mf_pred<-predict(X_GLM$mf_glm,newx = CV_set$med_freq_test_cv,type = "response",s="lambda.min")
    sm_pred<-predict(X_GLM$sm_glm,newx = CV_set$sma_test_cv,type = "response",s="lambda.min")
    wv_pred<-predict(X_GLM$wv_glm,newx = CV_set$wavelets_test_cv,type = "response",s="lambda.min")
  }

  tempDF<-cbind(Y[IdxSet],as.data.frame(cbind(bp_pred,fe_pred,fp_pred,mf_pred,sm_pred,wv_pred)))
  colnames(tempDF)<-c("ground_truth","bp_pred","fe_pred","fp_pred","mf_pred","sm_pred","wv_pred")
  
  return(tempDF)
} 