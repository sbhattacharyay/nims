predict_SVM<-function(X_SVM,CV_set,Y,IdxSet,train_logical){
  
  if (train_logical == TRUE){
    
    bp_pred<-X_SVM$bp_svm%>%predict(CV_set$band_power_train_cv,type = "prob")%>%dplyr::select(2)
    fe_pred<-X_SVM$fe_svm%>%predict(CV_set$freq_entropy_train_cv,type = "prob")%>%dplyr::select(2)
    fp_pred<-X_SVM$fp_svm%>%predict(CV_set$freq_pairs_train_cv,type = "prob")%>%dplyr::select(2)
    mf_pred<-X_SVM$mf_svm%>%predict(CV_set$med_freq_train_cv,type = "prob")%>%dplyr::select(2)
    sm_pred<-X_SVM$sm_svm%>%predict(CV_set$sma_train_cv,type = "prob")%>%dplyr::select(2)
    wv_pred<-X_SVM$wv_svm%>%predict(CV_set$wavelets_train_cv,type = "prob")%>%dplyr::select(2)
    
  } else {
    
    bp_pred<-X_SVM$bp_svm%>%predict(CV_set$band_power_test_cv,type = "prob")%>%dplyr::select(2)
    fe_pred<-X_SVM$fe_svm%>%predict(CV_set$freq_entropy_test_cv,type = "prob")%>%dplyr::select(2)
    fp_pred<-X_SVM$fp_svm%>%predict(CV_set$freq_pairs_test_cv,type = "prob")%>%dplyr::select(2)
    mf_pred<-X_SVM$mf_svm%>%predict(CV_set$med_freq_test_cv,type = "prob")%>%dplyr::select(2)
    sm_pred<-X_SVM$sm_svm%>%predict(CV_set$sma_test_cv,type = "prob")%>%dplyr::select(2)
    wv_pred<-X_SVM$wv_svm%>%predict(CV_set$wavelets_test_cv,type = "prob")%>%dplyr::select(2)
    
  }
  
  tempDF<-cbind(Y[IdxSet],as.data.frame(cbind(bp_pred,fe_pred,fp_pred,mf_pred,sm_pred,wv_pred)))
  colnames(tempDF)<-c("ground_truth","bp_pred","fe_pred","fp_pred","mf_pred","sm_pred","wv_pred")
  
  return(tempDF)
} 