prepare_testing_covariates<-function(motion_features,Y_LOL,test_tf_covariates,currTestIdx,gose_logical){
  
  bp_reduced<-motion_features$band_powerFeats[currTestIdx,]%*%Y_LOL$band_powerLOL$A[,1:3]
  fe_reduced<-motion_features$freq_entropyFeats[currTestIdx,]%*%Y_LOL$freq_entropyLOL$A[,1:3]
  fp_reduced<-motion_features$freq_pairsFeats[currTestIdx,]%*%Y_LOL$freq_pairsLOL$A[,1:3]
  mf_reduced<-motion_features$med_freqFeats[currTestIdx,]%*%Y_LOL$med_freqLOL$A[,1:3]
  sm_reduced<-motion_features$smaFeats[currTestIdx,]%*%Y_LOL$smaLOL$A[,1:3]
  wv_reduced<-motion_features$waveletsFeats[currTestIdx,]%*%Y_LOL$waveletsLOL$A[,1:3]

  if (gose_logical == TRUE){
    band_power_test_cv<-cbind(test_tf_covariates$Age,test_tf_covariates$Sex,test_tf_covariates$GCS_en,bp_reduced)
    freq_entropy_test_cv<-cbind(test_tf_covariates$Age,test_tf_covariates$Sex,test_tf_covariates$GCS_en,fe_reduced)
    freq_pairs_test_cv<-cbind(test_tf_covariates$Age,test_tf_covariates$Sex,test_tf_covariates$GCS_en,fp_reduced)
    med_freq_test_cv<-cbind(test_tf_covariates$Age,test_tf_covariates$Sex,test_tf_covariates$GCS_en,mf_reduced)
    sma_test_cv<-cbind(test_tf_covariates$Age,test_tf_covariates$Sex,test_tf_covariates$GCS_en,sm_reduced)
    wavelets_test_cv<-cbind(test_tf_covariates$Age,test_tf_covariates$Sex,test_tf_covariates$GCS_en,wv_reduced)
    
    colnames(band_power_test_cv)<-c("Age","Sex","GCS_en","MF.1","MF.2","MF.3")
    colnames(freq_entropy_test_cv)<-c("Age","Sex","GCS_en","MF.1","MF.2","MF.3")
    colnames(freq_pairs_test_cv)<-c("Age","Sex","GCS_en","MF.1","MF.2","MF.3")
    colnames(med_freq_test_cv)<-c("Age","Sex","GCS_en","MF.1","MF.2","MF.3")
    colnames(sma_test_cv)<-c("Age","Sex","GCS_en","MF.1","MF.2","MF.3")
    colnames(wavelets_test_cv)<-c("Age","Sex","GCS_en","MF.1","MF.2","MF.3")
  } else {
    band_power_test_cv<-cbind(test_tf_covariates$Age,test_tf_covariates$Sex,test_tf_covariates$APACHE,bp_reduced)
    freq_entropy_test_cv<-cbind(test_tf_covariates$Age,test_tf_covariates$Sex,test_tf_covariates$APACHE,fe_reduced)
    freq_pairs_test_cv<-cbind(test_tf_covariates$Age,test_tf_covariates$Sex,test_tf_covariates$APACHE,fp_reduced)
    med_freq_test_cv<-cbind(test_tf_covariates$Age,test_tf_covariates$Sex,test_tf_covariates$APACHE,mf_reduced)
    sma_test_cv<-cbind(test_tf_covariates$Age,test_tf_covariates$Sex,test_tf_covariates$APACHE,sm_reduced)
    wavelets_test_cv<-cbind(test_tf_covariates$Age,test_tf_covariates$Sex,test_tf_covariates$APACHE,wv_reduced)
    
    colnames(band_power_test_cv)<-c("Age","Sex","APACHE","MF.1","MF.2","MF.3")
    colnames(freq_entropy_test_cv)<-c("Age","Sex","APACHE","MF.1","MF.2","MF.3")
    colnames(freq_pairs_test_cv)<-c("Age","Sex","APACHE","MF.1","MF.2","MF.3")
    colnames(med_freq_test_cv)<-c("Age","Sex","APACHE","MF.1","MF.2","MF.3")
    colnames(sma_test_cv)<-c("Age","Sex","APACHE","MF.1","MF.2","MF.3")
    colnames(wavelets_test_cv)<-c("Age","Sex","APACHE","MF.1","MF.2","MF.3")
  }
  
  tempList<-list(band_power_test_cv,freq_entropy_test_cv,freq_pairs_test_cv,med_freq_test_cv,sma_test_cv,wavelets_test_cv)
  names(tempList)<-c("band_power_test_cv","freq_entropy_test_cv","freq_pairs_test_cv","med_freq_test_cv","sma_test_cv","wavelets_test_cv")
  return(tempList)
} 