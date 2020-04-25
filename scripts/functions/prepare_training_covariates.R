prepare_training_covariates<-function(Y_LOL,train_tf_covariates,gose_logical){

  if (gose_logical == TRUE){
    band_power_train_cv<-cbind(train_tf_covariates$Age,train_tf_covariates$Sex,train_tf_covariates$GCS_en,Y_LOL$band_powerLOL$Xr[,1:3])
    freq_entropy_train_cv<-cbind(train_tf_covariates$Age,train_tf_covariates$Sex,train_tf_covariates$GCS_en,Y_LOL$freq_entropyLOL$Xr[,1:3])
    freq_pairs_train_cv<-cbind(train_tf_covariates$Age,train_tf_covariates$Sex,train_tf_covariates$GCS_en,Y_LOL$freq_pairsLOL$Xr[,1:3])
    med_freq_train_cv<-cbind(train_tf_covariates$Age,train_tf_covariates$Sex,train_tf_covariates$GCS_en,Y_LOL$med_freqLOL$Xr[,1:3])
    sma_train_cv<-cbind(train_tf_covariates$Age,train_tf_covariates$Sex,train_tf_covariates$GCS_en,Y_LOL$smaLOL$Xr[,1:3])
    wavelets_train_cv<-cbind(train_tf_covariates$Age,train_tf_covariates$Sex,train_tf_covariates$GCS_en,Y_LOL$waveletsLOL$Xr[,1:3])
    
    colnames(band_power_train_cv)<-c("Age","Sex","GCS_en","MF.1","MF.2","MF.3")
    colnames(freq_entropy_train_cv)<-c("Age","Sex","GCS_en","MF.1","MF.2","MF.3")
    colnames(freq_pairs_train_cv)<-c("Age","Sex","GCS_en","MF.1","MF.2","MF.3")
    colnames(med_freq_train_cv)<-c("Age","Sex","GCS_en","MF.1","MF.2","MF.3")
    colnames(sma_train_cv)<-c("Age","Sex","GCS_en","MF.1","MF.2","MF.3")
    colnames(wavelets_train_cv)<-c("Age","Sex","GCS_en","MF.1","MF.2","MF.3")
  } else {
    band_power_train_cv<-cbind(train_tf_covariates$Age,train_tf_covariates$Sex,train_tf_covariates$APACHE,Y_LOL$band_powerLOL$Xr[,1:3])
    freq_entropy_train_cv<-cbind(train_tf_covariates$Age,train_tf_covariates$Sex,train_tf_covariates$APACHE,Y_LOL$freq_entropyLOL$Xr[,1:3])
    freq_pairs_train_cv<-cbind(train_tf_covariates$Age,train_tf_covariates$Sex,train_tf_covariates$APACHE,Y_LOL$freq_pairsLOL$Xr[,1:3])
    med_freq_train_cv<-cbind(train_tf_covariates$Age,train_tf_covariates$Sex,train_tf_covariates$APACHE,Y_LOL$med_freqLOL$Xr[,1:3])
    sma_train_cv<-cbind(train_tf_covariates$Age,train_tf_covariates$Sex,train_tf_covariates$APACHE,Y_LOL$smaLOL$Xr[,1:3])
    wavelets_train_cv<-cbind(train_tf_covariates$Age,train_tf_covariates$Sex,train_tf_covariates$APACHE,Y_LOL$waveletsLOL$Xr[,1:3])
    
    colnames(band_power_train_cv)<-c("Age","Sex","APACHE","MF.1","MF.2","MF.3")
    colnames(freq_entropy_train_cv)<-c("Age","Sex","APACHE","MF.1","MF.2","MF.3")
    colnames(freq_pairs_train_cv)<-c("Age","Sex","APACHE","MF.1","MF.2","MF.3")
    colnames(med_freq_train_cv)<-c("Age","Sex","APACHE","MF.1","MF.2","MF.3")
    colnames(sma_train_cv)<-c("Age","Sex","APACHE","MF.1","MF.2","MF.3")
    colnames(wavelets_train_cv)<-c("Age","Sex","APACHE","MF.1","MF.2","MF.3")
  }
  
  tempList<-list(band_power_train_cv,freq_entropy_train_cv,freq_pairs_train_cv,med_freq_train_cv,sma_train_cv,wavelets_train_cv)
  names(tempList)<-c("band_power_train_cv","freq_entropy_train_cv","freq_pairs_train_cv","med_freq_train_cv","sma_train_cv","wavelets_train_cv")
  return(tempList)
} 