prepare_training_covariates<-function(clinicalVars,Y_LOL,train_tf_covariates,r){
  
  bp_clin_train_cv<-cbind(train_tf_covariates[clinicalVars],Y_LOL$band_powerLOL$Xr[,1:r])
  fe_clin_train_cv<-cbind(train_tf_covariates[clinicalVars],Y_LOL$freq_entropyLOL$Xr[,1:r])
  fp_clin_train_cv<-cbind(train_tf_covariates[clinicalVars],Y_LOL$freq_pairsLOL$Xr[,1:r])
  mf_clin_train_cv<-cbind(train_tf_covariates[clinicalVars],Y_LOL$med_freqLOL$Xr[,1:r])
  sm_clin_train_cv<-cbind(train_tf_covariates[clinicalVars],Y_LOL$smaLOL$Xr[,1:r])
  wv_clin_train_cv<-cbind(train_tf_covariates[clinicalVars],Y_LOL$waveletsLOL$Xr[,1:r])
  clin_only_train_cv <- train_tf_covariates[clinicalVars]
  bp_only_train_cv<-as.data.frame(Y_LOL$band_powerLOL$Xr[,1:r])
  fe_only_train_cv<-as.data.frame(Y_LOL$freq_entropyLOL$Xr[,1:r])
  fp_only_train_cv<-as.data.frame(Y_LOL$freq_pairsLOL$Xr[,1:r])
  mf_only_train_cv<-as.data.frame(Y_LOL$med_freqLOL$Xr[,1:r])
  sm_only_train_cv<-as.data.frame(Y_LOL$smaLOL$Xr[,1:r])
  wv_only_train_cv<-as.data.frame(Y_LOL$waveletsLOL$Xr[,1:r])

  colnames(bp_clin_train_cv)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
  colnames(fe_clin_train_cv)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
  colnames(fp_clin_train_cv)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
  colnames(mf_clin_train_cv)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
  colnames(sm_clin_train_cv)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
  colnames(wv_clin_train_cv)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
  
  colnames(bp_only_train_cv)<-paste0(rep("MF.",r),1:r)
  colnames(fe_only_train_cv)<-paste0(rep("MF.",r),1:r)
  colnames(fp_only_train_cv)<-paste0(rep("MF.",r),1:r)
  colnames(mf_only_train_cv)<-paste0(rep("MF.",r),1:r)
  colnames(sm_only_train_cv)<-paste0(rep("MF.",r),1:r)
  colnames(wv_only_train_cv)<-paste0(rep("MF.",r),1:r)
  
  tempList <-
    list(
      bp_clin_train_cv,
      fe_clin_train_cv,
      fp_clin_train_cv,
      mf_clin_train_cv,
      sm_clin_train_cv,
      wv_clin_train_cv,
      clin_only_train_cv,
      bp_only_train_cv,
      fe_only_train_cv,
      fp_only_train_cv,
      mf_only_train_cv,
      sm_only_train_cv,
      wv_only_train_cv
    )
  names(tempList)<-c(
    "bp_clin_train_cv",
    "fe_clin_train_cv",
    "fp_clin_train_cv",
    "mf_clin_train_cv",
    "sm_clin_train_cv",
    "wv_clin_train_cv",
    "clin_only_train_cv",
    "bp_only_train_cv",
    "fe_only_train_cv",
    "fp_only_train_cv",
    "mf_only_train_cv",
    "sm_only_train_cv",
    "wv_only_train_cv"
  )
  return(tempList)
} 