prepare_testing_covariates<-function(motion_features,Y_LOL,test_tf_covariates,currTestIdx,r){
  
  bp_reduced<-motion_features$band_powerFeats[currTestIdx,]%*%Y_LOL$band_powerLOL$A[,1:r]
  fe_reduced<-motion_features$freq_entropyFeats[currTestIdx,]%*%Y_LOL$freq_entropyLOL$A[,1:r]
  fp_reduced<-motion_features$freq_pairsFeats[currTestIdx,]%*%Y_LOL$freq_pairsLOL$A[,1:r]
  mf_reduced<-motion_features$med_freqFeats[currTestIdx,]%*%Y_LOL$med_freqLOL$A[,1:r]
  sm_reduced<-motion_features$smaFeats[currTestIdx,]%*%Y_LOL$smaLOL$A[,1:r]
  wv_reduced<-motion_features$waveletsFeats[currTestIdx,]%*%Y_LOL$waveletsLOL$A[,1:r]
  
  bp_clin_test_cv<-cbind(test_tf_covariates[clinicalVars],bp_reduced)
  fe_clin_test_cv<-cbind(test_tf_covariates[clinicalVars],fe_reduced)
  fp_clin_test_cv<-cbind(test_tf_covariates[clinicalVars],fp_reduced)
  mf_clin_test_cv<-cbind(test_tf_covariates[clinicalVars],mf_reduced)
  sm_clin_test_cv<-cbind(test_tf_covariates[clinicalVars],sm_reduced)
  wv_clin_test_cv<-cbind(test_tf_covariates[clinicalVars],wv_reduced)
  clin_only_test_cv <- test_tf_covariates[clinicalVars]
  bp_only_test_cv<-as.data.frame(bp_reduced)
  fe_only_test_cv<-as.data.frame(fe_reduced)
  fp_only_test_cv<-as.data.frame(fp_reduced)
  mf_only_test_cv<-as.data.frame(mf_reduced)
  sm_only_test_cv<-as.data.frame(sm_reduced)
  wv_only_test_cv<-as.data.frame(wv_reduced)
  
  colnames(bp_clin_test_cv)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
  colnames(fe_clin_test_cv)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
  colnames(fp_clin_test_cv)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
  colnames(mf_clin_test_cv)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
  colnames(sm_clin_test_cv)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
  colnames(wv_clin_test_cv)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
  
  colnames(bp_only_test_cv)<-paste0(rep("MF.",r),1:r)
  colnames(fe_only_test_cv)<-paste0(rep("MF.",r),1:r)
  colnames(fp_only_test_cv)<-paste0(rep("MF.",r),1:r)
  colnames(mf_only_test_cv)<-paste0(rep("MF.",r),1:r)
  colnames(sm_only_test_cv)<-paste0(rep("MF.",r),1:r)
  colnames(wv_only_test_cv)<-paste0(rep("MF.",r),1:r)
  
  tempList <-
    list(
      bp_clin_test_cv,
      fe_clin_test_cv,
      fp_clin_test_cv,
      mf_clin_test_cv,
      sm_clin_test_cv,
      wv_clin_test_cv,
      clin_only_test_cv,
      bp_only_test_cv,
      fe_only_test_cv,
      fp_only_test_cv,
      mf_only_test_cv,
      sm_only_test_cv,
      wv_only_test_cv
    )
  names(tempList)<-c(
    "bp_clin_test_cv",
    "fe_clin_test_cv",
    "fp_clin_test_cv",
    "mf_clin_test_cv",
    "sm_clin_test_cv",
    "wv_clin_test_cv",
    "clin_only_test_cv",
    "bp_only_test_cv",
    "fe_only_test_cv",
    "fp_only_test_cv",
    "mf_only_test_cv",
    "sm_only_test_cv",
    "wv_only_test_cv"
  )
  return(tempList)
} 