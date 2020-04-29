train_LDA<-function(train_CV,Y,currTrainIdx){
  
  tempList <- lapply(train_CV,function(x) lda(make.names(Y)~., data=cbind(as.data.frame(x),Y[currTrainIdx])%>% rename(Y = `Y[currTrainIdx]`)))
  names(tempList)<-c(
    "bp_clin_lda",
    "fe_clin_lda",
    "fp_clin_lda",
    "mf_clin_lda",
    "sm_clin_lda",
    "wv_clin_lda",
    "clin_only_lda",
    "bp_only_lda",
    "fe_only_lda",
    "fp_only_lda",
    "mf_only_lda",
    "sm_only_lda",
    "wv_only_lda"
  )
  return(tempList)
} 