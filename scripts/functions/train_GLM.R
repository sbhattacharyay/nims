train_GLM<-function(train_CV,Y,currTrainIdx,a,k){

    bp_glm<-cv.glmnet(train_CV$band_power_train_cv,Y[currTrainIdx],alpha = a,family = "binomial",type.measure = "auc",nfolds=k)
    fe_glm<-cv.glmnet(train_CV$band_power_train_cv,Y[currTrainIdx],alpha = a,family = "binomial",type.measure = "auc",nfolds=k)
    fp_glm<-cv.glmnet(train_CV$band_power_train_cv,Y[currTrainIdx],alpha = a,family = "binomial",type.measure = "auc",nfolds=k)
    mf_glm<-cv.glmnet(train_CV$band_power_train_cv,Y[currTrainIdx],alpha = a,family = "binomial",type.measure = "auc",nfolds=k)
    sm_glm<-cv.glmnet(train_CV$band_power_train_cv,Y[currTrainIdx],alpha = a,family = "binomial",type.measure = "auc",nfolds=k)
    wv_glm<-cv.glmnet(train_CV$band_power_train_cv,Y[currTrainIdx],alpha = a,family = "binomial",type.measure = "auc",nfolds=k)
    
    tempList<-list(bp_glm,fe_glm,fp_glm,mf_glm,sm_glm,wv_glm)
    names(tempList)<-c("bp_glm","fe_glm","fp_glm","mf_glm","sm_glm","wv_glm")
    return(tempList)
} 