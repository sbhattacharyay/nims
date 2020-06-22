# Requires `caret` package

createOuterCVFolds <-function(patient_clinical_data,outer_fold_count,seed=2020){
  
  avail_at_dis <- which(!is.na(patient_clinical_data$fav_mort_dis))
  avail_at_12m <- which(!is.na(patient_clinical_data$fav_mort_12m))
  
  set.seed(seed = seed)
  outer_cvIdx_mort_dis <- createFolds(patient_clinical_data$fav_mort_dis,outer_fold_count,list = TRUE)
  outer_cvIdx_mort_12m <- createFolds(patient_clinical_data$fav_mort_12m[avail_at_12m],outer_fold_count,list = TRUE)
  outer_cvIdx_GOSE_dis <- createFolds(patient_clinical_data$fav_GOSE_dis,outer_fold_count,list = TRUE)
  outer_cvIdx_GOSE_12m <- createFolds(patient_clinical_data$fav_GOSE_12m[avail_at_12m],outer_fold_count,list = TRUE)
  outer_cvIdx_mRS_dis  <- createFolds(patient_clinical_data$fav_mRS_dis,outer_fold_count,list = TRUE)
  outer_cvIdx_mRS_12m  <- createFolds(patient_clinical_data$fav_mRS_12m[avail_at_12m],outer_fold_count,list = TRUE)
  
  outerFolds <- list(outer_cvIdx_mort_dis,outer_cvIdx_mort_12m,outer_cvIdx_GOSE_dis,outer_cvIdx_GOSE_12m,outer_cvIdx_mRS_dis,outer_cvIdx_mRS_12m)
  names(outerFolds) <- c("fav_mort_dis","fav_mort_12m","fav_GOSE_dis","fav_GOSE_12m","fav_mRS_dis","fav_mRS_12m")
  
  return(outerFolds)
}