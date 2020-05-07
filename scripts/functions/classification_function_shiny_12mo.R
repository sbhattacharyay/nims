# Limit dataset to patients with 12 month information available
idx_for_12mo <- !is.na(patient_clinical_data$favorable_12mo)
seq_for_12mo <- seq(62)[idx_for_12mo]

patient_clinical_data_yr<-patient_clinical_data[seq_for_12mo,]
tod_tf_covariates_yr<-tod_tf_covariates[seq_for_12mo,]
tfr_tf_covariates_yr<-tfr_tf_covariates[seq_for_12mo,]
tod_motion_features_yr<-lapply(tod_motion_features, function(x) x[seq_for_12mo,])
tfr_motion_features_yr<-lapply(tfr_motion_features, function(x) x[seq_for_12mo,])

# Split for k-fold cross validation for discharge predictions:
k <- 5
cvIdx <- cross_val_splits(patient_clinical_data_yr, k)

classification_function_shiny_12mo <-function(time_choice,time_slide,classifier_choice,r,mf_choice,clinicalVars,sensor_loc){

  sensor_choice_Idxs <- match(sensor_loc,sensor_options)
  
  if (time_choice == 'tod'){
    sensor_idx_ranges<-lapply(tod_sensor_idx_ranges,function(x) x[tod_range >= time_slide[1] & tod_range <= time_slide[2]])
    fp2_idx_ranges   <-lapply(tod_fp2_idx_ranges,function(x) x[tod_range >= time_slide[1] & tod_range <= time_slide[2]])
    curr_sensor_Idx<-do.call(c, sensor_idx_ranges[sensor_choice_Idxs])
    curr_fp_Idx<-c(curr_sensor_Idx,do.call(c, fp2_idx_ranges[sensor_choice_Idxs]))
    curr_motion_features <- lapply(tod_motion_features_yr[mf_choice], function(x) x[,curr_sensor_Idx])
    
    if ('freq_pairsFeats' %in% mf_choice) curr_motion_features$freq_pairsFeats <- tod_motion_features_yr$freq_pairsFeats[,curr_fp_Idx]
    
    curr_motion_features <- do.call(cbind,curr_motion_features)
    tf_covariates <- tod_tf_covariates_yr
    
  } else {
    sensor_idx_ranges<-lapply(tfr_sensor_idx_ranges,function(x) x[tfr_range >= time_slide[1] & tfr_range <= time_slide[2]])
    fp2_idx_ranges   <-lapply(tfr_fp2_idx_ranges,function(x) x[tfr_range >= time_slide[1] & tfr_range <= time_slide[2]])
    curr_sensor_Idx<-do.call(c, sensor_idx_ranges[sensor_choice_Idxs])
    curr_fp_Idx<-c(curr_sensor_Idx,do.call(c, fp2_idx_ranges[sensor_choice_Idxs]))
    curr_motion_features <- lapply(tfr_motion_features_yr[mf_choice], function(x) x[,curr_sensor_Idx])
    
    if ('freq_pairsFeats' %in% mf_choice) curr_motion_features$freq_pairsFeats <- tfr_motion_features_yr$freq_pairsFeats[,curr_fp_Idx]
    
    curr_motion_features <- do.call(cbind,curr_motion_features)
    tf_covariates <- tfr_tf_covariates_yr
    
  }  

  GOSE_predictions  <- vector(mode = "list", length = k)
  fav_predictions   <- vector(mode = "list", length = k)
  death_predictions <- vector(mode = "list", length = k)
  
  for (i in 1:length(cvIdx)) {
    currTestIdx <- cvIdx[[as.character(i)]]
    currTrainIdx <- seq(nrow(patient_clinical_data_yr))[-currTestIdx]
    
    # Convert numeric to logical indexing
    logicalTest <- rep(FALSE, nrow(patient_clinical_data_yr))
    logicalTest[currTestIdx] <- TRUE
    logicalTrain <- !logicalTest
    
    train_tf_covariates <- tf_covariates[currTrainIdx,]
    test_tf_covariates  <- tf_covariates[currTestIdx,]
    
    # Perform LOL on training data
    
    # - GOSE LOL:
    Y <- as.factor(patient_clinical_data_yr$gose_12mo)
    logicalIdx <- logicalTrain & !is.na(Y)
    GOSE_LOL<-lol.project.lol(curr_motion_features[logicalIdx,], Y[logicalIdx], r)
    
    # - Favorable Outcome LOL:
    Y <- as.factor(patient_clinical_data_yr$favorable_12mo)
    logicalIdx <- logicalTrain & !is.na(Y)
    fav_LOL<-lol.project.lol(curr_motion_features[logicalIdx,], Y[logicalIdx], r)
    
    # - Mortality Outcome LOL:
    Y <- as.factor(patient_clinical_data_yr$death_12mo)
    logicalIdx <- logicalTrain & !is.na(Y)
    death_LOL<-lol.project.lol(curr_motion_features[logicalIdx,], Y[logicalIdx], r)
    
    # Prepare training covariates
    
    if ("diag" %in% clinicalVars) {
      clinicalVars <- c(clinicalVars[!clinicalVars %in% "diag"],"CVA","ICH","SAH","BT","SDH.TBI")
    }
    
    train_GOSE_CV <- cbind(train_tf_covariates[clinicalVars],GOSE_LOL$Xr[,1:r])
    train_fav_CV  <- cbind(train_tf_covariates[clinicalVars],fav_LOL$Xr[,1:r])
    train_death_CV<- cbind(train_tf_covariates[clinicalVars],death_LOL$Xr[,1:r])
    
    colnames(train_GOSE_CV)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
    colnames(train_fav_CV)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
    colnames(train_death_CV)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
    
    # Prepare testing covariates
    test_GOSE_CV <- cbind(test_tf_covariates[clinicalVars],curr_motion_features[currTestIdx,]%*%GOSE_LOL$A[,1:r])
    test_fav_CV  <- cbind(test_tf_covariates[clinicalVars],curr_motion_features[currTestIdx,]%*%fav_LOL$A[,1:r])
    test_death_CV<- cbind(test_tf_covariates[clinicalVars],curr_motion_features[currTestIdx,]%*%death_LOL$A[,1:r])
    
    colnames(test_GOSE_CV)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
    colnames(test_fav_CV)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
    colnames(test_death_CV)<-c(clinicalVars,paste0(rep("MF.",r),1:r))
    
    # Train models on training data
    Y <- as.factor(patient_clinical_data_yr$favorable_12mo)
    
    GOSE_models <- train_caret_models(train_GOSE_CV,Y,currTrainIdx,4,classifier_choice)
    fav_models  <- train_caret_models(train_fav_CV,Y,currTrainIdx,4,classifier_choice)
    
    Y <- as.factor(patient_clinical_data_yr$death_12mo)
    
    death_models  <- train_caret_models(train_death_CV,Y,currTrainIdx,4,classifier_choice)
    
    # Predict outcomes on testing data using trained models
    
    Y <- as.factor(patient_clinical_data_yr$favorable_12mo)
    
    GOSE_test_val <-predict_caret_models(GOSE_models, test_GOSE_CV, Y, currTestIdx)
    fav_test_val  <-predict_caret_models(fav_models, test_fav_CV, Y, currTestIdx)
    
    Y <- as.factor(patient_clinical_data_yr$death_12mo)
    
    death_test_val<-predict_caret_models(death_models, test_death_CV, Y, currTestIdx)
    
    GOSE_predictions[[i]] <- GOSE_test_val
    fav_predictions[[i]] <- fav_test_val
    death_predictions[[i]] <- death_test_val
  }
  
  outList <- list(GOSE_predictions,fav_predictions,death_predictions)
  names(outList) <- c("GOSE_predictions","fav_predictions","death_predictions")
  return(outList)
}