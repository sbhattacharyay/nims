source('./functions/load_patient_clinical_data.R')
source('./functions/update_clinicalVariableList.R')
source("./functions/load_tf_patient_covariates.R")
source("./functions/train_caret_models.R")
source("./functions/predict_caret_models.R")

# Load patient clinical data (sorts by PY numbering and corrects variable types)
patient_clinical_data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv')

# Load and update clinical variable list
clinicalVariableList <- update_clinicalVariableList('../clinical_data/clinicalVariableList.csv')

# Load transformed covariates for both TOD and TFR
tod_tf_covariates <-load_tf_patient_covariates('../tod_motion_feature_data/tf_patient_covariates.csv')
tfr_tf_covariates <-load_tf_patient_covariates('../tfr_motion_feature_data/tf_patient_covariates.csv')

tod_sensor_idx_ranges <- split(1:77754, ceiling(seq_along(1:77754)/12959))
tod_fp2_idx_ranges    <- lapply(tod_sensor_idx_ranges, function(x) x+77754)

tfr_sensor_idx_ranges <- split(1:34554, ceiling(seq_along(1:34554)/5759))
tfr_fp2_idx_ranges    <- lapply(tfr_sensor_idx_ranges, function(x) x+34554)

tfr_range<-seq(5/3600,8,length.out=5759)
tod_range<-seq(as.POSIXct("2020-05-05 18:00:05"),as.POSIXct("2020-05-06 12:00:00"),length.out=12959)

sensor_options <- c("left_ank","left_el","left_wr","right_ank","right_el","right_wr")

classification_function <-function(time_choice,time_slide,classifier_choice,outcomes,mf_choice,sensor_loc){
  
  avail_idx <- !is.na(outcomes)
  avail_seq <- seq(length(avail_idx))[avail_idx]
  
  patient_clinical_data_avail<-patient_clinical_data[avail_seq,]
  tod_tf_covariates_avail<-tod_tf_covariates[avail_seq,]
  tfr_tf_covariates_avail<-tfr_tf_covariates[avail_seq,]
  tod_motion_features_avail<-lapply(tod_motion_features, function(x) x[avail_seq,])
  tfr_motion_features_avail<-lapply(tfr_motion_features, function(x) x[avail_seq,])
  outcomes_avail<-outcomes[avail_seq]
  
  set.seed(123)
  
  # Nested stratified cross-validation: Create 5 folds for outer cross-validation (training vs. testing)
  k <- 5
  outer_cvIdx <- createFolds(outcomes_avail,k,list = TRUE);
  
  sensor_choice_Idxs <- match(sensor_loc,sensor_options)
  
  if (time_choice == 'tod'){
    sensor_idx_ranges<-lapply(tod_sensor_idx_ranges,function(x) x[tod_range >= time_slide[1] & tod_range <= time_slide[2]])
    fp2_idx_ranges   <-lapply(tod_fp2_idx_ranges,function(x) x[tod_range >= time_slide[1] & tod_range <= time_slide[2]])
    curr_sensor_Idx<-do.call(c, sensor_idx_ranges[sensor_choice_Idxs])
    curr_fp_Idx<-c(curr_sensor_Idx,do.call(c, fp2_idx_ranges[sensor_choice_Idxs]))
    curr_motion_features <- lapply(tod_motion_features_avail[mf_choice], function(x) x[,curr_sensor_Idx])
    
    if ('freq_pairsFeats' %in% mf_choice) curr_motion_features$freq_pairsFeats <- tod_motion_features_avail$freq_pairsFeats[,curr_fp_Idx]
    
    curr_motion_features <- do.call(cbind,curr_motion_features)
    tf_covariates <- tod_tf_covariates_avail
    
  } else {
    sensor_idx_ranges<-lapply(tfr_sensor_idx_ranges,function(x) x[tfr_range >= time_slide[1] & tfr_range <= time_slide[2]])
    fp2_idx_ranges   <-lapply(tfr_fp2_idx_ranges,function(x) x[tfr_range >= time_slide[1] & tfr_range <= time_slide[2]])
    curr_sensor_Idx<-do.call(c, sensor_idx_ranges[sensor_choice_Idxs])
    curr_fp_Idx<-c(curr_sensor_Idx,do.call(c, fp2_idx_ranges[sensor_choice_Idxs]))
    curr_motion_features <- lapply(tfr_motion_features_avail[mf_choice], function(x) x[,curr_sensor_Idx])
    
    if ('freq_pairsFeats' %in% mf_choice) curr_motion_features$freq_pairsFeats <- tfr_motion_features_avail$freq_pairsFeats[,curr_fp_Idx]
    
    curr_motion_features <- do.call(cbind,curr_motion_features)
    tf_covariates <- tfr_tf_covariates_avail
    
  }  
  
  combined_self_preds <- vector(mode = "list", length = k)
  clin_only_self_preds  <- vector(mode = "list", length = k)
  mf_only_self_preds  <- vector(mode = "list", length = k)

  combined_predictions <- vector(mode = "list", length = k)
  clin_only_predictions  <- vector(mode = "list", length = k)
  mf_only_predictions  <- vector(mode = "list", length = k)
  
  for (i in 1:k) {
    currTestIdx <- outer_cvIdx[[i]]
    currTrainIdx <- seq(nrow(patient_clinical_data_avail))[-currTestIdx]
   
    # Convert numeric to logical indexing
    logicalTest <- rep(FALSE, nrow(patient_clinical_data_avail))
    logicalTest[currTestIdx] <- TRUE
    logicalTrain <- !logicalTest

    train_tf_covariates <- tf_covariates[currTrainIdx,]
    test_tf_covariates  <- tf_covariates[currTestIdx,]
    
    # Nested stratified cross-validation: Create 4 folds for inner cross-validation (hyperparameter tuning)
    inner_cvIdx <- createFolds(outcomes_avail[currTrainIdx], 4, returnTrain = T)
    
    # Convert inner cross-validation indices to LOL-compatible cross-validation indices
    lol_sets <- list()
    for (q in 1:length(inner_cvIdx)){
      tempList <- list(train = inner_cvIdx[[q]], 
                       test = seq(length(outcomes_avail[currTrainIdx]))[!seq(length(outcomes_avail[currTrainIdx])) %in% inner_cvIdx[[q]]])
      lol_sets[[as.character(q)]] <- tempList
    }
    
    # Note: on all trials, we found that r_opt = 1.
    # # Tune LOL on the inner cross-validation set to find optimal dimension for reduction:
    # LOL <- lol.xval.optimal_dimselect(curr_motion_features[currTrainIdx,], outcomes_avail[currTrainIdx], rs=c(1:3,5:6,10), lol.project.lol, sets=lol_sets)
    # r <- LOL$optimal.r
    r <- 1
    LOL <- lol.project.lol(curr_motion_features[currTrainIdx,], outcomes_avail[currTrainIdx], r)
    
    # Prepare training covariates
    combined_train_CV<- cbind(train_tf_covariates$APACHE_risk,LOL$Xr[,1:r])
    colnames(combined_train_CV)<-c('APACHE_risk',paste0(rep("MF.",r),1:r))
    clin_only_train_CV<-as.matrix(combined_train_CV[,'APACHE_risk'])
    mf_only_train_CV<-as.matrix(combined_train_CV[,paste0(rep("MF.",r),1:r)])
    
    # Prepare testing covariates
    combined_test_CV<- cbind(test_tf_covariates$APACHE_risk,curr_motion_features[currTestIdx,]%*%LOL$A[,1:r])
    colnames(combined_test_CV)<-c('APACHE_risk',paste0(rep("MF.",r),1:r))
    clin_only_test_CV<-as.matrix(combined_test_CV[,'APACHE_risk'])
    mf_only_test_CV<-as.matrix(combined_test_CV[,paste0(rep("MF.",r),1:r)])
    
    # Train models on training data
    combined_models  <- train_caret_models(combined_train_CV,outcomes_avail,currTrainIdx,inner_cvIdx,classifier_choice)
    if ('glmnet' %in% classifier_choice){
      temp_classifier_choice<-replace(classifier_choice, classifier_choice == 'glmnet', 'glm')
      clin_only_models <- train_caret_models(clin_only_train_CV,outcomes_avail,currTrainIdx,inner_cvIdx,temp_classifier_choice)
    } else {
      clin_only_models <- train_caret_models(clin_only_train_CV,outcomes_avail,currTrainIdx,inner_cvIdx,classifier_choice)
    }
    if ('glmnet' %in% classifier_choice & r==1){
      temp_classifier_choice<-replace(classifier_choice, classifier_choice == 'glmnet', 'glm')
      mf_only_models <- train_caret_models(mf_only_train_CV,outcomes_avail,currTrainIdx,inner_cvIdx,temp_classifier_choice)
    } else {
      mf_only_models <- train_caret_models(mf_only_train_CV,outcomes_avail,currTrainIdx,inner_cvIdx,classifier_choice)
    }
    
    # Determine optimal cutoff of each probability output by minimizing distance to optimal classifier (TPR = 1, FPR = 0) on training data
    combined_self_val<-predict_caret_models(combined_models, combined_train_CV, outcomes_avail, currTrainIdx)
    clin_only_self_val<-predict_caret_models(clin_only_models, clin_only_train_CV, outcomes_avail, currTrainIdx)
    mf_only_self_val<-predict_caret_models(mf_only_models, mf_only_train_CV, outcomes_avail, currTrainIdx)
    
    combined_self_preds[[i]] <- combined_self_val
    clin_only_self_preds[[i]]  <- clin_only_self_val
    mf_only_self_preds[[i]] <- mf_only_self_val
    
    # Predict outcomes on testing data using trained models

    combined_test_val<-predict_caret_models(combined_models, combined_test_CV, outcomes_avail, currTestIdx)
    clin_only_test_val<-predict_caret_models(clin_only_models, clin_only_test_CV, outcomes_avail, currTestIdx)
    mf_only_test_val<-predict_caret_models(mf_only_models, mf_only_test_CV, outcomes_avail, currTestIdx)

    combined_predictions[[i]] <- combined_test_val
    clin_only_predictions[[i]] <- clin_only_test_val
    mf_only_predictions[[i]] <- mf_only_test_val
  }
  
  outList1 <-
    list(
      combined_predictions,
      clin_only_predictions,
      mf_only_predictions
    )
  names(outList1) <-
    c(
      "combined_predictions",
      "clin_only_predictions",
      "mf_only_predictions"
    )
  outList2 <-
    list(
      combined_self_preds,
      clin_only_self_preds,
      mf_only_self_preds
    )
  names(outList2) <-
    c(
      "combined_self_preds",
      "clin_only_self_preds",
      "mf_only_self_preds"
    )
  
  finalOutList <- list(outList1,outList2)
  return(finalOutList)
}