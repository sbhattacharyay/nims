source('./functions/load_patient_clinical_data.R')
source('./functions/update_clinicalVariableList.R')
source('./functions/get_motion_features.R')
source("./functions/load_tf_patient_covariates.R")
source("./functions/train_caret_models.R")
source("./functions/predict_caret_models.R")

# Load patient clinical data (sorts by PY numbering and corrects variable types)
patient_clinical_data <- load_patient_clinical_data()

# Load and update clinical variable list
clinicalVariableList <- update_clinicalVariableList()

# Load Motion Features Organized by Time of Day (TOD)
if (!exists("tod_motion_features")) {
  tod_sensors <-
    readMat('../tod_motion_feature_data/bed_corrected_imputed_complete_sensor_data.mat')$bed.corrected.sensors
  tod_motion_features <- get_motion_features(do.call(rbind, tod_sensors))
}

# Load Motion Features Organized by Time from Recording (TFR)
if (!exists("tfr_motion_features")) {
  tfr_sensors <-
    readMat('../tfr_motion_feature_data/bed_corrected_imputed_complete_sensor_data.mat')$bed.corrected.sensors
  tfr_motion_features <- get_motion_features(do.call(rbind, tfr_sensors))
}

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

idx_for_12mo <- !is.na(patient_clinical_data$Death12Months)
seq_for_12mo <- seq(length(idx_for_12mo))[idx_for_12mo]

patient_clinical_data_yr<-patient_clinical_data[seq_for_12mo,]
tod_tf_covariates_yr<-tod_tf_covariates[seq_for_12mo,]
tfr_tf_covariates_yr<-tfr_tf_covariates[seq_for_12mo,]
tod_motion_features_yr<-lapply(tod_motion_features, function(x) x[seq_for_12mo,])
tfr_motion_features_yr<-lapply(tfr_motion_features, function(x) x[seq_for_12mo,])

classification_function_shiny_12mo <-function(time_choice,time_slide,classifier_choice,r,mf_choice,sensor_loc){
  
  set.seed(123)
  
  # Split for k-fold cross validation for discharge predictions:
  k <- 5
  death_cvIdx <- createFolds(patient_clinical_data_yr$Death12Months,k,list = TRUE);
  
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
  
  death_predictions <- vector(mode = "list", length = k)
  clin_only_death_predictions  <- vector(mode = "list", length = k)
  mf_only_death_predictions  <- vector(mode = "list", length = k)
  
  for (i in 1:k) {
    death_currTestIdx <- death_cvIdx[[i]]
    death_currTrainIdx <- seq(nrow(patient_clinical_data_yr))[-death_currTestIdx]
    
    # Convert numeric to logical indexing
    death_logicalTest <- rep(FALSE, nrow(patient_clinical_data_yr))
    death_logicalTest[death_currTestIdx] <- TRUE
    death_logicalTrain <- !death_logicalTest
    
    death_train_tf_covariates <- tf_covariates[death_currTrainIdx,]
    death_test_tf_covariates  <- tf_covariates[death_currTestIdx,]
    
    # Perform LOL on training data - Mortality Outcome LOL:
    Y <- as.factor(patient_clinical_data_yr$Death12Months)
    logicalIdx <- death_logicalTrain & !is.na(Y)
    
    death_LOL<-lol.project.lol(curr_motion_features[logicalIdx,], Y[logicalIdx], r)
    
    # Prepare training covariates
    train_death_CV<- cbind(death_train_tf_covariates$APACHE_risk,death_LOL$Xr[,1:r])
    colnames(train_death_CV)<-c('APACHE_risk',paste0(rep("MF.",r),1:r))
    clin_only_train_CV_death<-as.matrix(train_death_CV[,'APACHE_risk'])
    mf_only_train_death_CV<-train_death_CV[,paste0(rep("MF.",r),1:r)]
    
    # Prepare testing covariates
    test_death_CV<- cbind(death_test_tf_covariates$APACHE_risk,curr_motion_features[death_currTestIdx,]%*%death_LOL$A[,1:r])
    colnames(test_death_CV)<-c('APACHE_risk',paste0(rep("MF.",r),1:r))
    clin_only_test_CV_death<-as.matrix(test_death_CV[,'APACHE_risk'])
    mf_only_test_death_CV<-test_death_CV[,paste0(rep("MF.",r),1:r)]
    
    # Train models on training data
    Y <- as.factor(patient_clinical_data_yr$Death12Months)
    death_models  <- train_caret_models(train_death_CV,Y,death_currTrainIdx,4,classifier_choice)
    
    if ('glmnet' %in% classifier_choice){
      temp_classifier_choice<-replace(classifier_choice, classifier_choice == 'glmnet', 'glm')
      clin_only_death_models <- train_caret_models(clin_only_train_CV_death,Y,death_currTrainIdx,4,temp_classifier_choice)
    } else {
      clin_only_death_models <- train_caret_models(clin_only_train_CV_death,Y,death_currTrainIdx,4,classifier_choice)
    }
    
    mf_only_death_models <- train_caret_models(mf_only_train_death_CV,Y,death_currTrainIdx,4,classifier_choice)
    
    # Predict outcomes on testing data using trained models
    
    Y <- as.factor(patient_clinical_data_yr$Death12Months)
    
    death_test_val<-predict_caret_models(death_models, test_death_CV, Y, death_currTestIdx)
    
    clin_only_death_test_val<-predict_caret_models(clin_only_death_models, clin_only_test_CV_death, Y, death_currTestIdx)
    
    mf_only_death_test_val<-predict_caret_models(mf_only_death_models, mf_only_test_death_CV, Y, death_currTestIdx)
    
    death_predictions[[i]] <- death_test_val
    clin_only_death_predictions[[i]] <- clin_only_death_test_val
    mf_only_death_predictions[[i]] <- mf_only_death_test_val
  }
  
  outList <-
    list(
      death_predictions,
      clin_only_death_predictions,
      mf_only_death_predictions
    )
  names(outList) <-
    c(
      "death_predictions",
      "clin_only_death_predictions",
      "mf_only_death_predictions"
    )
  return(outList)
}