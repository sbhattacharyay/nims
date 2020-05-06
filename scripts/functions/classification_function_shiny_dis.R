# Load patient clinical data (sorts by PY numbering and corrects variable types)
patient_clinical_data <- load_patient_clinical_data()

# Load and update clinical variable list
clinicalVariableList <- update_clinicalVariableList()

# Load Motion Features Organized by Time of Day (TOD)
if (!exists("tod_sensors")) {
  tod_sensors <-
    readMat('../motion_feature_data/bed_corrected_imputed_complete_sensor_data.mat')$bed.corrected.sensors
}

# Load Motion Features Organized by Time from Recording (TFR)
if (!exists("tfr_sensors")) {
  tfr_sensors <-
    readMat('../tfr_motion_feature_data/bed_corrected_imputed_complete_sensor_data.mat')$bed.corrected.sensors
}

# Organizing Features for both motion feature categorizations
tod_motion_features <- get_motion_features(do.call(rbind, tod_sensors))
tfr_motion_features <- get_motion_features(do.call(rbind, tfr_sensors))

# Load transformed covariates for both TOD and TFR
tod_tf_covariates <-load_tf_patient_covariates('../motion_feature_data/tf_patient_covariates.csv')
tfr_tf_covariates <-load_tf_patient_covariates('../tfr_motion_feature_data/tf_patient_covariates.csv')

tod_sensor_idx_ranges <- split(1:77754, ceiling(seq_along(1:77754)/12959))
tod_fp2_idx_ranges    <- lapply(tod_sensor_idx_ranges, function(x) x+77754)

tfr_sensor_idx_ranges <- split(1:34554, ceiling(seq_along(1:34554)/5759))
tfr_fp2_idx_ranges    <- lapply(tfr_sensor_idx_ranges, function(x) x+34554)

tfr_range<-seq(5/3600,8,length.out=5759)
tod_range<-seq(as.POSIXct("2020-05-05 18:00:05"),as.POSIXct("2020-05-06 12:00:00"),length.out=12959)

sensor_options <- c("left_ank","left_el","left_wr","right_ank","right_el","right_wr")

classification_function_shiny_dis <-function(time_choice,time_slide,classifier_choice,r,mf_choice,clinicalVars,sensor_loc){
  
  sensor_choice_Idxs <- match(sensor_loc,sensor_options)
  
  if (time_choice == 'tod'){
    sensor_idx_ranges<-lapply(tod_sensor_idx_ranges,function(x) x[tod_range >= time_slide[1] & tod_range <= time_slide[2]])
    fp2_idx_ranges   <-lapply(tod_fp2_idx_ranges,function(x) x[tod_range >= time_slide[1] & tod_range <= time_slide[2]])
    curr_sensor_Idx<-do.call(c, sensor_idx_ranges[sensor_choice_Idxs])
    curr_fp_Idx<-c(curr_sensor_Idx,do.call(c, fp2_idx_ranges[sensor_choice_Idxs]))
    curr_motion_features <- lapply(tod_motion_features[mf_choice], function(x) x[,curr_sensor_Idx])
    
    if ('freq_pairsFeats' %in% mf_choice) curr_motion_features$freq_pairsFeats <- tod_motion_features$freq_pairsFeats[,curr_fp_Idx]
    
    curr_motion_features <- do.call(cbind,curr_motion_features)
    tf_covariates <- tod_tf_covariates
    
  } else {
    sensor_idx_ranges<-lapply(tfr_sensor_idx_ranges,function(x) x[tfr_range >= time_slide[1] & tfr_range <= time_slide[2]])
    fp2_idx_ranges   <-lapply(tfr_fp2_idx_ranges,function(x) x[tfr_range >= time_slide[1] & tfr_range <= time_slide[2]])
    curr_sensor_Idx<-do.call(c, sensor_idx_ranges[sensor_choice_Idxs])
    curr_fp_Idx<-c(curr_sensor_Idx,do.call(c, fp2_idx_ranges[sensor_choice_Idxs]))
    curr_motion_features <- lapply(tfr_motion_features[mf_choice], function(x) x[,curr_sensor_Idx])
    
    if ('freq_pairsFeats' %in% mf_choice) curr_motion_features$freq_pairsFeats <- tfr_motion_features$freq_pairsFeats[,curr_fp_Idx]
    
    curr_motion_features <- do.call(cbind,curr_motion_features)
    tf_covariates <- tfr_tf_covariates
    
  }  
  
  # Split for k-fold cross validation for discharge predictions:
  k <- 5
  cvIdx <- cross_val_splits(patient_clinical_data, k)
  
  GOSE_predictions  <- vector(mode = "list", length = k)
  fav_predictions   <- vector(mode = "list", length = k)
  death_predictions <- vector(mode = "list", length = k)
  
  for (i in 1:length(cvIdx)) {
    currTestIdx <- cvIdx[[as.character(i)]]
    currTrainIdx <- seq(nrow(patient_clinical_data))[-currTestIdx]
    
    # Convert numeric to logical indexing
    logicalTest <- rep(FALSE, nrow(patient_clinical_data))
    logicalTest[currTestIdx] <- TRUE
    logicalTrain <- !logicalTest
    
    train_tf_covariates <- tf_covariates[currTrainIdx,]
    test_tf_covariates  <- tf_covariates[currTestIdx,]
    
    # Perform LOL on training data
    
    # - GOSE LOL:
    Y <- as.factor(patient_clinical_data$gose)
    logicalIdx <- logicalTrain & !is.na(Y)
    GOSE_LOL<-lol.project.lol(curr_motion_features[logicalIdx,], Y[logicalIdx], r)
    
    # - Favorable Outcome LOL:
    Y <- as.factor(patient_clinical_data$favorable)
    logicalIdx <- logicalTrain & !is.na(Y)
    fav_LOL<-lol.project.lol(curr_motion_features[logicalIdx,], Y[logicalIdx], r)

    # - Mortality Outcome LOL:
    Y <- as.factor(patient_clinical_data$death)
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
    Y <- as.factor(patient_clinical_data$favorable)
    
    GOSE_models <- train_caret_models(train_GOSE_CV,Y,currTrainIdx,4,classifier_choice)
    fav_models  <- train_caret_models(train_fav_CV,Y,currTrainIdx,4,classifier_choice)

    Y <- as.factor(patient_clinical_data$death)
    
    death_models  <- train_caret_models(train_death_CV,Y,currTrainIdx,4,classifier_choice)

    # Predict outcomes on testing data using trained models

    Y <- as.factor(patient_clinical_data$favorable)
    
    GOSE_test_val <-predict_caret_models(GOSE_models, test_GOSE_CV, Y, currTestIdx)
    fav_test_val  <-predict_caret_models(fav_models, test_fav_CV, Y, currTestIdx)
   
    Y <- as.factor(patient_clinical_data$death)
    
    death_test_val<-predict_caret_models(death_models, test_death_CV, Y, currTestIdx)
    
    GOSE_predictions[[i]] <- GOSE_test_val
    fav_predictions[[i]] <- fav_test_val
    death_predictions[[i]] <- death_test_val
  }
  
  outList <- list(GOSE_predictions,fav_predictions,death_predictions)
  names(outList) <- c("GOSE_predictions","fav_predictions","death_predictions")
  return(outList)
}