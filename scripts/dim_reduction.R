#### LOL Embedding for class imbalance ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# University of Cambridge
# Johns Hopkins University
# email address: sb2406@cam.ac.uk

# Add library path for R packages and set working directory:
.libPaths(c("~/Rpackages" , .libPaths()))
setwd('~/work/nims/scripts')

# Load necessary packages
library(tidyverse)
library(readxl)
library(caret)
library(lolR)
library(UBL)
library(viridis)
library(shadowtext)

# Load clinical patient data:
source('./functions/load_patient_clinical_data.R')
patient_clinical_data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv') %>% arrange(AccelPatientNo_) %>% mutate(ptIdx = 1:nrow(.))

### Format feature spaces into matrices and save matrices

dir.create("../all_motion_feature_data/04_formatted_predictor_matrices",showWarnings = FALSE)
impFiles <- list.files(path="../all_motion_feature_data/03_bed_corrected_imputed_features/",pattern = "*.RData")

# (a) Detection model matrices

# Load detection labels and partitions
load('../validation_resampling/detection_labels.RData')
load('../validation_resampling/detection_partitions.RData')

for (impNo in 1:length(impFiles)){
  print(paste("Imputation no.",impNo,"out of",length(impFiles),"started."))
  dir.create(paste0("../all_motion_feature_data/04_formatted_predictor_matrices/imp",impNo),showWarnings = FALSE)
  load(paste0('../all_motion_feature_data/03_bed_corrected_imputed_features/',impFiles[impNo]))
  gc()
  for (i in 1:length(det_gcs_labels)) {
    print(paste("Parameter combination no.",i,"out of",length(det_gcs_labels),"started."))
    curr_obs_window <- det_parameters$obs_windows[i]
    curr_secs_window <- curr_obs_window*3600
    
    curr_det_gcs_labels <- det_gcs_labels[[i]]
    motor_indices <- which(!is.na(curr_det_gcs_labels$Best.Motor.Response))
    eye_indices <- which(!is.na(curr_det_gcs_labels$Eye.Opening))
    curr_det_motor_train_idx <- det_motor_train_idx[[i]]
    curr_det_eye_train_idx <- det_eye_train_idx[[i]]
    curr_det_motor_test_idx <- det_motor_test_idx[[i]]
    curr_det_eye_test_idx <- det_eye_test_idx[[i]]
    dir.create(paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/detection_window_',curr_obs_window),showWarnings = FALSE)
    
    # Motor train matrix:
    curr_motor_train_matrix <- as.data.frame(matrix(nrow = 0, ncol = ((curr_secs_window/5)+1)*7*6))
    curr.long.motor.train.df <- as.data.frame(matrix(ncol=8, nrow = 0))
    print("Motor training matrix started:")
    for (j in 1:length(curr_det_motor_train_idx)){
      curr_GCS_row <- curr_det_motor_train_idx[j]
      currTimeStamp <- curr_det_gcs_labels$TakenInstant[curr_GCS_row]
      curr_GCSm <- curr_det_gcs_labels$Best.Motor.Response[curr_GCS_row]
      currFeatSet <- currImp %>% 
        filter(patient_clinical_data$AccelPatientNo_[ptIdx] == curr_det_gcs_labels$AccelPatientNo_[curr_GCS_row], 
               TakenInstant <= currTimeStamp, 
               TakenInstant >= currTimeStamp - curr_secs_window)
      long.form.curr.feat.set <- currFeatSet %>% 
        pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>%
        select(-c(timeStamps)) 
      long.form.curr.feat.set <- long.form.curr.feat.set %>%
        mutate(timeCount = plyr::mapvalues(timeCount, 
                                           from = unique(long.form.curr.feat.set$timeCount), 
                                           to = 1:length(unique(long.form.curr.feat.set$timeCount))),
               GCSm = curr_GCSm,
               sampleIdx = j)
      curr.long.motor.train.df <- rbind(curr.long.motor.train.df,long.form.curr.feat.set)
      tempReshape <- long.form.curr.feat.set %>% 
        pivot_wider(names_from=c("featureType","sensor"),values_from="value") %>% 
        select(-c(ptIdx,TakenInstant,GCSm,sampleIdx)) %>%
        pivot_longer(cols = -timeCount) %>%
        pivot_wider(names_from = c(name,timeCount),values_from = value)
      curr_motor_train_matrix <- rbind(curr_motor_train_matrix,tempReshape)
      print(paste("Row no.",j,"out of",length(curr_det_motor_train_idx),"complete."))
    }
    saveRDS(curr.long.motor.train.df,file = paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/detection_window_',curr_obs_window,'/long_form_motor_train_values.rds'))
    saveRDS(curr_motor_train_matrix,file = paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/detection_window_',curr_obs_window,'/motor_train_matrix.rds'))
    
    # Eye train matrix:
    curr_eye_train_matrix <- as.data.frame(matrix(nrow = 0, ncol = ((curr_secs_window/5)+1)*7*6))
    curr.long.eye.train.df <- as.data.frame(matrix(ncol=8, nrow = 0))
    print("Eye training matrix started:")
    for (j in 1:length(curr_det_eye_train_idx)){
      curr_GCS_row <- curr_det_eye_train_idx[j]
      currTimeStamp <- curr_det_gcs_labels$TakenInstant[curr_GCS_row]
      curr_GCSe <- curr_det_gcs_labels$Eye.Opening[curr_GCS_row]
      currFeatSet <- currImp %>% 
        filter(patient_clinical_data$AccelPatientNo_[ptIdx] == curr_det_gcs_labels$AccelPatientNo_[curr_GCS_row], 
               TakenInstant <= currTimeStamp, 
               TakenInstant >= currTimeStamp - curr_secs_window)
      long.form.curr.feat.set <- currFeatSet %>% 
        pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>%
        select(-c(timeStamps)) 
      long.form.curr.feat.set <- long.form.curr.feat.set %>%
        mutate(timeCount = plyr::mapvalues(timeCount, 
                                           from = unique(long.form.curr.feat.set$timeCount), 
                                           to = 1:length(unique(long.form.curr.feat.set$timeCount))),
               GCSe = curr_GCSe,
               sampleIdx = j)
      curr.long.eye.train.df <- rbind(curr.long.eye.train.df,long.form.curr.feat.set)
      tempReshape <- long.form.curr.feat.set %>% 
        pivot_wider(names_from=c("featureType","sensor"),values_from="value") %>% 
        select(-c(ptIdx,TakenInstant,GCSe,sampleIdx)) %>%
        pivot_longer(cols = -timeCount) %>%
        pivot_wider(names_from = c(name,timeCount),values_from = value)
      curr_eye_train_matrix <- rbind(curr_eye_train_matrix,tempReshape)
      print(paste("Row no.",j,"out of",length(curr_det_eye_train_idx),"complete."))
    }
    saveRDS(curr.long.eye.train.df,file = paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/detection_window_',curr_obs_window,'/long_form_eye_train_values.rds'))
    saveRDS(curr_eye_train_matrix,file = paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/detection_window_',curr_obs_window,'/eye_train_matrix.rds'))

    # Motor test matrix:
    curr_motor_test_matrix <- as.data.frame(matrix(nrow = 0, ncol = ((curr_secs_window/5)+1)*7*6))
    curr.long.motor.test.df <- as.data.frame(matrix(ncol=8, nrow = 0))
    print("Motor testing matrix started:")
    for (j in 1:length(curr_det_motor_test_idx)){
      curr_GCS_row <- curr_det_motor_test_idx[j]
      currTimeStamp <- curr_det_gcs_labels$TakenInstant[curr_GCS_row]
      curr_GCSm <- curr_det_gcs_labels$Best.Motor.Response[curr_GCS_row]
      currFeatSet <- currImp %>% 
        filter(patient_clinical_data$AccelPatientNo_[ptIdx] == curr_det_gcs_labels$AccelPatientNo_[curr_GCS_row], 
               TakenInstant <= currTimeStamp, 
               TakenInstant >= currTimeStamp - curr_secs_window)
      long.form.curr.feat.set <- currFeatSet %>% 
        pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>%
        select(-c(timeStamps)) 
      long.form.curr.feat.set <- long.form.curr.feat.set %>%
        mutate(timeCount = plyr::mapvalues(timeCount, 
                                           from = unique(long.form.curr.feat.set$timeCount), 
                                           to = 1:length(unique(long.form.curr.feat.set$timeCount))),
               GCSm = curr_GCSm,
               sampleIdx = j)
      curr.long.motor.test.df <- rbind(curr.long.motor.test.df,long.form.curr.feat.set)
      tempReshape <- long.form.curr.feat.set %>% 
        pivot_wider(names_from=c("featureType","sensor"),values_from="value") %>% 
        select(-c(ptIdx,TakenInstant,GCSm,sampleIdx)) %>%
        pivot_longer(cols = -timeCount) %>%
        pivot_wider(names_from = c(name,timeCount),values_from = value)
      curr_motor_test_matrix <- rbind(curr_motor_test_matrix,tempReshape)
      print(paste("Row no.",j,"out of",length(curr_det_motor_test_idx),"complete."))
    }
    saveRDS(curr.long.motor.test.df,file = paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/detection_window_',curr_obs_window,'/long_form_motor_test_values.rds'))
    saveRDS(curr_motor_test_matrix,file = paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/detection_window_',curr_obs_window,'/motor_test_matrix.rds'))
    
    # Eye test matrix:
    curr_eye_test_matrix <- as.data.frame(matrix(nrow = 0, ncol = ((curr_secs_window/5)+1)*7*6))
    curr.long.eye.test.df <- as.data.frame(matrix(ncol=8, nrow = 0))
    print("Eye testing matrix started:")
    for (j in 1:length(curr_det_eye_test_idx)){
      curr_GCS_row <- curr_det_eye_test_idx[j]
      currTimeStamp <- curr_det_gcs_labels$TakenInstant[curr_GCS_row]
      curr_GCSe <- curr_det_gcs_labels$Eye.Opening[curr_GCS_row]
      currFeatSet <- currImp %>% 
        filter(patient_clinical_data$AccelPatientNo_[ptIdx] == curr_det_gcs_labels$AccelPatientNo_[curr_GCS_row], 
               TakenInstant <= currTimeStamp, 
               TakenInstant >= currTimeStamp - curr_secs_window)
      long.form.curr.feat.set <- currFeatSet %>% 
        pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>%
        select(-c(timeStamps)) 
      long.form.curr.feat.set <- long.form.curr.feat.set %>%
        mutate(timeCount = plyr::mapvalues(timeCount, 
                                           from = unique(long.form.curr.feat.set$timeCount), 
                                           to = 1:length(unique(long.form.curr.feat.set$timeCount))),
               GCSe = curr_GCSe,
               sampleIdx = j)
      curr.long.eye.test.df <- rbind(curr.long.eye.test.df,long.form.curr.feat.set)
      tempReshape <- long.form.curr.feat.set %>% 
        pivot_wider(names_from=c("featureType","sensor"),values_from="value") %>% 
        select(-c(ptIdx,TakenInstant,GCSe,sampleIdx)) %>%
        pivot_longer(cols = -timeCount) %>%
        pivot_wider(names_from = c(name,timeCount),values_from = value)
      curr_eye_test_matrix <- rbind(curr_eye_test_matrix,tempReshape)
      print(paste("Row no.",j,"out of",length(curr_det_eye_test_idx),"complete."))
    }
    saveRDS(curr.long.eye.test.df,file = paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/detection_window_',curr_obs_window,'/long_form_eye_test_values.rds'))
    saveRDS(curr_eye_test_matrix,file = paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/detection_window_',curr_obs_window,'/eye_test_matrix.rds'))
    
    print(paste("Parameter combination no.",i,"out of",length(det_gcs_labels),"completed."))
  }
  print(paste("Imputation no.",impNo,"out of",length(impFiles),"completed."))
}

# temp: correct GCSm train and test matrices
impDirs <- list.files('../all_motion_feature_data/04_formatted_predictor_matrices',pattern = 'imp*',include.dirs = TRUE, full.names = TRUE)
for (i in 1:length(impDirs)){
  print(paste("Imputation no.",i,"out of",length(impDirs),"started."))
  currImpDir <- impDirs[i]
  detection_folders <- list.files(path = currImpDir,pattern = 'detection_*',include.dirs = TRUE, full.names = TRUE)
  for (j in 1:length(detection_folders)){
    print(paste("Detection folder no.",j,"out of",length(detection_folders),"started."))
    curr.motor.train.long.form <- readRDS(file.path(detection_folders[j],'long_form_motor_train_values.rds'))
    curr.motor.test.long.form <- readRDS(file.path(detection_folders[j],'long_form_motor_test_values.rds'))
    train.matrix <- curr.motor.train.long.form %>% 
      select(-c(ptIdx,TakenInstant,GCSm)) %>%
      pivot_wider(id_cols = sampleIdx,names_from=c("featureType","sensor","timeCount"),values_from="value") %>%
      select(-sampleIdx)
    test.matrix <- curr.motor.test.long.form %>% 
      select(-c(ptIdx,TakenInstant,GCSm)) %>%
      pivot_wider(id_cols = sampleIdx,names_from=c("featureType","sensor","timeCount"),values_from="value") %>%
      select(-sampleIdx)
    saveRDS(train.matrix,file = file.path(detection_folders[j],'motor_train_matrix.rds'))
    saveRDS(test.matrix,file = file.path(detection_folders[j],'motor_test_matrix.rds'))
    rm(train.matrix,test.matrix)
    gc()
  }
}


# (b) Prediction model matrices

# Load prediction labels and partitions
load('../validation_resampling/prediction_labels.RData')
load('../validation_resampling/prediction_partitions.RData')

for (impNo in 1:length(impFiles)){
  print(paste("Imputation no.",impNo,"out of",length(impFiles),"started."))
  dir.create(paste0("../all_motion_feature_data/04_formatted_predictor_matrices/imp",impNo),showWarnings = FALSE)
  load(paste0('../all_motion_feature_data/03_bed_corrected_imputed_features/',impFiles[impNo]))
  gc()
  for (i in 1:length(pre_gcs_labels)) {
    print(paste("Parameter combination no.",i,"out of",length(pre_gcs_labels),"started."))
    curr_obs_window <- pre_parameters$obs_windows[i]
    curr_lead_time <- pre_parameters$lead_times[i]
    curr_secs_window <- curr_obs_window*3600
    curr_secs_lead <- curr_lead_time*3600
    
    curr_pre_gcs_labels <- pre_gcs_labels[[i]]
    motor_indices <- which(!is.na(curr_pre_gcs_labels$Best.Motor.Response))
    eye_indices <- which(!is.na(curr_pre_gcs_labels$Eye.Opening))
    curr_pre_motor_train_idx <- pre_motor_train_idx[[i]]
    curr_pre_eye_train_idx <- pre_eye_train_idx[[i]]
    curr_pre_motor_test_idx <- pre_motor_test_idx[[i]]
    curr_pre_eye_test_idx <- pre_eye_test_idx[[i]]
    
    dir.create(paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/prediction_window_',curr_obs_window,'_lead_',curr_lead_time),showWarnings = FALSE)
    
    # Motor train matrix:
    curr_motor_train_matrix <- as.data.frame(matrix(nrow = 0, ncol = ((curr_secs_window/5)+1)*7*6))
    curr.long.motor.train.df <- as.data.frame(matrix(ncol=8, nrow = 0))
    print("Motor training matrix started:")
    for (j in 1:length(curr_pre_motor_train_idx)){
      curr_GCS_row <- curr_pre_motor_train_idx[j]
      currTimeStamp <- curr_pre_gcs_labels$TakenInstant[curr_GCS_row]
      curr_label <- curr_pre_gcs_labels$Net.GCSm.Change[curr_GCS_row]
      currFeatSet <- currImp %>% 
        filter(patient_clinical_data$AccelPatientNo_[ptIdx] == curr_pre_gcs_labels$AccelPatientNo_[curr_GCS_row], 
               TakenInstant <= currTimeStamp, 
               TakenInstant >= currTimeStamp - curr_secs_window)
      long.form.curr.feat.set <- currFeatSet %>% 
        pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>%
        select(-c(timeStamps)) 
      long.form.curr.feat.set <- long.form.curr.feat.set %>%
        mutate(timeCount = plyr::mapvalues(timeCount, 
                                           from = unique(long.form.curr.feat.set$timeCount), 
                                           to = 1:length(unique(long.form.curr.feat.set$timeCount))),
               label = curr_label,
               sampleIdx = j)
      curr.long.motor.train.df <- rbind(curr.long.motor.train.df,long.form.curr.feat.set)
      tempReshape <- long.form.curr.feat.set %>% 
        pivot_wider(names_from=c("featureType","sensor"),values_from="value") %>% 
        select(-c(ptIdx,TakenInstant,GCSe,sampleIdx)) %>%
        pivot_longer(cols = -timeCount) %>%
        pivot_wider(names_from = c(name,timeCount),values_from = value)
      curr_motor_train_matrix <- rbind(curr_motor_train_matrix,tempReshape)
      print(paste("Row no.",j,"out of",length(curr_pre_motor_train_idx),"complete."))
    }
    saveRDS(curr.long.motor.train.df,file = paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/prediction_window_',curr_obs_window,'/long_form_motor_train_values.rds'))
    saveRDS(curr_motor_train_matrix,file = paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/prediction_window_',curr_obs_window,'/motor_train_matrix.rds'))
    
    # Eye train matrix:
    curr_eye_train_matrix <- as.data.frame(matrix(nrow = 0, ncol = ((curr_secs_window/5)+1)*7*6))
    curr.long.eye.train.df <- as.data.frame(matrix(ncol=8, nrow = 0))
    print("Eye training matrix started:")
    for (j in 1:length(curr_pre_eye_train_idx)){
      curr_GCS_row <- curr_pre_eye_train_idx[j]
      currTimeStamp <- curr_pre_gcs_labels$TakenInstant[curr_GCS_row]
      curr_label <- curr_pre_gcs_labels$Net.GCSe.Change[curr_GCS_row]
      currFeatSet <- currImp %>% 
        filter(patient_clinical_data$AccelPatientNo_[ptIdx] == curr_pre_gcs_labels$AccelPatientNo_[curr_GCS_row], 
               TakenInstant <= currTimeStamp, 
               TakenInstant >= currTimeStamp - curr_secs_window)
      long.form.curr.feat.set <- currFeatSet %>% 
        pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>%
        select(-c(timeStamps)) 
      long.form.curr.feat.set <- long.form.curr.feat.set %>%
        mutate(timeCount = plyr::mapvalues(timeCount, 
                                           from = unique(long.form.curr.feat.set$timeCount), 
                                           to = 1:length(unique(long.form.curr.feat.set$timeCount))),
               label = curr_label,
               sampleIdx = j)
      curr.long.eye.train.df <- rbind(curr.long.eye.train.df,long.form.curr.feat.set)
      tempReshape <- long.form.curr.feat.set %>% 
        pivot_wider(names_from=c("featureType","sensor"),values_from="value") %>% 
        select(-c(ptIdx,TakenInstant,GCSe,sampleIdx)) %>%
        pivot_longer(cols = -timeCount) %>%
        pivot_wider(names_from = c(name,timeCount),values_from = value)
      curr_eye_train_matrix <- rbind(curr_eye_train_matrix,tempReshape)
      print(paste("Row no.",j,"out of",length(curr_pre_eye_train_idx),"complete."))
    }
    saveRDS(curr.long.eye.train.df,file = paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/prediction_window_',curr_obs_window,'/long_form_eye_train_values.rds'))
    saveRDS(curr_eye_train_matrix,file = paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/prediction_window_',curr_obs_window,'/eye_train_matrix.rds'))
    
    # Motor test matrix:
    curr_motor_test_matrix <- as.data.frame(matrix(nrow = 0, ncol = ((curr_secs_window/5)+1)*7*6))
    curr.long.motor.test.df <- as.data.frame(matrix(ncol=8, nrow = 0))
    print("Motor testing matrix started:")
    for (j in 1:length(curr_pre_motor_test_idx)){
      curr_GCS_row <- curr_pre_motor_test_idx[j]
      currTimeStamp <- curr_pre_gcs_labels$TakenInstant[curr_GCS_row]
      curr_label <- curr_pre_gcs_labels$Net.GCSm.Change[curr_GCS_row]
      currFeatSet <- currImp %>% 
        filter(patient_clinical_data$AccelPatientNo_[ptIdx] == curr_pre_gcs_labels$AccelPatientNo_[curr_GCS_row], 
               TakenInstant <= currTimeStamp, 
               TakenInstant >= currTimeStamp - curr_secs_window)
      long.form.curr.feat.set <- currFeatSet %>% 
        pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>%
        select(-c(timeStamps)) 
      long.form.curr.feat.set <- long.form.curr.feat.set %>%
        mutate(timeCount = plyr::mapvalues(timeCount, 
                                           from = unique(long.form.curr.feat.set$timeCount), 
                                           to = 1:length(unique(long.form.curr.feat.set$timeCount))),
               label = curr_label,
               sampleIdx = j)
      curr.long.motor.test.df <- rbind(curr.long.motor.test.df,long.form.curr.feat.set)
      tempReshape <- long.form.curr.feat.set %>% 
        pivot_wider(names_from=c("featureType","sensor"),values_from="value") %>% 
        select(-c(ptIdx,TakenInstant,GCSe,sampleIdx)) %>%
        pivot_longer(cols = -timeCount) %>%
        pivot_wider(names_from = c(name,timeCount),values_from = value)
      curr_motor_test_matrix <- rbind(curr_motor_test_matrix,tempReshape)
      print(paste("Row no.",j,"out of",length(curr_pre_motor_test_idx),"complete."))
    }
    saveRDS(curr.long.motor.test.df,file = paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/prediction_window_',curr_obs_window,'/long_form_motor_test_values.rds'))
    saveRDS(curr_motor_test_matrix,file = paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/prediction_window_',curr_obs_window,'/motor_test_matrix.rds'))
    
    # Eye test matrix:
    curr_eye_test_matrix <- as.data.frame(matrix(nrow = 0, ncol = ((curr_secs_window/5)+1)*7*6))
    curr.long.eye.test.df <- as.data.frame(matrix(ncol=8, nrow = 0))
    print("Eye testing matrix started:")
    for (j in 1:length(curr_pre_eye_test_idx)){
      curr_GCS_row <- curr_pre_eye_test_idx[j]
      currTimeStamp <- curr_pre_gcs_labels$TakenInstant[curr_GCS_row]
      curr_label <- curr_pre_gcs_labels$Net.GCSe.Change[curr_GCS_row]
      currFeatSet <- currImp %>% 
        filter(patient_clinical_data$AccelPatientNo_[ptIdx] == curr_pre_gcs_labels$AccelPatientNo_[curr_GCS_row], 
               TakenInstant <= currTimeStamp, 
               TakenInstant >= currTimeStamp - curr_secs_window)
      long.form.curr.feat.set <- currFeatSet %>% 
        pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>%
        select(-c(timeStamps)) 
      long.form.curr.feat.set <- long.form.curr.feat.set %>%
        mutate(timeCount = plyr::mapvalues(timeCount, 
                                           from = unique(long.form.curr.feat.set$timeCount), 
                                           to = 1:length(unique(long.form.curr.feat.set$timeCount))),
               label = curr_label,
               sampleIdx = j)
      curr.long.eye.test.df <- rbind(curr.long.eye.test.df,long.form.curr.feat.set)
      tempReshape <- long.form.curr.feat.set %>% 
        pivot_wider(names_from=c("featureType","sensor"),values_from="value") %>% 
        select(-c(ptIdx,TakenInstant,GCSe,sampleIdx)) %>%
        pivot_longer(cols = -timeCount) %>%
        pivot_wider(names_from = c(name,timeCount),values_from = value)
      curr_eye_test_matrix <- rbind(curr_eye_test_matrix,tempReshape)
      print(paste("Row no.",j,"out of",length(curr_pre_eye_test_idx),"complete."))
    }
    saveRDS(curr.long.eye.test.df,file = paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/prediction_window_',curr_obs_window,'/long_form_eye_test_values.rds'))
    saveRDS(curr_eye_test_matrix,file = paste0('../all_motion_feature_data/04_formatted_predictor_matrices/imp',impNo,'/prediction_window_',curr_obs_window,'/eye_test_matrix.rds'))
    
    print(paste("Parameter combination no.",i,"out of",length(pre_gcs_labels),"completed"))
  }
  print(paste("Imputation no.",impNo,"out of",length(impFiles),"completed."))
}

### Train and perform linear optimal low-rank projection:
rm(list = ls())
gc()

impDirs <- list.files('../all_motion_feature_data/04_formatted_predictor_matrices',pattern = 'imp*',include.dirs = TRUE, full.names = TRUE)

# (a) LOL on Detection Cases

# Load detection labels and partitions
load('../validation_resampling/detection_labels.RData')
load('../validation_resampling/detection_partitions.RData')

for (i in 1:length(impDirs)){
  print(paste("Imputation no.",i,"out of",length(impDirs),"started."))
  currImpDir <- impDirs[i]
  detection_folders <- list.files(path = currImpDir,pattern = 'detection_*',include.dirs = TRUE, full.names = TRUE)
  for (j in 1:length(detection_folders)){
    
    print(paste("Detection folder no.",j,"out of",length(detection_folders),"started."))
    
    curr_window_size <- as.numeric(sub(".*window_", "", detection_folders[j]))
    curr_window_idx <- which(det_parameters$obs_windows == curr_window_size)
    
    curr_label_set <- det_gcs_labels[[curr_window_idx]]
    
    tryCatch({
      curr_motor_train_labels <- curr_label_set$Best.Motor.Response[det_motor_train_idx[[curr_window_idx]]]
      
      curr_motor_train_matrix <- as.matrix(readRDS(file.path(detection_folders[j],'motor_train_matrix.rds')) %>% scale())
      curr_motor_test_matrix <- as.matrix(readRDS(file.path(detection_folders[j],'motor_test_matrix.rds')) %>% scale())
      
      curr_motor_lol <- lol.project.lol(curr_motor_train_matrix,curr_motor_train_labels,r = 200)
      
      saveRDS(curr_motor_lol,file.path(detection_folders[j],'motor_lol2.rds'))
      saveRDS(curr_motor_lol$Xr,file.path(detection_folders[j],'motor_train_matrix_lol2.rds'))
      saveRDS(curr_motor_test_matrix %*% curr_motor_lol$A,file.path(detection_folders[j],'motor_test_matrix_lol2.rds'))
      
      rm(curr_motor_train_matrix,curr_motor_test_matrix,curr_motor_lol)
      gc()
    }, error=function(e){cat("ERROR ON MOTOR, WINDOW SIZE", curr_window_size,"HOURS:",conditionMessage(e), "\n")})

    # tryCatch({
    #   curr_eye_train_labels <- curr_label_set$Eye.Opening[det_eye_train_idx[[curr_window_idx]]]
    #   
    #   curr_eye_train_matrix <- readRDS(file.path(detection_folders[j],'eye_train_matrix.rds'))
    #   curr_eye_test_matrix <- readRDS(file.path(detection_folders[j],'eye_test_matrix.rds'))
    #   
    #   curr_eye_lol <- lol.project.lol(abs(curr_eye_train_matrix),curr_eye_train_labels,r = 200)
    #   
    #   saveRDS(curr_eye_lol,file.path(detection_folders[j],'eye_lol.rds'))
    #   saveRDS(curr_eye_lol$Xr,file.path(detection_folders[j],'eye_train_matrix_lol.rds'))
    #   saveRDS(curr_eye_test_matrix %*% curr_eye_lol$A,file.path(detection_folders[j],'eye_test_matrix_lol.rds'))
    #   
    # }, error=function(e){cat("ERROR ON EYE, WINDOW SIZE", curr_window_size,"HOURS:",conditionMessage(e), "\n")})
    # 
    print(paste("Detection folder no.",j,"out of",length(detection_folders),"completed."))
  }
  print(paste("Imputation no.",i,"out of",length(impDirs),"completed."))
}

# (b) LOL on Prediction Cases

# Load prediction labels and partitions
load('../validation_resampling/prediction_labels.RData')
load('../validation_resampling/prediction_partitions.RData')

for (i in 1:length(impDirs)){
  print(paste("Imputation no.",i,"out of",length(impDirs),"started."))
  currImpDir <- impDirs[i]
  prediction_folders <- list.files(path = currImpDir,pattern = 'prediction_*',include.dirs = TRUE, full.names = TRUE)
  for (j in 1:length(prediction_folders)){
    
    gc()
    
    print(paste("Prediction folder no.",j,"out of",length(prediction_folders),"started."))
    
    pattern <- "window_\\s*(.*?)\\s*_lead"
    curr_window_size <- as.numeric(regmatches(prediction_folders[j], regexec(pattern, prediction_folders[j]))[[1]][2])
    curr_lead_time <- as.numeric(sub(".*lead_", "", prediction_folders[j]))
      
    curr_window_idx <- which(pre_parameters$obs_windows == curr_window_size & pre_parameters$lead_times == curr_lead_time)
    
    curr_label_set <- pre_gcs_labels[[curr_window_idx]]
    
    tryCatch({
      curr_motor_train_labels <- curr_label_set$Net.GCSm.Change[pre_motor_train_idx[[curr_window_idx]]]
      
      curr_motor_train_matrix <- readRDS(file.path(prediction_folders[j],'motor_train_matrix.rds'))
      curr_motor_test_matrix <- readRDS(file.path(prediction_folders[j],'motor_test_matrix.rds'))
      
      curr_motor_lol <- lol.project.lol(abs(curr_motor_train_matrix),curr_motor_train_labels,r = 200)
      
      saveRDS(curr_motor_lol,file.path(prediction_folders[j],'motor_lol.rds'))
      saveRDS(curr_motor_lol$Xr,file.path(prediction_folders[j],'motor_train_matrix_lol.rds'))
      saveRDS(curr_motor_test_matrix %*% curr_motor_lol$A,file.path(prediction_folders[j],'motor_test_matrix_lol.rds'))
      
    }, error=function(e){cat("ERROR ON MOTOR, WINDOW SIZE", curr_window_size,"HOURS:",conditionMessage(e), "\n")})
    
    tryCatch({
      curr_eye_train_labels <- curr_label_set$Net.GCSe.Change[pre_eye_train_idx[[curr_window_idx]]]
      
      curr_eye_train_matrix <- readRDS(file.path(prediction_folders[j],'eye_train_matrix.rds'))
      curr_eye_test_matrix <- readRDS(file.path(prediction_folders[j],'eye_test_matrix.rds'))
      
      curr_eye_lol <- lol.project.lol(abs(curr_eye_train_matrix),curr_eye_train_labels,r = 200)
      
      saveRDS(curr_eye_lol,file.path(prediction_folders[j],'eye_lol.rds'))
      saveRDS(curr_eye_lol$Xr,file.path(prediction_folders[j],'eye_train_matrix_lol.rds'))
      saveRDS(curr_eye_test_matrix %*% curr_eye_lol$A,file.path(prediction_folders[j],'eye_test_matrix_lol.rds'))
      
    }, error=function(e){cat("ERROR ON EYE, WINDOW SIZE", curr_window_size,"HOURS:",conditionMessage(e), "\n")})
    
    print(paste("Prediction folder no.",j,"out of",length(prediction_folders),"completed."))
  }
  print(paste("Imputation no.",i,"out of",length(impDirs),"completed."))
}

#### Perform SMOTE to repair class imbalances
# rm(list = ls())
# gc()

# impDirs <- list.files('../all_motion_feature_data/04_formatted_predictor_matrices',include.dirs = TRUE, full.names = TRUE)
# 
# # (a) SMOTE on detection cases
# 
# # Load detection labels and partitions
# load('../validation_resampling/detection_labels.RData')
# load('../validation_resampling/detection_partitions.RData')
# 
# for (i in 1:length(impDirs)){
#   print(paste("Imputation no.",i,"out of",length(impDirs),"started."))
#   currImpDir <- impDirs[i]
#   detection_folders <- list.files(path = currImpDir,pattern = 'detection_*',include.dirs = TRUE, full.names = TRUE)
#   for (j in 1:length(detection_folders)){
#     
#     gc()
#     
#     print(paste("Detection folder no.",j,"out of",length(detection_folders),"started."))
#     
#     curr_window_size <- as.numeric(sub(".*window_", "", detection_folders[j]))
#     curr_window_idx <- which(det_parameters$obs_windows == curr_window_size)
#     
#     curr_label_set <- det_gcs_labels[[curr_window_idx]]
#     
#     tryCatch({
#       curr_motor_train_matrix_lol <- as.data.frame(readRDS(file.path(detection_folders[j],'motor_train_matrix_lol.rds')))
#       curr_motor_train_matrix_lol$labels <- as.factor(curr_label_set$Best.Motor.Response[det_motor_train_idx[[curr_window_idx]]])
#       curr_motor_train_smote <- SmoteClassif(labels ~., dat = curr_motor_train_matrix_lol)
#       saveRDS(curr_motor_train_smote,file.path(detection_folders[j],'motor_train_smote.rds'))
#     }, error=function(e){cat("ERROR ON MOTOR, WINDOW SIZE", curr_window_size,"HOURS:",conditionMessage(e), "\n")})
#     
#     tryCatch({
#       curr_eye_train_matrix_lol <- as.data.frame(readRDS(file.path(detection_folders[j],'eye_train_matrix_lol.rds')))
#       curr_eye_train_matrix_lol$labels <- as.factor(curr_label_set$Eye.Opening[det_eye_train_idx[[curr_window_idx]]])
#       curr_eye_train_smote <- SmoteClassif(labels ~., dat = curr_eye_train_matrix_lol)
#       saveRDS(curr_eye_train_smote,file.path(detection_folders[j],'eye_train_smote.rds'))
#     }, error=function(e){cat("ERROR ON EYE, WINDOW SIZE", curr_window_size,"HOURS:",conditionMessage(e), "\n")})
#     
#     print(paste("Detection folder no.",j,"out of",length(detection_folders),"completed."))
#   }
#   print(paste("Imputation no.",i,"out of",length(impDirs),"completed."))
# }
# 
# # (b) SMOTE on prediction cases
# 
# # Load prediction labels and partitions
# load('../all_motion_feature_data/gcs_labels/prediction_labels.RData')
# load('../all_motion_feature_data/gcs_labels/prediction_partitions.RData')
# 
# for (i in 1:length(impDirs)){
#   print(paste("Imputation no.",i,"out of",length(impDirs),"started."))
#   currImpDir <- impDirs[i]
#   prediction_folders <- list.files(path = currImpDir,pattern = 'prediction_*',include.dirs = TRUE, full.names = TRUE)
#   for (j in 1:length(prediction_folders)){
#     
#     gc()
#     
#     print(paste("prediction folder no.",j,"out of",length(prediction_folders),"started."))
#  
#     pattern <- "window_\\s*(.*?)\\s*_lead"
#     curr_window_size <- as.numeric(regmatches(prediction_folders[j], regexec(pattern, prediction_folders[j]))[[1]][2])
#     curr_lead_time <- as.numeric(sub(".*lead_", "", prediction_folders[j]))
#     curr_window_idx <- which(pre_parameters$obs_windows == curr_window_size & pre_parameters$lead_times == curr_lead_time)
#     
#     curr_label_set <- pre_gcs_labels[[curr_window_idx]]
#     
#     tryCatch({
#       curr_motor_train_matrix_lol <- as.data.frame(readRDS(file.path(prediction_folders[j],'motor_train_matrix_lol.rds')))
#       curr_motor_train_matrix_lol$GCS.m <- curr_label_set$Best.Motor.Response[pre_motor_train_idx[[curr_window_idx]]]
#       curr_motor_train_matrix_lol$labels <- as.factor(curr_label_set$Net.GCSm.Change[pre_motor_train_idx[[curr_window_idx]]])
#       curr_motor_train_smote <- SmoteClassif(labels ~., dat = curr_motor_train_matrix_lol)
#       saveRDS(curr_motor_train_smote,file.path(prediction_folders[j],'motor_train_smote.rds'))
#     }, error=function(e){cat("ERROR ON MOTOR, WINDOW SIZE", curr_window_size,"HOURS:",conditionMessage(e), "\n")})
#     
#     tryCatch({
#       curr_eye_train_matrix_lol <- as.data.frame(readRDS(file.path(prediction_folders[j],'eye_train_matrix_lol.rds')))
#       curr_eye_train_matrix_lol$GCS.e <- curr_label_set$Eye.Opening[pre_eye_train_idx[[curr_window_idx]]]
#       curr_eye_train_matrix_lol$labels <- as.factor(curr_label_set$Net.GCSe.Change[pre_eye_train_idx[[curr_window_idx]]])
#       curr_eye_train_smote <- SmoteClassif(labels ~., dat = curr_eye_train_matrix_lol)
#       saveRDS(curr_eye_train_smote,file.path(prediction_folders[j],'eye_train_smote.rds'))
#     }, error=function(e){cat("ERROR ON EYE, WINDOW SIZE", curr_window_size,"HOURS:",conditionMessage(e), "\n")})
#     
#     print(paste("Prediction folder no.",j,"out of",length(prediction_folders),"completed."))
#   }
#   print(paste("Imputation no.",i,"out of",length(impDirs),"completed."))
# }

#### Analyzing LOL coefficients for feature importance
### Train and perform linear optimal low-rank projection:
rm(list = ls())
gc()

impDirs <- list.files('../all_motion_feature_data/04_formatted_predictor_matrices',pattern = 'imp*',include.dirs = TRUE, full.names = TRUE)

# (a) LOL - Detection Cases

# Load detection labels and partitions
load('../validation_resampling/detection_labels.RData')
load('../validation_resampling/detection_partitions.RData')

variables.detection.0.5.hr <- readRDS('../all_motion_feature_data/04_formatted_predictor_matrices/imp1/detection_window_0.5/motor_test_matrix.rds') %>% names()
variables.detection.1.hr <- readRDS('../all_motion_feature_data/04_formatted_predictor_matrices/imp1/detection_window_1/motor_test_matrix.rds') %>% names()
variables.detection.3.hr <- readRDS('../all_motion_feature_data/04_formatted_predictor_matrices/imp1/detection_window_3/motor_test_matrix.rds') %>% names()
variables.detection.6.hr <- readRDS('../all_motion_feature_data/04_formatted_predictor_matrices/imp1/detection_window_6/motor_test_matrix.rds') %>% names()

variables.list <- list(variables.detection.0.5.hr,variables.detection.1.hr,variables.detection.3.hr,variables.detection.6.hr)

feature.from.names <- c("band_power","freq_entropy","freq_pairs1","freq_pairs2","med_freq","sma","wavelets")
feature.to.names <- c("BPW","FDE","HLF (h)","HLF (l)","MFR","SMA","WVL")

sensor.names <- c("RE","LE","RW","LW","RA","LA")

lol.coeff.df <- as.data.frame(matrix(ncol = 9,nrow = 0))

for (i in 1:length(impDirs)){
  print(paste("Imputation no.",i,"out of",length(impDirs),"started."))
  currImpDir <- impDirs[i]
  detection_folders <- list.files(path = currImpDir,pattern = 'detection_*',include.dirs = TRUE, full.names = TRUE)
  for (j in 1:length(detection_folders)){
    print(paste("Detection folder no.",j,"out of",length(detection_folders),"started."))
    
    curr_window_size <- as.numeric(sub(".*window_", "", detection_folders[j]))
    curr.motor.lol <- readRDS(file.path(detection_folders[j],'motor_lol2.rds'))
    
    temp.df <- data.frame(full.var.name = variables.list[[j]],
                          first.lol.coeff = curr.motor.lol$A[,1],
                          second.lol.coeff = curr.motor.lol$A[,2],
                          third.lol.coeff = curr.motor.lol$A[,3],
                          imp = i,
                          observation.window = curr_window_size) %>%
      extract(full.var.name, c("temp.var", "time.count"), "(.*)_(.*)",remove = FALSE) %>%
      dplyr::select(-temp.var) %>%
      mutate(time.count = as.integer(time.count))
    
    # Assign feature types based on variable names
    temp.df$feature.type <- NA
    for (k in 1:length(feature.from.names)){
      temp.df$feature.type[grep(feature.from.names[k], temp.df$full.var.name)] <- feature.to.names[k]
    }
    
    # Assign sensor names based on variable names
    temp.df$sensor <- NA
    for (k in 1:length(sensor.names)){
      temp.df$sensor[grep(sensor.names[k], temp.df$full.var.name)] <- sensor.names[k]
    }
    
    lol.coeff.df <- rbind(lol.coeff.df,temp.df)
  }
}
    
# Save LOL compiled coefficient dataframe
saveRDS(lol.coeff.df, '../all_motion_feature_data/04_formatted_predictor_matrices/compiled_lol_coefficients.rds')

# Load dataframe of compiled LOL coefficients
lol.coeff.df <- readRDS('../all_motion_feature_data/04_formatted_predictor_matrices/compiled_lol_coefficients.rds')
lol.coeff.df$coeff <- (abs(lol.coeff.df$first.lol.coeff) + abs(lol.coeff.df$second.lol.coeff) + abs(lol.coeff.df$third.lol.coeff))/3
lol.coeff.df$extrem <- str_sub(lol.coeff.df$sensor,2,2)
lol.coeff.df$side <- str_sub(lol.coeff.df$sensor,1,1)
lol.coeff.df$upper.lower <- plyr::mapvalues(lol.coeff.df$sensor,
                                            from = c("RE","LE","RW","LW","RA","LA"),
                                            to = c('U','U','U','U','L','L'))
lol.coeff.df <- lol.coeff.df[,-c(3,4,5)]

# Contralateral test
contralat.test.df <-
  lol.coeff.df %>% 
  pivot_wider(id_cols = c("feature.type", "observation.window", "imp", "time.count","extrem"),
              names_from = "side",
              values_from = "coeff")

p.values.contralateral <- c()
for (imp.no in unique(contralat.test.df$imp)){
  curr.imp.temp.df <- contralat.test.df %>% filter(imp == imp.no)
  curr.wx.test <- wilcox.test(curr.imp.temp.df$L, curr.imp.temp.df$R, paired = TRUE,p.adjust.method = "BH",alternative = "less")
  p.values.contralateral <- c(p.values.contralateral, curr.wx.test$p.value)
}

# Upper Extremities vs. Lower Extremities Test
features.test.df <-
  lol.coeff.df %>%
  group_by(feature.type, observation.window, imp, time.count,side,upper.lower) %>%
  summarise(meanCoeff = mean(coeff)) %>%
  pivot_wider(id_cols = c("feature.type", "observation.window", "imp", "time.count","side"),
              names_from = "upper.lower",
              values_from = "meanCoeff")

p.values.up.lo <- c()
for (imp.no in unique(upper.downer.test.df$imp)){
  curr.imp.temp.df <- upper.downer.test.df %>% filter(imp == imp.no)
  curr.wx.test <- wilcox.test(curr.imp.temp.df$L, curr.imp.temp.df$U, paired = TRUE,p.adjust.method = "BH",alternative = "less")
  p.values.up.lo <- c(p.values.up.lo, curr.wx.test$p.value)
}

## Test of importance across the features
# Paired pairwise wilcoxon tests
greater.wx.test <- pairwise.wilcox.test(lol.coeff.df$coeff, lol.coeff.df$feature.type,p.adjust.method = "BH",paired = TRUE,alternative = 'greater')
lesser.wx.test <- pairwise.wilcox.test(lol.coeff.df$coeff, lol.coeff.df$feature.type,p.adjust.method = "BH",paired = TRUE,alternative = 'less')
two.sided.wx.test <- pairwise.wilcox.test(lol.coeff.df$coeff, lol.coeff.df$feature.type,p.adjust.method = "BH",paired = TRUE)

temp.df.1 <- as.data.frame(greater.wx.test$p.value)
temp.df.1$first.class <- rownames(temp.df.1)
temp.df.1 <- temp.df.1 %>% pivot_longer(cols = -first.class, names_to = 'second.class',values_to = 'p.value') %>% 
  drop_na(p.value) %>% mutate(test = 'greater')

temp.df.2 <- as.data.frame(lesser.wx.test$p.value)
temp.df.2$first.class <- rownames(temp.df.2)
temp.df.2 <- temp.df.2 %>% pivot_longer(cols = -first.class, names_to = 'second.class',values_to = 'p.value') %>% 
  drop_na(p.value) %>% mutate(test = 'less')

temp.df.3 <- as.data.frame(two.sided.wx.test$p.value)
temp.df.3$first.class <- rownames(temp.df.3)
temp.df.3 <- temp.df.3 %>% pivot_longer(cols = -first.class, names_to = 'second.class',values_to = 'p.value') %>% 
  drop_na(p.value) %>% mutate(test = 'two.tailed')

feature.compiled.wx.p.values <- rbind(temp.df.1,temp.df.2,temp.df.3)
write.csv(feature.compiled.wx.p.values,'../results/lol_coefficient_stats/feature_wx_tests.csv',row.names = FALSE)

# Kruskall-Wallis tests across different, relevant subsets
compiled.kw.p.values <- as.data.frame(matrix(ncol = 4, nrow = 0))
for (i in 1:length(unique(lol.coeff.df$imp))){
  curr.imp <- unique(lol.coeff.df$imp)[i]
  for (j in 1:length(unique(lol.coeff.df$observation.window))){
    curr.o.w <- unique(lol.coeff.df$observation.window)[j]
    for (k in 1:length(unique(lol.coeff.df$sensor))){
      curr.sensor <- unique(lol.coeff.df$sensor)[k]
      curr.coeff.set <- lol.coeff.df %>% filter(imp == curr.imp &
                                                  observation.window ==  curr.o.w &
                                                  sensor == curr.sensor)
      curr.kw.test <- kruskal.test(coeff ~ feature.type, data = curr.coeff.set)
      compiled.kw.p.values <- rbind(compiled.kw.p.values,
                                    data.frame(imp = curr.imp,
                                               observation.window = curr.o.w,
                                               sensor = curr.sensor,
                                               p.value = curr.kw.test$p.value))
    }
  }
}
write.csv(compiled.kw.p.values,'../results/lol_coefficient_stats/feature_kw_tests.csv',row.names = FALSE)

## Test of importance across the sensors
# Paired pairwise wilcoxon tests
sensor.greater.wx.test <- pairwise.wilcox.test(lol.coeff.df$coeff, lol.coeff.df$sensor,p.adjust.method = "BH",paired = TRUE,alternative = 'greater')
sensor.lesser.wx.test <- pairwise.wilcox.test(lol.coeff.df$coeff, lol.coeff.df$sensor,p.adjust.method = "BH",paired = TRUE,alternative = 'less')
sensor.two.sided.wx.test <- pairwise.wilcox.test(lol.coeff.df$coeff, lol.coeff.df$sensor.type,p.adjust.method = "BH",paired = TRUE)

temp.df.1 <- as.data.frame(sensor.greater.wx.test$p.value)
temp.df.1$first.class <- rownames(temp.df.1)
temp.df.1 <- temp.df.1 %>% pivot_longer(cols = -first.class, names_to = 'second.class',values_to = 'p.value') %>% 
  drop_na(p.value) %>% mutate(test = 'greater')

temp.df.2 <- as.data.frame(sensor.lesser.wx.test$p.value)
temp.df.2$first.class <- rownames(temp.df.2)
temp.df.2 <- temp.df.2 %>% pivot_longer(cols = -first.class, names_to = 'second.class',values_to = 'p.value') %>% 
  drop_na(p.value) %>% mutate(test = 'less')

temp.df.3 <- as.data.frame(sensor.lesser.wx.test$p.value)
temp.df.3$first.class <- rownames(temp.df.3)
temp.df.3 <- temp.df.3 %>% pivot_longer(cols = -first.class, names_to = 'second.class',values_to = 'p.value') %>% 
  drop_na(p.value) %>% mutate(test = 'two.tailed')

sensor.compiled.wx.p.values <- rbind(temp.df.1,temp.df.2,temp.df.3)
write.csv(sensor.compiled.wx.p.values,'../results/lol_coefficient_stats/sensor_wx_tests.csv',row.names = FALSE)

# Kruskall-Wallis tests across different, relevant subsets
sensors.compiled.kw.p.values <- as.data.frame(matrix(ncol = 4, nrow = 0))
for (i in 1:length(unique(lol.coeff.df$imp))){
  curr.imp <- unique(lol.coeff.df$imp)[i]
  for (j in 1:length(unique(lol.coeff.df$observation.window))){
    curr.o.w <- unique(lol.coeff.df$observation.window)[j]
    for (k in 1:length(unique(lol.coeff.df$feature.type))){
      curr.feature.type <- unique(lol.coeff.df$feature.type)[k]
      curr.coeff.set <- lol.coeff.df %>% filter(imp == curr.imp &
                                                  observation.window ==  curr.o.w &
                                                  feature.type == curr.feature.type)
      curr.kw.test <- kruskal.test(coeff ~ sensor, data = curr.coeff.set)
      sensors.compiled.kw.p.values <- rbind(sensors.compiled.kw.p.values,
                                    data.frame(imp = curr.imp,
                                               observation.window = curr.o.w,
                                               feature.type = curr.feature.type,
                                               p.value = curr.kw.test$p.value))
    }
  }
}
write.csv(sensors.compiled.kw.p.values,'../results/lol_coefficient_stats/sensor_kw_tests.csv',row.names = FALSE)

## Test of importance across the time bins before GCSm evaluation
time.before.compiled.kw.p.values <- as.data.frame(matrix(ncol = 5, nrow = 0))
time.before.compiled.wx.p.values <- as.data.frame(matrix(ncol = 8, nrow = 0))

for (i in 1:length(unique(lol.coeff.df$imp))){
  print(paste('imputation',i,'started'))
  curr.imp <- unique(lol.coeff.df$imp)[i]
  for (j in 1:length(unique(lol.coeff.df$observation.window))){
    print(paste('observation.window',j,'started'))
    curr.o.w <- unique(lol.coeff.df$observation.window)[j]
    if (curr.o.w == 0.5){
      next
    }
    for (k in 1:length(unique(lol.coeff.df$feature.type))){
      print(paste('feature.type',k,'started'))
      curr.feature.type <- unique(lol.coeff.df$feature.type)[k]
      for (m in  1:length(unique(lol.coeff.df$sensor))){
        print(paste('sensor',m,'started'))
        curr.sensor <- unique(lol.coeff.df$sensor)[m]
        curr.coeff.set <- lol.coeff.df %>% filter(imp == curr.imp &
                                                    observation.window ==  curr.o.w &
                                                    feature.type == curr.feature.type &
                                                    sensor == curr.sensor)
        curr.coeff.set$time.bin = NA
        uniq.time.counts <- sort(unique(curr.coeff.set$time.count)) %>% .[-length(.)]
        for (timeBin in 1:(curr.o.w/.5)){
          curr.time.idx <- ((curr.o.w/.5 - timeBin)*360) + 1:360
          curr.coeff.set$time.bin[curr.coeff.set$time.count %in% curr.time.idx] = timeBin
        }
        
        curr.coeff.set <- curr.coeff.set %>% drop_na(time.bin)
        
        curr.kw.test <- kruskal.test(coeff ~ time.bin, data = curr.coeff.set)
        time.before.compiled.kw.p.values <- rbind(time.before.compiled.kw.p.values,
                                                  data.frame(imp = curr.imp,
                                                             observation.window = curr.o.w,
                                                             feature.type = curr.feature.type,
                                                             sensor = curr.sensor,
                                                             p.value = curr.kw.test$p.value))
        
        curr.time.bin.wx.greater.test <- pairwise.wilcox.test(curr.coeff.set$coeff, curr.coeff.set$time.bin,p.adjust.method = "BH",alternative = 'greater')
        curr.time.bin.wx.less.test <- pairwise.wilcox.test(curr.coeff.set$coeff, curr.coeff.set$time.bin,p.adjust.method = "BH",alternative = 'less')
        curr.time.bin.wx.two.tailed.test <- pairwise.wilcox.test(curr.coeff.set$coeff, curr.coeff.set$time.bin,p.adjust.method = "BH")
        
        temp.df.1 <- as.data.frame(curr.time.bin.wx.greater.test$p.value)
        temp.df.1$first.class <- rownames(temp.df.1)
        temp.df.1 <- temp.df.1 %>% pivot_longer(cols = -first.class, names_to = 'second.class',values_to = 'p.value') %>% 
          drop_na(p.value) %>% mutate(test = 'greater',
                                      imp = curr.imp,
                                      observation.window = curr.o.w,
                                      feature.type = curr.feature.type,
                                      sensor = curr.sensor)

        temp.df.2 <- as.data.frame(curr.time.bin.wx.less.test$p.value)
        temp.df.2$first.class <- rownames(temp.df.2)
        temp.df.2 <- temp.df.2 %>% pivot_longer(cols = -first.class, names_to = 'second.class',values_to = 'p.value') %>% 
          drop_na(p.value) %>% mutate(test = 'less',
                                      imp = curr.imp,
                                      observation.window = curr.o.w,
                                      feature.type = curr.feature.type,
                                      sensor = curr.sensor)

        temp.df.3 <- as.data.frame(curr.time.bin.wx.two.tailed.test$p.value)
        temp.df.3$first.class <- rownames(temp.df.3)
        temp.df.3 <- temp.df.3 %>% pivot_longer(cols = -first.class, names_to = 'second.class',values_to = 'p.value') %>% 
          drop_na(p.value) %>% mutate(test = 'two.tailed',
                                      imp = curr.imp,
                                      observation.window = curr.o.w,
                                      feature.type = curr.feature.type,
                                      sensor = curr.sensor)
        
        time.before.compiled.wx.p.values <- rbind(time.before.compiled.wx.p.values,
                                                  temp.df.1,
                                                  temp.df.2,
                                                  temp.df.3)
      }
    }
  }
}
write.csv(time.before.compiled.kw.p.values,'../results/lol_coefficient_stats/time_kw_tests.csv',row.names = FALSE)
write.csv(time.before.compiled.wx.p.values,'../results/lol_coefficient_stats/time_wx_tests.csv',row.names = FALSE)

# Pool p-values across imputations with z-transform
grouped.t.b.compiled.kw.p.values <- time.before.compiled.kw.p.values %>% 
  mutate(z.norm = qnorm(p.value, lower.tail=FALSE)) %>%
  group_by(observation.window, feature.type, sensor) %>% 
  summarise(mean.p.value = mean(p.value),
            std.p.value = sd(p.value),
            mean.z.norm = mean(z.norm),
            var.z.norm = 1 + (1 + (1/n()))*var(z.norm),
            r.m = (1 + (1/n()))*var(z.norm),
            m = n()) %>%
  mutate(d.o.f. = (m-1)*((1 + (1/r.m))^2)) %>%
  mutate(p.m = pt(q = mean.z.norm, df = d.o.f., lower.tail = FALSE)) %>%
  mutate(p.value.imp = sprintf("%.2f",p.m))

time.before.compiled.wx.p.values$p.value[time.before.compiled.wx.p.values$p.value == 1] <- .999
grouped.t.b.compiled.wx.p.values <- time.before.compiled.wx.p.values %>% 
  mutate(z.norm = qnorm(p.value, lower.tail=FALSE)) %>%
  group_by(observation.window, first.class,second.class,test,feature.type, sensor) %>% 
  summarise(mean.p.value = mean(p.value),
            std.p.value = sd(p.value),
            mean.z.norm = mean(z.norm),
            var.z.norm = 1 + (1 + (1/n()))*var(z.norm),
            r.m = (1 + (1/n()))*var(z.norm),
            m = n()) %>%
  mutate(d.o.f. = (m-1)*((1 + (1/r.m))^2)) %>%
  mutate(p.m = pt(q = mean.z.norm, df = d.o.f., lower.tail = FALSE)) %>%
  mutate(p.value.imp = sprintf("%.2f",p.m))

# filter of wilcoxon values
filt.t.b.wx.p.values <- grouped.t.b.compiled.wx.p.values %>%
  filter(test != 'two.tailed') %>%
  filter(p.m <= 0.05)

less.filt.t.b.wx.p.values <- filt.t.b.wx.p.values %>%
  filter(test == 'less') %>%
  group_by(observation.window, first.class,second.class) %>%
  summarise(n = n())

greater.filt.t.b.wx.p.values <- filt.t.b.wx.p.values %>%
  filter(test == 'greater') %>%
  group_by(observation.window, first.class,second.class) %>%
  summarise(n = n())