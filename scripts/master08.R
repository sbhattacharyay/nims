#### Master Script 8: Embedding and Over-Sampling ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

# Add library path for R packages and set working directory:
.libPaths(c("~/Rpackages" , .libPaths()))
setwd('~/work/nims/scripts')

# Load necessary packages
library(tidyverse)
library(readxl)
library(caret)
library(lolR)

# Load clinical patient data:
source('./functions/load_patient_clinical_data.R')
patient_clinical_data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv') %>% arrange(AccelPatientNo_) %>% mutate(ptIdx = 1:nrow(.))

#  Format feature spaces into matrices and save matrices

dir.create("~/scratch/all_motion_feature_data/formatted_matrices",showWarnings = FALSE)
impFiles <- list.files(path="~/scratch/all_motion_feature_data/final_imputed_features/",pattern = "*.RData")

# (a) Detection model matrices

# Load detection labels and partitions
load('~/scratch/all_motion_feature_data/gcs_labels/detection_labels.RData')
load('~/scratch/all_motion_feature_data/gcs_labels/detection_partitions.RData')

for (impNo in 1:length(impFiles)){
  print(paste("Imputation no.",impNo,"out of",length(impFiles),"started."))
  dir.create(paste0("~/scratch/all_motion_feature_data/formatted_matrices/imp",impNo),showWarnings = FALSE)
  load(paste0('~/scratch/all_motion_feature_data/final_imputed_features/',impFiles[impNo]))
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
    curr_det_motor_test_idx <- setdiff(motor_indices,curr_det_motor_train_idx)
    curr_det_eye_test_idx <- setdiff(eye_indices,curr_det_eye_train_idx)
    dir.create(paste0('~/scratch/all_motion_feature_data/formatted_matrices/imp',impNo,'/detection_window_',curr_obs_window),showWarnings = FALSE)
    
    # Motor train matrix:
    curr_motor_train_matrix <- matrix(nrow = length(curr_det_motor_train_idx), ncol = ((curr_secs_window/5)+1)*7*6)
    print("Motor training matrix started:")
    for (j in 1:length(curr_det_motor_train_idx)){
      curr_GCS_row <- curr_det_motor_train_idx[j]
      currTimeStamp <- curr_det_gcs_labels$TakenInstant[curr_GCS_row]
      currFeatSet <- currImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == curr_det_gcs_labels$AccelPatientNo_[curr_GCS_row], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - curr_secs_window)
      tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
      tempReshape <- as.vector(as.matrix(tempReshape[,-1]))
      curr_motor_train_matrix[j,] <- tempReshape
      print(paste("Row no.",j,"out of",length(curr_det_motor_train_idx),"complete."))
    }
    saveRDS(curr_motor_train_matrix,file = paste0('~/scratch/all_motion_feature_data/formatted_matrices/imp',impNo,'/detection_window_',curr_obs_window,'/motor_train_matrix.rds'))
    
    # Eye train matrix:
    curr_eye_train_matrix <- matrix(nrow = length(curr_det_eye_train_idx), ncol = ((curr_secs_window/5)+1)*7*6)
    print("Eye training matrix started:")
    for (j in 1:length(curr_det_eye_train_idx)){
      curr_GCS_row <- curr_det_eye_train_idx[j]
      currTimeStamp <- curr_det_gcs_labels$TakenInstant[curr_GCS_row]
      currFeatSet <- currImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == curr_det_gcs_labels$AccelPatientNo_[curr_GCS_row], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - curr_secs_window)
      tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
      tempReshape <- as.vector(as.matrix(tempReshape[,-1]))
      curr_eye_train_matrix[j,] <- tempReshape
      print(paste("Row no.",j,"out of",length(curr_det_eye_train_idx),"complete."))
    }
    saveRDS(curr_eye_train_matrix,file = paste0('~/scratch/all_motion_feature_data/formatted_matrices/imp',impNo,'/detection_window_',curr_obs_window,'/eye_train_matrix.rds'))
    
    # Motor test matrix:
    curr_motor_test_matrix <- matrix(nrow = length(curr_det_motor_test_idx), ncol = ((curr_secs_window/5)+1)*7*6)
    print("Motor testing matrix started:")
    for (j in 1:length(curr_det_motor_test_idx)){
      curr_GCS_row <- curr_det_motor_test_idx[j]
      currTimeStamp <- curr_det_gcs_labels$TakenInstant[curr_GCS_row]
      currFeatSet <- currImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == curr_det_gcs_labels$AccelPatientNo_[curr_GCS_row], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - curr_secs_window)
      tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
      tempReshape <- as.vector(as.matrix(tempReshape[,-1]))
      curr_motor_test_matrix[j,] <- tempReshape
      print(paste("Row no.",j,"out of",length(curr_det_motor_test_idx),"complete."))
    }
    saveRDS(curr_motor_test_matrix,file = paste0('~/scratch/all_motion_feature_data/formatted_matrices/imp',impNo,'/detection_window_',curr_obs_window,'/motor_test_matrix.rds'))
    
    # Eye test matrix:
    curr_eye_test_matrix <- matrix(nrow = length(curr_det_eye_test_idx), ncol = ((curr_secs_window/5)+1)*7*6)
    print("Eye testing matrix started:")
    for (j in 1:length(curr_det_eye_test_idx)){
      curr_GCS_row <- curr_det_eye_test_idx[j]
      currTimeStamp <- curr_det_gcs_labels$TakenInstant[curr_GCS_row]
      currFeatSet <- currImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == curr_det_gcs_labels$AccelPatientNo_[curr_GCS_row], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - curr_secs_window)
      tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
      tempReshape <- as.vector(as.matrix(tempReshape[,-1]))
      curr_eye_test_matrix[j,] <- tempReshape
      print(paste("Row no.",j,"out of",length(curr_det_eye_test_idx),"complete."))
    }
    saveRDS(curr_eye_test_matrix,file = paste0('~/scratch/all_motion_feature_data/formatted_matrices/imp',impNo,'/detection_window_',curr_obs_window,'/eye_test_matrix.rds'))
    
    print(paste("Parameter combination no.",i,"out of",length(det_gcs_labels),"started."))
  }
  print(paste("Imputation no.",impNo,"out of",length(impFiles),"completed."))
}

# (b) Prediction model matrices

# Load prediction labels and partitions
load('~/scratch/all_motion_feature_data/gcs_labels/prediction_labels.RData')
load('~/scratch/all_motion_feature_data/gcs_labels/prediction_partitions.RData')

for (impNo in 1:length(impFiles)){
  print(paste("Imputation no.",impNo,"out of",length(impFiles),"started."))
  dir.create(paste0("~/scratch/all_motion_feature_data/formatted_matrices/imp",impNo),showWarnings = FALSE)
  load(paste0('~/scratch/all_motion_feature_data/final_imputed_features/',impFiles[impNo]))
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
    curr_pre_motor_test_idx <- setdiff(motor_indices,curr_pre_motor_train_idx)
    curr_pre_eye_test_idx <- setdiff(eye_indices,curr_pre_eye_train_idx)
    dir.create(paste0('~/scratch/all_motion_feature_data/formatted_matrices/imp',impNo,'/prediction_window_',curr_obs_window),showWarnings = FALSE)
    
    # Motor train matrix:
    curr_motor_train_matrix <- matrix(nrow = length(curr_pre_motor_train_idx), ncol = ((curr_secs_window/5)+1)*7*6)
    print("Motor training matrix started:")
    for (j in 1:length(curr_pre_motor_train_idx)){
      curr_GCS_row <- curr_pre_motor_train_idx[j]
      currTimeStamp <- curr_pre_gcs_labels$TakenInstant[curr_GCS_row]
      currFeatSet <- currImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == curr_pre_gcs_labels$AccelPatientNo_[curr_GCS_row], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - curr_secs_window)
      tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
      tempReshape <- as.vector(as.matrix(tempReshape[,-1]))
      curr_motor_train_matrix[j,] <- tempReshape
      print(paste("Row no.",j,"out of",length(curr_pre_motor_train_idx),"complete."))
    }
    saveRDS(curr_motor_train_matrix,file = paste0('~/scratch/all_motion_feature_data/formatted_matrices/imp',impNo,'/prediction_window_',curr_obs_window,'/motor_train_matrix.rds'))
    
    # Eye train matrix:
    curr_eye_train_matrix <- matrix(nrow = length(curr_pre_eye_train_idx), ncol = ((curr_secs_window/5)+1)*7*6)
    print("Eye training matrix started:")
    for (j in 1:length(curr_pre_eye_train_idx)){
      curr_GCS_row <- curr_pre_eye_train_idx[j]
      currTimeStamp <- curr_pre_gcs_labels$TakenInstant[curr_GCS_row]
      currFeatSet <- currImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == curr_pre_gcs_labels$AccelPatientNo_[curr_GCS_row], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - curr_secs_window)
      tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
      tempReshape <- as.vector(as.matrix(tempReshape[,-1]))
      curr_eye_train_matrix[j,] <- tempReshape
      print(paste("Row no.",j,"out of",length(curr_pre_eye_train_idx),"complete."))
    }
    saveRDS(curr_eye_train_matrix,file = paste0('~/scratch/all_motion_feature_data/formatted_matrices/imp',impNo,'/prediction_window_',curr_obs_window,'/eye_train_matrix.rds'))
    
    # Motor test matrix:
    curr_motor_test_matrix <- matrix(nrow = length(curr_pre_motor_test_idx), ncol = ((curr_secs_window/5)+1)*7*6)
    print("Motor testing matrix started:")
    for (j in 1:length(curr_pre_motor_test_idx)){
      curr_GCS_row <- curr_pre_motor_test_idx[j]
      currTimeStamp <- curr_pre_gcs_labels$TakenInstant[curr_GCS_row]
      currFeatSet <- currImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == curr_pre_gcs_labels$AccelPatientNo_[curr_GCS_row], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - curr_secs_window)
      tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
      tempReshape <- as.vector(as.matrix(tempReshape[,-1]))
      curr_motor_test_matrix[j,] <- tempReshape
      print(paste("Row no.",j,"out of",length(curr_pre_motor_test_idx),"complete."))
    }
    saveRDS(curr_motor_test_matrix,file = paste0('~/scratch/all_motion_feature_data/formatted_matrices/imp',impNo,'/prediction_window_',curr_obs_window,'/motor_test_matrix.rds'))
    
    # Eye test matrix:
    curr_eye_test_matrix <- matrix(nrow = length(curr_pre_eye_test_idx), ncol = ((curr_secs_window/5)+1)*7*6)
    print("Eye testing matrix started:")
    for (j in 1:length(curr_pre_eye_test_idx)){
      curr_GCS_row <- curr_pre_eye_test_idx[j]
      currTimeStamp <- curr_pre_gcs_labels$TakenInstant[curr_GCS_row]
      currFeatSet <- currImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == curr_pre_gcs_labels$AccelPatientNo_[curr_GCS_row], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - curr_secs_window)
      tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
      tempReshape <- as.vector(as.matrix(tempReshape[,-1]))
      curr_eye_test_matrix[j,] <- tempReshape
      print(paste("Row no.",j,"out of",length(curr_pre_eye_test_idx),"complete."))
    }
    saveRDS(curr_eye_test_matrix,file = paste0('~/scratch/all_motion_feature_data/formatted_matrices/imp',impNo,'/prediction_window_',curr_obs_window,'/eye_test_matrix.rds'))
    
    print(paste("Parameter combination no.",i,"out of",length(pre_gcs_labels),"started."))
  }
  print(paste("Imputation no.",impNo,"out of",length(impFiles),"completed."))
}
