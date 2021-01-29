#### Resampling of GCS data for classification ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu
#
# Load necessary packages
.libPaths(c("~/Rpackages" , .libPaths()))
setwd("~/work/nims/scripts")

library(tidyverse)
library(readxl)
library(caret)
source('./functions/load_patient_clinical_data.R')

# Load clinical patient data:
source('./functions/load_patient_clinical_data.R')
patient.clinical.data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv') %>% arrange(AccelPatientNo_) %>% mutate(ptIdx = 1:nrow(.))

# Load missing accelerometry information:
missing.data.info <- read_xlsx('~/scratch/all_motion_feature_data/MissingPercentTable.xlsx',.name_repair = "universal") %>% mutate(Start.Timestamp = as.POSIXct(Start.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"), End.Timestamp = as.POSIXct(End.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"))

# Load automatically extracted GCS labels:
gcs_data <- read.csv('../clinical_data/02_clean_auto_GCS_table.csv') %>% select(-X) %>% mutate(TakenInstant = as.POSIXct(TakenInstant, tz = "America/New_York"))

# Concurrent GCS detection model parameters:
det_obs_windows <- c(.5, 1, 3, 6) #in hours
det_secs_windows <- det_obs_windows*3600
det_no_pts <- det_obs_windows*60*12 + 1
default_det_idx <- 2

# Future GCS prediction model parameters:
pre_obs_windows <- c(.5, 1, 3, 6) #in hours
pre_lead_times <- c(0,1,2,6) #in hours
pre_pre_window <- 24 #in hours
pre_secs_windows <- pre_obs_windows*3600
pre_secs_leads <- pre_lead_times*3600
pre_secs_pred <- pre_pre_window*3600
pre_no_pts <- pre_obs_windows*60*12 + 1
default_pre_window_idx <- 2
default_pre_leads_idx <- 1

# Isolate and save candidate (coincide with accelerometry) GCS observations:

dir.create('../validation_resampling',showWarnings = FALSE)

# a) Detection models:

det_parameters <- as.data.frame(matrix(nrow = length(det_obs_windows),ncol = 0))
det_gcs_labels <- vector(mode = "list")
for (i in 1:length(det_obs_windows)) {
  det_parameters$obs_windows[i] <- det_obs_windows[i]
  print(paste("Observation Window Size:",det_obs_windows[i],"hours started"))
  curr_secs_window <- det_secs_windows[i]
  gcsLabels <- as.data.frame(matrix(nrow = 0,ncol = 8))
  for (patIdx in unique(patient.clinical.data$ptIdx)) {
    print(paste("Patient no.",patIdx,"started"))
    currStudyNo <- patient.clinical.data$StudyPatientNo_[patIdx]
    currAccelNo <- patient.clinical.data$AccelPatientNo_[patIdx]
    currStartTm <- missing.data.info$Start.Timestamp[missing.data.info$Accel.Patient.No. == currAccelNo]
    currEndTm <- missing.data.info$End.Timestamp[missing.data.info$Accel.Patient.No. == currAccelNo]
    currGCS <- gcs_data %>% filter(StudyPatientNo_ == currStudyNo, AccelPatientNo_ == currAccelNo)
    gcsFilter <- currGCS$TakenInstant >= currStartTm+curr_secs_window & currGCS$TakenInstant <= currEndTm
    gcsLabels <- rbind(gcsLabels,currGCS[gcsFilter,])
  }
  det_gcs_labels[[i]] <- gcsLabels
}
save(det_parameters, det_gcs_labels, file = '../validation_resampling/detection_labels.RData')

# b) Prediction models:

pre_parameters <- as.data.frame(matrix(nrow = length(pre_obs_windows)*length(pre_lead_times),ncol = 0))
pre_gcs_labels <- vector(mode = "list")
counter <- 0

for (i in 1:length(pre_lead_times)){
  curr_secs_lead <- pre_secs_leads[i]
  for (j in 1:length(pre_obs_windows)) {
    counter <- counter + 1
    pre_parameters$lead_times[counter] <- pre_lead_times[i]
    pre_parameters$obs_windows[counter] <- pre_obs_windows[j]
    print(paste("Observation Lead Time:",pre_lead_times[i],"hours, Window Size:",pre_obs_windows[j],"hours started"))
    curr_secs_window <- pre_secs_windows[j]
    gcsLabels <- as.data.frame(matrix(nrow = 0,ncol = 14))
    for (patIdx in unique(patient.clinical.data$ptIdx)) {
      print(paste("Patient no.",patIdx,"started"))
      currStudyNo <- patient.clinical.data$StudyPatientNo_[patIdx]
      currAccelNo <- patient.clinical.data$AccelPatientNo_[patIdx]
      currStartTm <- missing.data.info$Start.Timestamp[missing.data.info$Accel.Patient.No. == currAccelNo]
      currEndTm <- missing.data.info$End.Timestamp[missing.data.info$Accel.Patient.No. == currAccelNo]
      currGCS <- gcs_data %>% filter(StudyPatientNo_ == currStudyNo, AccelPatientNo_ == currAccelNo)
      gcsFilter <- currGCS$TakenInstant >= currStartTm+curr_secs_window+curr_secs_lead & currGCS$TakenInstant <= currEndTm+curr_secs_lead
      
      if (all(!gcsFilter)){
        print(paste("Patient no.",patIdx,"Empty: completed"))
        next
      }
      
      currFiltGCS <- currGCS[gcsFilter,]
      predWindowLabels <- as.data.frame(matrix(nrow = 0,ncol = 6),stringsAsFactors = FALSE)
      for (l in 1:nrow(currFiltGCS)){
        currTimeStamp <- currFiltGCS$TakenInstant[l]
        startGCSm <- currFiltGCS$Best.Motor.Response[l]
        startGCSe <- currFiltGCS$Eye.Opening[l]
        predWindowFilter <- currGCS$TakenInstant >= currTimeStamp & currGCS$TakenInstant <= currTimeStamp+pre_secs_pred
        currPredWindowGCS <- currGCS[predWindowFilter,]
        
        if (!is.na(startGCSm)) {
          maxPosGCSm <- -min(startGCSm - currPredWindowGCS$Best.Motor.Response, na.rm = TRUE)
          maxNegGCSm <- -max(startGCSm - currPredWindowGCS$Best.Motor.Response, na.rm = TRUE)
          
          if (is.na(maxPosGCSm) & is.na(maxNegGCSm)){
            GCSm_change <- NA
          } else if (maxPosGCSm == maxNegGCSm & maxNegGCSm == 0){
            GCSm_change <- "no.change"
          } else if (abs(maxPosGCSm) == abs(maxNegGCSm)){
            GCSm_change <- "decrease" # tiebreaker case: assume decrease is more significant
          } else if (abs(maxPosGCSm) > abs(maxNegGCSm)){
            GCSm_change <- "increase"
          } else if (abs(maxPosGCSm) < abs(maxNegGCSm)){
            GCSm_change <- "decrease"
          }
        } else {
          maxPosGCSm <- NA
          maxNegGCSm <- NA
          GCSm_change <- NA
        }
        
        if (!is.na(startGCSe)) {
          maxPosGCSe <- -min(startGCSe - currPredWindowGCS$Eye.Opening, na.rm = TRUE)
          maxNegGCSe <- -max(startGCSe - currPredWindowGCS$Eye.Opening, na.rm = TRUE)
          if (is.na(maxPosGCSe) & is.na(maxNegGCSe)){
            GCSe_change <- NA
          } else if (maxPosGCSe == maxNegGCSe & maxNegGCSe == 0){
            GCSe_change <- "no.change"
          } else if (abs(maxPosGCSe) == abs(maxNegGCSe)){
            GCSe_change <- "decrease" # tiebreaker case: assume decrease is more significant
          } else if (abs(maxPosGCSe) > abs(maxNegGCSe)){
            GCSe_change <- "increase"
          } else if (abs(maxPosGCSe) < abs(maxNegGCSe)){
            GCSe_change <- "decrease"
          }
        } else {
          maxPosGCSe <- NA
          maxNegGCSe <- NA
          GCSe_change <- NA
        }
        
        predWindowLabels <- rbind(predWindowLabels,list(maxPosGCSm,maxNegGCSm,maxPosGCSe,maxNegGCSe,GCSm_change,GCSe_change))
        predWindowLabels[,5] <- as.character(predWindowLabels[,5])
        predWindowLabels[,6] <- as.character(predWindowLabels[,6])
      }
      names(predWindowLabels) <- c("Max.GCSm.Increase","Max.GCSm.Decrease","Max.GCSe.Increase","Max.GCSe.Decrease","Net.GCSm.Change","Net.GCSe.Change")
      gcsLabels <- rbind(gcsLabels,cbind(currFiltGCS,predWindowLabels))
    }
    pre_gcs_labels[[counter]] <- gcsLabels
  }
}
save(pre_parameters, pre_gcs_labels, file = '../validation_resampling/prediction_labels.RData')

# Partition data into training, validation, and testing sets based on preserving class imbalance:
rm(list = ls())
gc()

p_train <- 0.8 # set proportion of observations in training (including validation) sets

# (a) Detection GCS labels:

load('../validation_resampling/detection_labels.RData')

det_motor_train_idx <- vector(mode = "list")
det_eye_train_idx <- vector(mode = "list")

det_motor_test_idx <- vector(mode = "list")
det_eye_test_idx <- vector(mode = "list")

set.seed(2020)
for (i in 1:length(det_gcs_labels)){
  currGCS <- det_gcs_labels[[i]]
  
  motor_nonmissingIdx <- which(!is.na(currGCS$Best.Motor.Response))
  eye_nonmissingIdx <- which(!is.na(currGCS$Eye.Opening))
  
  currGCSm <- currGCS[motor_nonmissingIdx,]
  currGCSe <- currGCS[eye_nonmissingIdx,]
  
  currGCSm_trainIdx <- createDataPartition(currGCSm$Best.Motor.Response,p=p_train,list = FALSE,times = 1)
  currGCSe_trainIdx <- createDataPartition(currGCSe$Eye.Opening,p=p_train,list = FALSE,times = 1)
  
  det_motor_train_idx[[i]] <- motor_nonmissingIdx[currGCSm_trainIdx]
  det_eye_train_idx[[i]] <- eye_nonmissingIdx[currGCSe_trainIdx]
  
  det_motor_test_idx[[i]] <- motor_nonmissingIdx[!motor_nonmissingIdx %in% motor_nonmissingIdx[currGCSm_trainIdx]]
  det_eye_test_idx[[i]] <- eye_nonmissingIdx[!eye_nonmissingIdx %in% eye_nonmissingIdx[currGCSe_trainIdx]]
}
save(det_motor_train_idx, det_eye_train_idx,det_motor_test_idx,det_eye_test_idx, file = '../validation_resampling/detection_partitions.RData')

# (b) Prediction GCS labels:

load('../validation_resampling/prediction_labels.RData')

pre_motor_train_idx <- vector(mode = "list")
pre_eye_train_idx <- vector(mode = "list")

pre_motor_test_idx <- vector(mode = "list")
pre_eye_test_idx <- vector(mode = "list")

set.seed(2020)
for (i in 1:length(pre_gcs_labels)){
  currGCS <- pre_gcs_labels[[i]]
  
  motor_nonmissingIdx <- which(!is.na(currGCS$Net.GCSm.Change))
  eye_nonmissingIdx <- which(!is.na(currGCS$Net.GCSe.Change))
  
  currGCSm <- currGCS[motor_nonmissingIdx,]
  currGCSe <- currGCS[eye_nonmissingIdx,]
  
  currGCSm_trainIdx <- createDataPartition(currGCSm$Net.GCSm.Change,p=p_train,list = FALSE,times = 1)
  currGCSe_trainIdx <- createDataPartition(currGCSe$Net.GCSe.Change,p=p_train,list = FALSE,times = 1)
  
  pre_motor_train_idx[[i]] <- motor_nonmissingIdx[currGCSm_trainIdx]
  pre_eye_train_idx[[i]] <- eye_nonmissingIdx[currGCSe_trainIdx]
  
  pre_motor_test_idx[[i]] <- motor_nonmissingIdx[!motor_nonmissingIdx %in% motor_nonmissingIdx[currGCSm_trainIdx]]
  pre_eye_test_idx[[i]] <- eye_nonmissingIdx[!eye_nonmissingIdx %in% eye_nonmissingIdx[currGCSe_trainIdx]]
}
save(pre_motor_train_idx, pre_eye_train_idx,pre_motor_test_idx,pre_eye_test_idx, file = '../validation_resampling/prediction_partitions.RData')
