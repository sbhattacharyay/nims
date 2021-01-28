#### Master Script 7: Preparation and resampling of GCS data for classification ####
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

# Load automatically extracted GCS data:
gcs_table <- read_xlsx('../clinical_data/automatic_GCS_table.xlsx') %>% drop_na(Value) %>% pivot_wider(id_cols = c(TakenInstant,PrimaryMrn,BirthDate,LastName,FirstName), names_from = DisplayName, values_from = Value)
patient_clinical_data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv')
accel_no_match <- patient_clinical_data %>% select(StudyPatientNo_,AccelPatientNo_,DOB)
id_gcs_table <- left_join(gcs_table,accel_no_match,by = c("BirthDate" = "DOB")) %>% drop_na(StudyPatientNo_)
missingTable <- read_xlsx('../all_motion_feature_data/MissingPercentTable.xlsx',.name_repair = "universal") %>% mutate(Start.Timestamp = as.POSIXct(Start.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"), End.Timestamp = as.POSIXct(End.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"))

# Add indicator of whether the data point coincides with accelerometry recording:
id_gcs_table$Coincide.with.Accel.Recording <- NA
for (i in 1:nrow(missingTable)){
  currAccelPatientNo <- missingTable$Study.Patient.No.[i]
  accelPatientIdx <- which(id_gcs_table$AccelPatientNo_ == currAccelPatientNo)
  filtGCS <- id_gcs_table[accelPatientIdx,]
  id_gcs_table$Coincide.with.Accel.Recording[accelPatientIdx] = filtGCS$TakenInstant <= missingTable$End.Timestamp[i] & filtGCS$TakenInstant >= missingTable$Start.Timestamp[i]
}

finalGCSTable <- id_gcs_table[,c(10,11,1,6,7,8,9,12)]
names(finalGCSTable) <- make.names(names(finalGCSTable))
write.csv(finalGCSTable,'../clinical_data/clean_auto_GCS_table.csv')

# Load manually extracted GCS data:
man_table <- read_xlsx('../clinical_data/manual_GCS_table.xlsx') %>% mutate(Coincide.with.Accel.Recording = `Coincide with accel` == 1,Time...3 = sprintf("%04d", Time...3)) %>% filter(`Study Patient Number` != 1) %>% mutate(combTime = paste(as.character(Date),Time...3)) %>% mutate(TakenInstant = as.POSIXct(combTime,format = '%Y-%m-%d %H%M'))
man_table <- man_table[,c(1,4:7,10,12)]
names(man_table)[1:5] <- c("StudyPatientNo_","Eye.Opening","Best.Verbal.Response","Best.Motor.Response","Glasgow.Coma.Scale.Score")
accel_no_match <- patient_clinical_data %>% select(StudyPatientNo_,AccelPatientNo_)

comb_man_gcs_table <- left_join(man_table,accel_no_match,by = "StudyPatientNo_") %>% drop_na(AccelPatientNo_)
final_man_gcs_table <- comb_man_gcs_table[,c(1,8,7,5,4,3,2,6)]
final_man_gcs_table <- final_man_gcs_table[order(final_man_gcs_table$TakenInstant),]
write.csv(final_man_gcs_table,'../clinical_data/clean_manu_GCS_table.csv')

# Clear all existing variables from environment to begin next step:
rm(list = ls())
gc()

# CHECKPOINT 1: beginning of GCS label resampling

# Load clinical patient data:
source('./functions/load_patient_clinical_data.R')
patient_clinical_data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv') %>% arrange(AccelPatientNo_) %>% mutate(ptIdx = 1:nrow(.))

# Load missing accelerometry information:
missingTimeInfo <- read_xlsx('~/scratch/all_motion_feature_data/MissingPercentTable.xlsx',.name_repair = "universal") %>% mutate(Start.Timestamp = as.POSIXct(Start.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"), End.Timestamp = as.POSIXct(End.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"))

# Load automatically extracted GCS labels:
gcs_data <- read.csv('../clinical_data/clean_auto_GCS_table.csv') %>% select(-X) %>% mutate(TakenInstant = as.POSIXct(TakenInstant, tz = "America/New_York"))

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

dir.create('~/scratch/all_motion_feature_data/gcs_labels',showWarnings = FALSE)

# a) Detection models:

det_parameters <- as.data.frame(matrix(nrow = length(det_obs_windows),ncol = 0))
det_gcs_labels <- vector(mode = "list")
for (i in 1:length(det_obs_windows)) {
  det_parameters$obs_windows[i] <- det_obs_windows[i]
  print(paste("Observation Window Size:",det_obs_windows[i],"hours started"))
  curr_secs_window <- det_secs_windows[i]
  gcsLabels <- as.data.frame(matrix(nrow = 0,ncol = 8))
  for (patIdx in unique(patient_clinical_data$ptIdx)) {
    print(paste("Patient no.",patIdx,"started"))
    currStudyNo <- patient_clinical_data$StudyPatientNo_[patIdx]
    currAccelNo <- patient_clinical_data$AccelPatientNo_[patIdx]
    currStartTm <- missingTimeInfo$Start.Timestamp[missingTimeInfo$Accel.Patient.No. == currAccelNo]
    currEndTm <- missingTimeInfo$End.Timestamp[missingTimeInfo$Accel.Patient.No. == currAccelNo]
    currGCS <- gcs_data %>% filter(StudyPatientNo_ == currStudyNo, AccelPatientNo_ == currAccelNo)
    gcsFilter <- currGCS$TakenInstant >= currStartTm+curr_secs_window & currGCS$TakenInstant <= currEndTm
    gcsLabels <- rbind(gcsLabels,currGCS[gcsFilter,])
  }
  det_gcs_labels[[i]] <- gcsLabels
}
save(det_parameters, det_gcs_labels, file = '~/scratch/all_motion_feature_data/gcs_labels/detection_labels.RData')

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
    for (patIdx in unique(patient_clinical_data$ptIdx)) {
      print(paste("Patient no.",patIdx,"started"))
      currStudyNo <- patient_clinical_data$StudyPatientNo_[patIdx]
      currAccelNo <- patient_clinical_data$AccelPatientNo_[patIdx]
      currStartTm <- missingTimeInfo$Start.Timestamp[missingTimeInfo$Accel.Patient.No. == currAccelNo]
      currEndTm <- missingTimeInfo$End.Timestamp[missingTimeInfo$Accel.Patient.No. == currAccelNo]
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
save(pre_parameters, pre_gcs_labels, file = '~/scratch/all_motion_feature_data/gcs_labels/prediction_labels.RData')

# Partition data into training, validation, and testing sets based on preserving class imbalance:
rm(list = ls())
gc()

p_train <- 0.8 # set proportion of observations in training (including validation) sets

# (a) Detection GCS labels:

load('~/scratch/all_motion_feature_data/gcs_labels/detection_labels.RData')

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
save(det_motor_train_idx, det_eye_train_idx,det_motor_test_idx,det_eye_test_idx, file = '~/scratch/all_motion_feature_data/gcs_labels/detection_partitions.RData')

# (b) Prediction GCS labels:

load('~/scratch/all_motion_feature_data/gcs_labels/prediction_labels.RData')

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
save(pre_motor_train_idx, pre_eye_train_idx,pre_motor_test_idx,pre_eye_test_idx, file = '~/scratch/all_motion_feature_data/gcs_labels/prediction_partitions.RData')

# Examnine GCS scores of each patient 
# Load automatically extracted GCS labels:
gcs_data <- read.csv('../clinical_data/clean_auto_GCS_table.csv') %>% select(-X) %>% mutate(TakenInstant = as.POSIXct(TakenInstant, tz = "America/New_York"))

# Merge GCS.data and patient clinical data
merged.gcs.data <- left_join(gcs_data, patient_clinical_data, by = "AccelPatientNo_") %>% mutate(TakenDay = as.Date(TakenInstant, tz = "America/New_York")) %>% arrange(AccelPatientNo_,TakenInstant)

## Steps:
#-calculate number of evals per day per patient
summ.gcs.stats <- merged.gcs.data %>% filter(TakenDay >= NCCUAdmissionDate & TakenDay <= NCCUDischargeDate) %>% group_by(AccelPatientNo_, TakenDay) %>% summarise(no.evals = n())
#-filter out GCS scores that coincide with accelerometry recording time
accel.merged.gcs.data <- left_join(gcs_data, missingTimeInfo, by = "AccelPatientNo_") %>% arrange(AccelPatientNo_,TakenInstant) %>% filter(TakenInstant >= Start.Timestamp & TakenInstant <= End.Timestamp)
#-worst GCSm within 24 hours of admission
worst.gcsm.nccu.admission <- merged.gcs.data %>% filter(TakenDay >= NCCUAdmissionDate & TakenDay <= NCCUAdmissionDate+1) %>% group_by(AccelPatientNo_) %>% summarise(worst.GCSm = min(Best.Motor.Response,na.rm = TRUE))