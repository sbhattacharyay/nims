#### Master Script 4: Resampling of GCS data for classification ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# University of Cambridge
# Johns Hopkins University
# email address: sb2406@cam.ac.uk
#
### Contents:
# I. Initialization
# II. Isolate and save candidate (coincide with accelerometry) GCS observations
# III. Partition data with stratified, repeated cross-validation

### I. Initialization
## Import necessary packages
library(tidyverse)
library(readxl)
library(caret)

## Load neurological assessment data
gcs.data <- read.csv('../clinical_data/neurological_assessments.csv')

## Load temporal clinical metadata
patient.temporal.data <- read.csv('../clinical_data/patient_temporal_info.csv')

## Load patient outcome information for GOSE-stratification
patient.outcomes <-read.csv('../clinical_data/patient_outcomes.csv')

## Set accelerometry observation window parameters
obs.windows.hours <- c(seq(0.05,0.5,by = .05),1:3,seq(6,18,by = 3),24) #in hours

## Set parameters for repeated cross-validation
NUM.REPEATS <- 5
NUM.FOLDS <- 5

### II. Isolate and save candidate (coincide with accelerometry) GCS observations:
obs.gcs.data <- vector(mode = "list")
for (curr.obs.window in obs.windows.hours){
  curr.obs.window.df <- as.data.frame(matrix(nrow = 0,ncol=8))
  for (curr.UPI in patient.temporal.data$UPI){
    curr.UPI.start.hour <- patient.temporal.data$HoursFromICUAdmissionAccelRecordingStart[patient.temporal.data$UPI == curr.UPI]
    curr.UPI.end.hour <- patient.temporal.data$HoursFromICUAdmissionAccelRecordingEnd[patient.temporal.data$UPI == curr.UPI]
    curr.gcs.data <- gcs.data %>%
      filter(UPI == curr.UPI,
             HoursFromICUAdmission <= curr.UPI.end.hour,
             HoursFromICUAdmission >= curr.UPI.start.hour + curr.obs.window)
    curr.obs.window.df <- rbind(curr.obs.window.df,curr.gcs.data)
  }
  curr.obs.window.df <- curr.obs.window.df %>% drop_na(GCSm)
  obs.gcs.data[[sprintf('%05.2f',curr.obs.window)]] <- curr.obs.window.df
}

### III. Partition data with stratified, repeated cross-validation
dir.create('../validation_sampling',showWarnings = F)
for (curr.obs.idx in 1:length(obs.gcs.data)){
  # Get median GCSm of each patient and add GOSE score at discharge
  grouped.curr.obs.window.df <- obs.gcs.data[[curr.obs.idx]] %>%
    group_by(UPI) %>%
    summarise(GCSm = floor(median(GCSm))) %>%
    left_join(patient.outcomes %>% select(UPI,GOSEDischarge,GOSE12Months),by='UPI')
  # Use `caret` package to create stratified repeated cross-validation folds
  curr.obs.window.GCSm.fold.idx <- createMultiFolds(grouped.curr.obs.window.df$GCSm, k = NUM.FOLDS, times = NUM.REPEATS)
  curr.obs.window.GOSE.fold.idx <- createMultiFolds(grouped.curr.obs.window.df$GOSEDischarge, k = NUM.FOLDS, times = NUM.REPEATS)
  
  # Isolate cases within the current observation window with non-missing GOSE outcomes at 12 months
  grouped.12m.curr.obs.window.df <- grouped.curr.obs.window.df %>% drop_na(GOSE12Months)
  
  # Use `caret` package to create stratified repeated cross-validation folds for 12-month outcomes
  curr.obs.window.GOSE.12m.fold.idx <- createMultiFolds(grouped.12m.curr.obs.window.df$GOSE12Months, k = NUM.FOLDS, times = NUM.REPEATS)
  
  # Create and store repeated cross-validation splits for GCSm detection
  curr.obs.window.GCSm.fold.df <- as.data.frame(matrix(ncol = 11,nrow = 0))
  for (fold.idx in 1:length(curr.obs.window.GCSm.fold.idx)){
    curr.fold <- as.integer(str_match(names(curr.obs.window.GCSm.fold.idx)[fold.idx], "Fold\\s*(.*?)\\s*.Rep")[,2])
    curr.repeat <- as.integer(sub(".*Rep", "", names(curr.obs.window.GCSm.fold.idx)[fold.idx]))
    train.UPI <- grouped.curr.obs.window.df$UPI[curr.obs.window.GCSm.fold.idx[[fold.idx]]]
    train.df <- obs.gcs.data[[curr.obs.idx]] %>% 
      filter(UPI %in% train.UPI) %>% 
      mutate(Fold = curr.fold,
             Repeat = curr.repeat,
             Split = 'Train')
    test.df <- obs.gcs.data[[curr.obs.idx]] %>% 
      filter(!UPI %in% train.UPI) %>% 
      mutate(Fold = curr.fold,
             Repeat = curr.repeat,
             Split = 'Test')
    curr.obs.window.GCSm.fold.df <- rbind(curr.obs.window.GCSm.fold.df,train.df,test.df)
  }
  # Get current observation window length and save current GCSm splits information
  curr.obs.window.hours <- names(obs.gcs.data)[curr.obs.idx]
  write.csv(curr.obs.window.GCSm.fold.df,paste0('../validation_sampling/',curr.obs.window.hours,'_h_GCSm_folds.csv'),row.names = F)
  
  # Create and store repeated cross-validation splits for GOSE prediction at discharge
  curr.obs.window.GOSE.fold.df <- as.data.frame(matrix(ncol = 11,nrow = 0))
  for (fold.idx in 1:length(curr.obs.window.GOSE.fold.idx)){
    curr.fold <- as.integer(str_match(names(curr.obs.window.GOSE.fold.idx)[fold.idx], "Fold\\s*(.*?)\\s*.Rep")[,2])
    curr.repeat <- as.integer(sub(".*Rep", "", names(curr.obs.window.GOSE.fold.idx)[fold.idx]))
    train.UPI <- grouped.curr.obs.window.df$UPI[curr.obs.window.GOSE.fold.idx[[fold.idx]]]
    train.df <- obs.gcs.data[[curr.obs.idx]] %>% 
      filter(UPI %in% train.UPI) %>% 
      mutate(Fold = curr.fold,
             Repeat = curr.repeat,
             Split = 'Train')
    test.df <- obs.gcs.data[[curr.obs.idx]] %>% 
      filter(!UPI %in% train.UPI) %>% 
      mutate(Fold = curr.fold,
             Repeat = curr.repeat,
             Split = 'Test')
    curr.obs.window.GOSE.fold.df <- rbind(curr.obs.window.GOSE.fold.df,train.df,test.df)
  }
  # Get current observation window length and save current GOSE splits information
  curr.obs.window.hours <- names(obs.gcs.data)[curr.obs.idx]
  write.csv(curr.obs.window.GOSE.fold.df,paste0('../validation_sampling/',curr.obs.window.hours,'_h_GOSE_folds.csv'),row.names = F)
  
  # Create and store repeated cross-validation splits for GOSE prediction at 12 months
  curr.obs.window.GOSE.12m.fold.df <- as.data.frame(matrix(ncol = 11,nrow = 0))
  for (fold.idx in 1:length(curr.obs.window.GOSE.12m.fold.idx)){
    curr.fold <- as.integer(str_match(names(curr.obs.window.GOSE.12m.fold.idx)[fold.idx], "Fold\\s*(.*?)\\s*.Rep")[,2])
    curr.repeat <- as.integer(sub(".*Rep", "", names(curr.obs.window.GOSE.12m.fold.idx)[fold.idx]))
    train.UPI <- grouped.12m.curr.obs.window.df$UPI[curr.obs.window.GOSE.12m.fold.idx[[fold.idx]]]
    test.UPI <- grouped.12m.curr.obs.window.df$UPI[-curr.obs.window.GOSE.12m.fold.idx[[fold.idx]]]
    train.df <- obs.gcs.data[[curr.obs.idx]] %>% 
      filter(UPI %in% train.UPI) %>% 
      mutate(Fold = curr.fold,
             Repeat = curr.repeat,
             Split = 'Train')
    test.df <- obs.gcs.data[[curr.obs.idx]] %>% 
      filter(UPI %in% test.UPI) %>% 
      mutate(Fold = curr.fold,
             Repeat = curr.repeat,
             Split = 'Test')
    curr.obs.window.GOSE.12m.fold.df <- rbind(curr.obs.window.GOSE.12m.fold.df,train.df,test.df)
  }
  # Get current observation window length and save current GOSE splits information
  curr.obs.window.hours <- names(obs.gcs.data)[curr.obs.idx]
  write.csv(curr.obs.window.GOSE.12m.fold.df,paste0('../validation_sampling/',curr.obs.window.hours,'_h_GOSE12m_folds.csv'),row.names = F)
}