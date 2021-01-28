#### Master Script 1: Clinical Data Extraction and Examination ####
#
# Eshan Joshi, Shubhayu Bhattacharyay
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

library(tidyverse)

# Load patient clinical data
source('./functions/load_patient_clinical_data.R')
patient.clinical.data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv') %>% arrange(AccelPatientNo_) %>% mutate(ptIdx = 1:nrow(.))

# Load automatically extracted GCS labels:
gcs.data <- read.csv('../clinical_data/clean_auto_GCS_table.csv') %>% 
    select(-X, -Coincide.with.Accel.Recording) %>% 
    mutate(TakenInstant = as.POSIXct(TakenInstant, tz = "America/New_York")) %>% 
    mutate(TakenDay = as.Date(TakenInstant, tz = "America/New_York"),
           det.indicator = FALSE,
           pre.indicator.lead.0  = FALSE,
           pre.indicator.lead.1  = FALSE,
           pre.indicator.lead.2  = FALSE,
           pre.indicator.lead.6  = FALSE) %>% 
    arrange(AccelPatientNo_,TakenInstant)

# Load missing time info (includes information on recording start and stop times)
missing.time.info <- read_xlsx('../all_motion_feature_data/MissingPercentTable.xlsx',.name_repair = "universal") %>% 
    mutate(Start.Timestamp = as.POSIXct(Start.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"), End.Timestamp = as.POSIXct(End.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York")) %>% 
    mutate(Recording.Duration = End.Timestamp - Start.Timestamp) %>%
    rename(AccelPatientNo_=Accel.Patient.No.) %>%
    mutate(Bed.miss.hours = Bed*Recording.Duration,LA.miss.hours = LA*Recording.Duration,LE.miss.hours = LE*Recording.Duration,LW.miss.hours = LW*Recording.Duration,RA.miss.hours = RA*Recording.Duration,RE.miss.hours = RE*Recording.Duration,RW.miss.hours = RW*Recording.Duration)

# Find GCS data that corresponds to detection and prediction visualization exercises 
for (i in 1:nrow(missing.time.info)){
    curr.start.time <- missing.time.info$Start.Timestamp[i]
    curr.end.time <- missing.time.info$End.Timestamp[i]
    
    keep.detection.idx <- which(gcs.data$AccelPatientNo_ == missing.time.info$AccelPatientNo_[i] & (gcs.data$TakenInstant>=curr.start.time & gcs.data$TakenInstant<=curr.end.time+(6*60*60)))
    keep.prediction.lead.1.idx <- which(gcs.data$AccelPatientNo_ == missing.time.info$AccelPatientNo_[i] & (gcs.data$TakenInstant>=curr.start.time+(1*60*60) & gcs.data$TakenInstant<=curr.end.time+(7*60*60)))
    keep.prediction.lead.2.idx <- which(gcs.data$AccelPatientNo_ == missing.time.info$AccelPatientNo_[i] & (gcs.data$TakenInstant>=curr.start.time+(2*60*60) & gcs.data$TakenInstant<=curr.end.time+(8*60*60)))
    keep.prediction.lead.6.idx <- which(gcs.data$AccelPatientNo_ == missing.time.info$AccelPatientNo_[i] & (gcs.data$TakenInstant>=curr.start.time+(6*60*60) & gcs.data$TakenInstant<=curr.end.time+(12*60*60)))
    
    gcs.data$det.indicator[keep.detection.idx] <- TRUE
    gcs.data$pre.indicator.lead.0[keep.detection.idx] <- TRUE
    gcs.data$pre.indicator.lead.1[keep.prediction.lead.1.idx] <- TRUE
    gcs.data$pre.indicator.lead.2[keep.prediction.lead.2.idx] <- TRUE
    gcs.data$pre.indicator.lead.6[keep.prediction.lead.6.idx] <- TRUE
}

# Save GCS file with indicator variables
write.csv(gcs.data,'../clinical_data/gcs_data_w_indicators.csv')