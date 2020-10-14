#### Temporary Master Script ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# Department of Biomedical Engineering
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

# Add library path for R packages and set working directory:
.libPaths(c("~/Rpackages" , .libPaths()))
setwd('~/work/nims/scripts')

# Load necessary libraries
library(tidyverse)
library(R.matlab)
library(readxl)
library(caret)

# Load clinical patient data:
source('./functions/load_patient_clinical_data.R')
patient_clinical_data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv') %>% arrange(AccelPatientNo_) %>% mutate(ptIdx = 1:nrow(.))

# Load missing accelerometry information:
missingTimeInfo <- read_xlsx('~/data/all_motion_feature_data/MissingPercentTable.xlsx',.name_repair = "universal") %>% mutate(Start.Timestamp = as.POSIXct(Start.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"), End.Timestamp = as.POSIXct(End.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"))

# Load GCS labels:
gcs_data <- read.csv('../clinical_data/clean_auto_GCS_table.csv') %>% select(-X) %>% mutate(TakenInstant = as.POSIXct(TakenInstant, tz = "America/New_York"))

# Load imputation 1:
load('~/scratch/all_motion_feature_data/final_imputed_features/imp1.RData')
imp1 <- currImp
rm(currImp)
gc()

# Set observation window time:
obs_window <- 1 # in hours
secs_window <- obs_window*3600
no_pts <- obs_window*60*12 + 1

# Isolate and save candidate GCS observations
gcsLabels <- as.data.frame(matrix(nrow = 0,ncol = 8))
for (patIdx in unique(imp1$ptIdx)) {
  print(paste("Patient no.",patIdx,"started"))
  currStudyNo <- patient_clinical_data$StudyPatientNo_[patIdx]
  currAccelNo <- patient_clinical_data$AccelPatientNo_[patIdx]
  currStartTm <- missingTimeInfo$Start.Timestamp[missingTimeInfo$Accel.Patient.No. == currAccelNo]
  currEndTm <- missingTimeInfo$End.Timestamp[missingTimeInfo$Accel.Patient.No. == currAccelNo]
  currGCS <- gcs_data %>% filter(StudyPatientNo_ == currStudyNo, AccelPatientNo_ == currAccelNo)
  gcsFilter <- currGCS$TakenInstant >= currStartTm+secs_window & currGCS$TakenInstant <= currEndTm
  gcsLabels <- rbind(gcsLabels,currGCS[gcsFilter,])
}
rm(gcs_data,currGCS,currStartTm,currEndTm)
gc()
write_csv(gcsLabels,'~/data/all_motion_feature_data/gcs_labels.csv')

# Partition data into training (including validation) and testing sets:
set.seed(2020)
GCSe <- gcsLabels %>% drop_na(Eye.Opening)
GCSm <- gcsLabels %>% drop_na(Best.Motor.Response)
GCSv <- gcsLabels %>% drop_na(Best.Verbal.Response)
GCST <- gcsLabels %>% drop_na(Glasgow.Coma.Scale.Score)

p_train <- 0.8 # proportion of observations in training (including validation) sets
GCSe_trainIdx <- createDataPartition(GCSe$Eye.Opening,p=p_train,list = FALSE,times = 1)
GCSm_trainIdx <- createDataPartition(GCSm$Best.Motor.Response,p=p_train,list = FALSE,times = 1)
GCSv_trainIdx <- createDataPartition(GCSv$Best.Verbal.Response,p=p_train,list = FALSE,times = 1)
GCST_trainIdx <- createDataPartition(GCST$Glasgow.Coma.Scale.Score,p=p_train,list = FALSE,times = 1)

save(GCSe_trainIdx,GCSm_trainIdx,GCSv_trainIdx,GCST_trainIdx,file="~/data/all_motion_feature_data/train_test_sets/train_idxs.RData")

GCSe_train <- GCSe[GCSe_trainIdx,]
GCSm_train <- GCSm[GCSm_trainIdx,]
GCSv_train <- GCSv[GCSv_trainIdx,]
GCST_train <- GCST[GCST_trainIdx,]

GCSe_test <- GCSe[-GCSe_trainIdx,]
GCSm_test <- GCSm[-GCSm_trainIdx,]
GCSv_test <- GCSv[-GCSv_trainIdx,]
GCST_test <- GCST[-GCST_trainIdx,]

# Feature combination names
featureTypes <- unique(imp1$featureType)
sensorLocs <- c("LE","LW","LA","RE","RW","RA")
dataPtIdx <- 1:no_pts
nameDF <- expand.grid(featureTypes, sensorLocs,dataPtIdx) %>% arrange(Var2, Var1, Var3)
nameVector <- paste0(nameDF$Var2,".",nameDF$Var1,".",nameDF$Var3)
no_channels <- length(sensorLocs)*length(featureTypes)

# Save training and testing data from each imputation into NumPy compatible arrays:
np <- import("numpy")

GCSe_trainLabels <- GCSe_train$Eye.Opening
GCSm_trainLabels <- GCSm_train$Best.Motor.Response
GCSv_trainLabels <- GCSv_train$Best.Verbal.Response
GCST_trainLabels <- GCST_train$Glasgow.Coma.Scale.Score

GCSe_testLabels <- GCSe_test$Eye.Opening
GCSm_testLabels <- GCSm_test$Best.Motor.Response
GCSv_testLabels <- GCSv_test$Best.Verbal.Response
GCST_testLabels <- GCST_test$Glasgow.Coma.Scale.Score

for (impNo in 1:length(list.files('~/data/all_motion_feature_data/final_imputed_features/'))){
  print(paste("Imputation no.",impNo,"started"))
  load(file.path('~/data/all_motion_feature_data/final_imputed_features',paste0('imp',impNo,'.RData')))
  
  GCSe_trainArray <- array(dim=c(length(GCSe_trainIdx),1,no_pts,no_channels))
  GCSm_trainArray <- array(dim=c(length(GCSm_trainIdx),1,no_pts,no_channels))
  GCSv_trainArray <- array(dim=c(length(GCSv_trainIdx),1,no_pts,no_channels))
  GCST_trainArray <- array(dim=c(length(GCST_trainIdx),1,no_pts,no_channels))
  
  GCSe_testArray <- array(dim=c(nrow(GCSe_test),1,no_pts,no_channels))
  GCSm_testArray <- array(dim=c(nrow(GCSm_test),1,no_pts,no_channels))
  GCSv_testArray <- array(dim=c(nrow(GCSv_test),1,no_pts,no_channels))
  GCST_testArray <- array(dim=c(nrow(GCST_test),1,no_pts,no_channels))
  
  # Eye Training
  for (i in 1:nrow(GCSe_trainIdx)){
    print(paste("Case no.",i,"started"))
    currTimeStamp <- GCSe_train$TakenInstant[i]
    currFeatSet <- currImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCSe_train$AccelPatientNo_[i], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
    tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
    tempReshape <- as.matrix(tempReshape[,-1])
    GCSe_trainArray[i,1,,] <- tempReshape
  }
  
  # Motor Training
  for (i in 1:nrow(GCSm_trainIdx)){
    print(paste("Case no.",i,"started"))
    currTimeStamp <- GCSm_train$TakenInstant[i]
    currFeatSet <- currImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCSm_train$AccelPatientNo_[i], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
    tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
    tempReshape <- as.matrix(tempReshape[,-1])
    GCSm_trainArray[i,1,,] <- tempReshape
  }
  
  # Verbal Training
  for (i in 1:nrow(GCSv_trainIdx)){
    print(paste("Case no.",i,"started"))
    currTimeStamp <- GCSv_train$TakenInstant[i]
    currFeatSet <- currImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCSv_train$AccelPatientNo_[i], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
    tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
    tempReshape <- as.matrix(tempReshape[,-1])
    GCSv_trainArray[i,1,,] <- tempReshape
  }
  
  # Total Training
  for (i in 1:nrow(GCST_trainIdx)){
    print(paste("Case no.",i,"started"))
    currTimeStamp <- GCST_train$TakenInstant[i]
    currFeatSet <- currImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCST_train$AccelPatientNo_[i], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
    tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
    tempReshape <- as.matrix(tempReshape[,-1])
    GCST_trainArray[i,1,,] <- tempReshape
  }
  
  # Eye Testing
  for (i in 1:nrow(GCSe_test)){
    print(paste("Case no.",i,"started"))
    currTimeStamp <- GCSe_test$TakenInstant[i]
    currFeatSet <- currImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCSe_test$AccelPatientNo_[i], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
    tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
    tempReshape <- as.matrix(tempReshape[,-1])
    GCSe_testArray[i,1,,] <- tempReshape
  }
  
  # Motor Testing
  for (i in 1:nrow(GCSm_test)){
    print(paste("Case no.",i,"started"))
    currTimeStamp <- GCSm_test$TakenInstant[i]
    currFeatSet <- currImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCSm_test$AccelPatientNo_[i], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
    tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
    tempReshape <- as.matrix(tempReshape[,-1])
    GCSm_testArray[i,1,,] <- tempReshape
  }
  
  # Verbal Testing
  for (i in 1:nrow(GCSv_test)){
    print(paste("Case no.",i,"started"))
    currTimeStamp <- GCSv_test$TakenInstant[i]
    currFeatSet <- currImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCSv_test$AccelPatientNo_[i], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
    tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
    tempReshape <- as.matrix(tempReshape[,-1])
    GCSv_testArray[i,1,,] <- tempReshape
  }
  
  # Total Testing
  for (i in 1:nrow(GCST_test)){
    print(paste("Case no.",i,"started"))
    currTimeStamp <- GCST_test$TakenInstant[i]
    currFeatSet <- currImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCST_test$AccelPatientNo_[i], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
    tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
    tempReshape <- as.matrix(tempReshape[,-1])
    GCST_testArray[i,1,,] <- tempReshape
  }
  
  # Save all training and testing sets into a single RData file:
  fileName <- file.path("~/data/all_motion_feature_data/train_test_sets",paste0('imp',impNo,'_test_train.RData'))
  save(GCSe_trainArray,GCSm_trainArray,GCSv_trainArray,GCST_trainArray,GCSe_testArray,GCSm_testArray,GCSv_testArray,GCST_testArray,GCSe_trainLabels,GCSm_trainLabels,GCSv_trainLabels,GCST_trainLabels,GCSe_testLabels,GCSm_testLabels,GCSv_testLabels,GCST_testLabels,file=fileName)
  
  # Save training and testing set arrays into NumPy files
  baseDir <- path.expand('~/data/all_motion_feature_data/train_test_sets')
  dir.create(file.path(baseDir,paste0('imp',impNo)))
  
  np$save(file.path(baseDir,paste0('imp',impNo),'GCSe_trainArray'),GCSe_trainArray)
  np$save(file.path(baseDir,paste0('imp',impNo),'GCSm_trainArray'),GCSm_trainArray)
  np$save(file.path(baseDir,paste0('imp',impNo),'GCSv_trainArray'),GCSv_trainArray)
  np$save(file.path(baseDir,paste0('imp',impNo),'GCST_trainArray'),GCST_trainArray)
  
  np$save(file.path(baseDir,paste0('imp',impNo),'GCSe_testArray'),GCSe_testArray)
  np$save(file.path(baseDir,paste0('imp',impNo),'GCSm_testArray'),GCSm_testArray)
  np$save(file.path(baseDir,paste0('imp',impNo),'GCSv_testArray'),GCSv_testArray)
  np$save(file.path(baseDir,paste0('imp',impNo),'GCST_testArray'),GCST_testArray)
  
  np$save(file.path(baseDir,paste0('imp',impNo),'GCSe_trainLabels'),GCSe_trainLabels)
  np$save(file.path(baseDir,paste0('imp',impNo),'GCSm_trainLabels'),GCSm_trainLabels)
  np$save(file.path(baseDir,paste0('imp',impNo),'GCSv_trainLabels'),GCSv_trainLabels)
  np$save(file.path(baseDir,paste0('imp',impNo),'GCST_trainLabels'),GCST_trainLabels)
  
  np$save(file.path(baseDir,paste0('imp',impNo),'GCSe_testLabels'),GCSe_testLabels)
  np$save(file.path(baseDir,paste0('imp',impNo),'GCSm_testLabels'),GCSm_testLabels)
  np$save(file.path(baseDir,paste0('imp',impNo),'GCSv_testLabels'),GCSv_testLabels)
  np$save(file.path(baseDir,paste0('imp',impNo),'GCST_testLabels'),GCST_testLabels)
}

############################################################################################################################################################################################################

# GCSe <- gcsLabels %>% drop_na(Eye.Opening) %>% group_by(AccelPatientNo_) %>% summarise(medScore = floor(median(Eye.Opening)), stdScores = sd(Eye.Opening), n = n())
# GCSm <- gcsLabels %>%drop_na(Best.Motor.Response) %>% group_by(AccelPatientNo_) %>% summarise(medScore = floor(median(Best.Motor.Response)), stdScores = sd(Best.Motor.Response), n = n())
# GCSv<- gcsLabels %>%drop_na(Best.Verbal.Response) %>%group_by(AccelPatientNo_) %>% summarise(medScore = floor(median(Best.Verbal.Response)), stdScores = sd(Best.Verbal.Response), n = n())
# GCST <- gcsLabels %>%drop_na(Glasgow.Coma.Scale.Score) %>%group_by(AccelPatientNo_) %>% summarise(medScore = floor(median(Glasgow.Coma.Scale.Score)), stdScores = sd(Glasgow.Coma.Scale.Score), n = n())
# 
# set.seed(2020)
# GCSe_trainIdx <- createDataPartition(GCSe$medScore,p=.8,list = FALSE,times = 1)
# GCSm_trainIdx <- createDataPartition(GCSm$medScore,p=.8,list = FALSE,times = 1)
# GCSv_trainIdx <- createDataPartition(GCSv$medScore,p=.8,list = FALSE,times = 1)
# GCST_trainIdx <- createDataPartition(GCST$medScore,p=.8,list = FALSE,times = 1)
# 
# GCSe_trainLabels <- gcsLabels[gcsLabels$StudyPatientNo_ %in% GCSe_trainIdx,]
# GCSm_trainLabels <- gcsLabels[gcsLabels$StudyPatientNo_ %in% GCSm_trainIdx,]
# GCSv_trainLabels <- gcsLabels[gcsLabels$StudyPatientNo_ %in% GCSv_trainIdx,]
# GCST_trainLabels <- gcsLabels[gcsLabels$StudyPatientNo_ %in% GCST_trainIdx,]
# 
# GCSe_testLabels <- gcsLabels[!gcsLabels$StudyPatientNo_ %in% GCSe_trainIdx,]
# GCSm_testLabels <- gcsLabels[!gcsLabels$StudyPatientNo_ %in% GCSm_trainIdx,]
# GCSv_testLabels <- gcsLabels[!gcsLabels$StudyPatientNo_ %in% GCSv_trainIdx,]
# GCST_testLabels <- gcsLabels[!gcsLabels$StudyPatientNo_ %in% GCST_trainIdx,]

# Prepare training set
for (j in 1:nrow(GCSe_trainLabels)){
  print(paste("Case no.",j,"started"))
  currTimeStamp <- GCSe_trainLabels$TakenInstant[j]
  currFeatSet <- imp1 %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCSe_trainLabels$AccelPatientNo_[j], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
  tempVector <- as.matrix(currFeatSet[,c(2,3,1,5,6,4)]) %>% as.vector()
  compiledTrainingMatrixGCSe <- rbind(compiledTrainingMatrixGCSe, tempVector)
  if (is_empty(tempVector)){print(paste("Case no.",j,"empty"))}
}

names(compiledTrainingMatrixGCSe) <- nameVector

compiledTestingMatrixGCSe <- as.data.frame(matrix(nrow = 0,ncol = length(nameVector)))
for (j in 1:nrow(GCSe_testLabels)){
  print(paste("Case no.",j,"started"))
  currTimeStamp <- GCSe_testLabels$TakenInstant[j]
  currFeatSet <- imp1 %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCSe_testLabels$AccelPatientNo_[j], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
  tempVector <- as.matrix(currFeatSet[,c(2,3,1,5,6,4)]) %>% as.vector()
  compiledTestingMatrixGCSe <- rbind(compiledTestingMatrixGCSe, tempVector)
  if (is_empty(tempVector)){print(paste("Case no.",j,"empty"))}
}
