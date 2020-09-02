library(tidyverse)
library(R.matlab)
library(readxl)
library(caret)
library(lolR)
source('./functions/load_patient_clinical_data.R')

patient_clinical_data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv') %>% arrange(AccelPatientNo_) 
patient_clinical_data <- patient_clinical_data %>% mutate(ptIdx = 1:nrow(patient_clinical_data))
gcs_data <- read.csv('../clinical_data/clean_auto_GCS_table.csv') %>% select(-X) %>% mutate(TakenInstant = as.POSIXct(TakenInstant, tz = "America/New_York"))
featureSet <- readRDS('../all_motion_feature_data/bc_i_c_dataset.rds')
indexedTimes <- read.csv('../all_motion_feature_data/indexed_times.csv') %>% mutate(times = as.POSIXct(times, format = "%d-%b-%Y %H:%M:%S",tz = "America/New_York"))
firstImp <- featureSet[[1]] %>% select(-Bed, -ImputNum)
firstImp$TakenInstant <- ""
missingTimeInfo <- read_xlsx('../all_motion_feature_data/MissingPercentTable.xlsx',.name_repair = "universal") %>% mutate(Start.Timestamp = as.POSIXct(Start.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"), End.Timestamp = as.POSIXct(End.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"))
rm(featureSet)
gc()

# Add times to motion feature data:
for (patIdx in unique(indexedTimes$ptIdx)){
  currPtIdx <- firstImp$ptIdx == patIdx
  currTimes <- indexedTimes %>% filter(ptIdx == patIdx) %>% rename(TakenInstant = times) %>% mutate(timeCount = 1:(sum(currPtIdx)/7)) %>% slice(rep(1:n(), times = 7))
  currTimes$TakenInstant <- as.character(currTimes$TakenInstant)
  firstImp[currPtIdx,c("TakenInstant","ptIdx","timeCount")] <- currTimes
}
firstImp$TakenInstant <- as.POSIXct(firstImp$TakenInstant, tz = "America/New_York")
rm(currPtIdx)
rm(currTimes)
rm(indexedTimes)
gc()

save(firstImp,file = '../all_motion_feature_data/firstImp.RData')

load('../all_motion_feature_data/firstImp.RData')

time_preceding <- 1 # in hours
secs_window <- time_preceding*3600
no_pts = time_preceding*60*12 + 1;

gcsLabels <- as.data.frame(matrix(nrow = 0,ncol = 8))
  
for (patIdx in unique(firstImp$ptIdx)) {
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

# Create Data Partitions based on GCS labels, ensuring the same patient is not in both training and testing

GCSe <- gcsLabels %>% drop_na(Eye.Opening) %>% group_by(AccelPatientNo_) %>% summarise(medScore = floor(median(Eye.Opening)), n = n())
GCSm <- gcsLabels %>%drop_na(Best.Motor.Response) %>% group_by(AccelPatientNo_) %>% summarise(medScore = floor(median(Best.Motor.Response)), n = n())
GCSv<- gcsLabels %>%drop_na(Best.Verbal.Response) %>%group_by(AccelPatientNo_) %>% summarise(medScore = floor(median(Best.Verbal.Response)), n = n())
GCST <- gcsLabels %>%drop_na(Glasgow.Coma.Scale.Score) %>%group_by(AccelPatientNo_) %>% summarise(medScore = floor(median(Glasgow.Coma.Scale.Score)), n = n())

set.seed(2020)
GCSe_trainIdx <- createDataPartition(GCSe$medScore,p=.8,list = FALSE,times = 1)
GCSm_trainIdx <- createDataPartition(GCSm$medScore,p=.8,list = FALSE,times = 1)
GCSv_trainIdx <- createDataPartition(GCSv$medScore,p=.8,list = FALSE,times = 1)
GCST_trainIdx <- createDataPartition(GCST$medScore,p=.8,list = FALSE,times = 1)

GCSe_trainLabels <- gcsLabels[gcsLabels$StudyPatientNo_ %in% GCSe_trainIdx,]
GCSm_trainLabels <- gcsLabels[gcsLabels$StudyPatientNo_ %in% GCSm_trainIdx,]
GCSv_trainLabels <- gcsLabels[gcsLabels$StudyPatientNo_ %in% GCSv_trainIdx,]
GCST_trainLabels <- gcsLabels[gcsLabels$StudyPatientNo_ %in% GCST_trainIdx,]

GCSe_testLabels <- gcsLabels[!gcsLabels$StudyPatientNo_ %in% GCSe_trainIdx,]
GCSm_testLabels <- gcsLabels[!gcsLabels$StudyPatientNo_ %in% GCSm_trainIdx,]
GCSv_testLabels <- gcsLabels[!gcsLabels$StudyPatientNo_ %in% GCSv_trainIdx,]
GCST_testLabels <- gcsLabels[!gcsLabels$StudyPatientNo_ %in% GCST_trainIdx,]

featureTypes <- unique(firstImp$featureType)
sensorLocs <- c("LE","LW","LA","RE","RW","RA")
dataPtIdx <- 1:no_pts

nameDF <- expand.grid(featureTypes, sensorLocs,dataPtIdx) %>% arrange(Var2, Var1, Var3)
nameVector <- paste0(nameDF$Var2,".",nameDF$Var1,".",nameDF$Var3)
compiledTrainingMatrixGCSe <- as.data.frame(matrix(nrow = 0,ncol = length(nameVector)))

# Prepare training set
for (j in 1:nrow(GCSe_trainLabels)){
  print(paste("Case no.",j,"started"))
  currTimeStamp <- GCSe_trainLabels$TakenInstant[j]
  currFeatSet <- firstImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCSe_trainLabels$AccelPatientNo_[j], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
  tempVector <- as.matrix(currFeatSet[,c(2,3,1,5,6,4)]) %>% as.vector()
  compiledTrainingMatrixGCSe <- rbind(compiledTrainingMatrixGCSe, tempVector)
  if (is_empty(tempVector)){print(paste("Case no.",j,"empty"))}
}

names(compiledTrainingMatrixGCSe) <- nameVector

compiledTestingMatrixGCSe <- as.data.frame(matrix(nrow = 0,ncol = length(nameVector)))
for (j in 1:nrow(GCSe_testLabels)){
  print(paste("Case no.",j,"started"))
  currTimeStamp <- GCSe_testLabels$TakenInstant[j]
  currFeatSet <- firstImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCSe_testLabels$AccelPatientNo_[j], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
  tempVector <- as.matrix(currFeatSet[,c(2,3,1,5,6,4)]) %>% as.vector()
  compiledTestingMatrixGCSe <- rbind(compiledTestingMatrixGCSe, tempVector)
  if (is_empty(tempVector)){print(paste("Case no.",j,"empty"))}
}

names(compiledTestingMatrixGCSe) <- nameVector

set.seed(9560)
down_train <- downSample(x = compiledTrainingMatrixGCSe,
                         y = as.factor(GCSe_trainLabels$Eye.Opening))

trainLOL <- lol.project.lol(as.matrix(compiledTrainingMatrixGCSe),GCSe_trainLabels$Eye.Opening,r = 200)
testLOL <- lol.embed(as.matrix(compiledTestingMatrixGCSe),trainLOL$A)

LDA_class_train <- MASS::lda(as.matrix(compiledTrainingMatrixGCSe), GCSe_trainLabels$Eye.Opening)
testPredicts <- predict(LDA_class_train,testLOL)

confusionMatrix(testPredicts$class,as.factor(GCSe_testLabels$Eye.Opening))

lhat <- 1 - sum(testPredicts$class == GCSe_trainLabels$Eye.Opening)/length(GCSe_trainLabels$Eye.Opening)

