# TEMPORARY WORK FILE!

.libPaths(c("~/Rpackages" , .libPaths()))

library(tidyverse)
library(R.matlab)
library(readxl)
library(caret)
library(lolR)
library(reticulate)
np <- import("numpy")
library(RcppCNPy)
source('./functions/load_patient_clinical_data.R')



patient_clinical_data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv') %>% arrange(AccelPatientNo_) 
patient_clinical_data <- patient_clinical_data %>% mutate(ptIdx = 1:nrow(patient_clinical_data))
gcs_data <- read.csv('../clinical_data/clean_auto_GCS_table.csv') %>% select(-X) %>% mutate(TakenInstant = as.POSIXct(TakenInstant, tz = "America/New_York"))
featureSet <- readRDS('~/data/all_motion_feature_data/bc_i_c_dataset.rds')
indexedTimes <- read.csv('~/data/all_motion_feature_data/indexed_times.csv') %>% mutate(times = as.POSIXct(times, format = "%d-%b-%Y %H:%M:%S",tz = "America/New_York"))
firstImp <- featureSet[[1]] %>% select(-Bed, -ImputNum)
firstImp$TakenInstant <- ""
missingTimeInfo <- read_xlsx('~/data/all_motion_feature_data/MissingPercentTable.xlsx',.name_repair = "universal") %>% mutate(Start.Timestamp = as.POSIXct(Start.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"), End.Timestamp = as.POSIXct(End.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"))
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

save(firstImp,file = '~/data/all_motion_feature_data/firstImp.RData')

## CHECKPOINT ONE: Load firstImp.RData

load('~/data/all_motion_feature_data/firstImp.RData')

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

write_csv(gcsLabels,'~/data/all_motion_feature_data/gcs_labels.csv')

# Create Data Partitions based on GCS labels, ensuring the same patient is not in both training and testing

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

## Alternative Resampling strategy: pure label split
set.seed(2020)
GCSe_new <- gcsLabels %>% drop_na(Eye.Opening)
GCSm_new <- gcsLabels %>%drop_na(Best.Motor.Response)
GCSv_new <- gcsLabels %>%drop_na(Best.Verbal.Response)
GCST_new <- gcsLabels %>%drop_na(Glasgow.Coma.Scale.Score)

GCSe_trainIdx <- createDataPartition(GCSe_new$Eye.Opening,p=.8,list = FALSE,times = 1)
GCSm_trainIdx <- createDataPartition(GCSm_new$Best.Motor.Response,p=.8,list = FALSE,times = 1)
GCSv_trainIdx <- createDataPartition(GCSv_new$Best.Verbal.Response,p=.8,list = FALSE,times = 1)
GCST_trainIdx <- createDataPartition(GCST_new$Glasgow.Coma.Scale.Score,p=.8,list = FALSE,times = 1)

GCSe_train <- GCSe_new[GCSe_trainIdx,]
GCSm_train <- GCSm_new[GCSm_trainIdx,]
GCSv_train <- GCSv_new[GCSv_trainIdx,]
GCST_train <- GCST_new[GCST_trainIdx,]

GCSe_test <- GCSe_new[-GCSe_trainIdx,]
GCSm_test <- GCSm_new[-GCSm_trainIdx,]
GCSv_test <- GCSv_new[-GCSv_trainIdx,]
GCST_test <- GCST_new[-GCST_trainIdx,]

featureTypes <- unique(firstImp$featureType)
sensorLocs <- c("LE","LW","LA","RE","RW","RA")
dataPtIdx <- 1:no_pts

nameDF <- expand.grid(featureTypes, sensorLocs,dataPtIdx) %>% arrange(Var2, Var1, Var3)
nameVector <- paste0(nameDF$Var2,".",nameDF$Var1,".",nameDF$Var3)

trainArrayGCSe <- array(dim=c(length(GCSe_trainIdx),1,721,42))
trainArrayGCSm <- array(dim=c(length(GCSm_trainIdx),1,721,42))
trainArrayGCSv <- array(dim=c(length(GCSv_trainIdx),1,721,42))
trainArrayGCST <- array(dim=c(length(GCST_trainIdx),1,721,42))

testArrayGCSe <- array(dim=c(nrow(GCSe_test),1,721,42))
testArrayGCSm <- array(dim=c(nrow(GCSm_test),1,721,42))
testArrayGCSv <- array(dim=c(nrow(GCSv_test),1,721,42))
testArrayGCST <- array(dim=c(nrow(GCST_test),1,721,42))

trainLabelsGCSe <- GCSe_train$Eye.Opening
trainLabelsGCSm <- GCSm_train$Best.Motor.Response
trainLabelsGCSv <- GCSv_train$Best.Verbal.Response
trainLabelsGCST <- GCST_train$Glasgow.Coma.Scale.Score

testLabelsGCSe <- GCSe_test$Eye.Opening
testLabelsGCSm <- GCSm_test$Best.Motor.Response
testLabelsGCSv <- GCSv_test$Best.Verbal.Response
testLabelsGCST <- GCST_test$Glasgow.Coma.Scale.Score

# Eye Training
for (i in 1:nrow(GCSe_trainIdx)){
  print(paste("Case no.",i,"started"))
  currTimeStamp <- GCSe_train$TakenInstant[i]
  currFeatSet <- firstImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCSe_train$AccelPatientNo_[i], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
  tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
  tempReshape <- as.matrix(tempReshape[,-1])
  trainArrayGCSe[i,1,,] <- tempReshape
}

# Motor Training
for (i in 1:nrow(GCSm_trainIdx)){
  print(paste("Case no.",i,"started"))
  currTimeStamp <- GCSm_train$TakenInstant[i]
  currFeatSet <- firstImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCSm_train$AccelPatientNo_[i], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
  tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
  tempReshape <- as.matrix(tempReshape[,-1])
  trainArrayGCSm[i,1,,] <- tempReshape
}

# Verbal Training
for (i in 1:nrow(GCSv_trainIdx)){
  print(paste("Case no.",i,"started"))
  currTimeStamp <- GCSv_train$TakenInstant[i]
  currFeatSet <- firstImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCSv_train$AccelPatientNo_[i], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
  tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
  tempReshape <- as.matrix(tempReshape[,-1])
  trainArrayGCSv[i,1,,] <- tempReshape
}

# Total Training
for (i in 1:nrow(GCST_trainIdx)){
  print(paste("Case no.",i,"started"))
  currTimeStamp <- GCST_train$TakenInstant[i]
  currFeatSet <- firstImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCST_train$AccelPatientNo_[i], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
  tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
  tempReshape <- as.matrix(tempReshape[,-1])
  trainArrayGCST[i,1,,] <- tempReshape
}

# Eye Testing
for (i in 1:nrow(GCSe_test)){
  print(paste("Case no.",i,"started"))
  currTimeStamp <- GCSe_test$TakenInstant[i]
  currFeatSet <- firstImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCSe_test$AccelPatientNo_[i], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
  tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
  tempReshape <- as.matrix(tempReshape[,-1])
  testArrayGCSe[i,1,,] <- tempReshape
}

# Motor Testing
for (i in 1:nrow(GCSm_test)){
  print(paste("Case no.",i,"started"))
  currTimeStamp <- GCSm_test$TakenInstant[i]
  currFeatSet <- firstImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCSm_test$AccelPatientNo_[i], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
  tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
  tempReshape <- as.matrix(tempReshape[,-1])
  testArrayGCSm[i,1,,] <- tempReshape
}

# Verbal Testing
for (i in 1:nrow(GCSv_test)){
  print(paste("Case no.",i,"started"))
  currTimeStamp <- GCSv_test$TakenInstant[i]
  currFeatSet <- firstImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCSv_test$AccelPatientNo_[i], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
  tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
  tempReshape <- as.matrix(tempReshape[,-1])
  testArrayGCSv[i,1,,] <- tempReshape
}

# Total Testing

for (i in 1:nrow(GCST_test)){
  print(paste("Case no.",i,"started"))
  currTimeStamp <- GCST_test$TakenInstant[i]
  currFeatSet <- firstImp %>% filter(patient_clinical_data$AccelPatientNo_[ptIdx] == GCST_test$AccelPatientNo_[i], TakenInstant <= currTimeStamp, TakenInstant >= currTimeStamp - secs_window)
  tempReshape <- currFeatSet %>% pivot_longer(cols=c(2,3,1,5,6,4),names_to = "sensor") %>% pivot_wider(id_cols = "timeCount",names_from=c("featureType","sensor"),values_from="value")
  tempReshape <- as.matrix(tempReshape[,-1])
  testArrayGCST[i,1,,] <- tempReshape
}

# Save all training and testing sets
save(trainArrayGCSe,trainArrayGCSm,trainArrayGCSv,trainArrayGCST,testArrayGCSe,testArrayGCSm,testArrayGCSv,testArrayGCST,trainLabelsGCSe,trainLabelsGCSm,trainLabelsGCSv,trainLabelsGCST,testLabelsGCSe,testLabelsGCSm,testLabelsGCSv,testLabelsGCST,file="~/data/all_motion_feature_data/train_and_test_sets.RData")
dir.create("~/data/all_motion_feature_data/train_test_sets")
save(GCSe_trainIdx,GCSm_trainIdx,GCSv_trainIdx,GCST_trainIdx,file="~/data/all_motion_feature_data/train_test_sets/train_idxs.RData")

# CHECKPOINT 2
load("~/data/all_motion_feature_data/train_test_sets/train_and_test_sets.RData")
load("~/data/all_motion_feature_data/train_test_sets/train_idxs.RData")

np$save('/home-4/sbhatt15@jhu.edu/data/all_motion_feature_data/train_test_sets/trainArrayGCSe',trainArrayGCSe)
np$save('/home-4/sbhatt15@jhu.edu/data/all_motion_feature_data/train_test_sets/trainArrayGCSm',trainArrayGCSm)
np$save('/home-4/sbhatt15@jhu.edu/data/all_motion_feature_data/train_test_sets/trainArrayGCSv',trainArrayGCSv)
np$save('/home-4/sbhatt15@jhu.edu/data/all_motion_feature_data/train_test_sets/trainArrayGCST',trainArrayGCST)

np$save('/home-4/sbhatt15@jhu.edu/data/all_motion_feature_data/train_test_sets/testArrayGCSe',testArrayGCSe)
np$save('/home-4/sbhatt15@jhu.edu/data/all_motion_feature_data/train_test_sets/testArrayGCSm',testArrayGCSm)
np$save('/home-4/sbhatt15@jhu.edu/data/all_motion_feature_data/train_test_sets/testArrayGCSv',testArrayGCSv)
np$save('/home-4/sbhatt15@jhu.edu/data/all_motion_feature_data/train_test_sets/testArrayGCST',testArrayGCST)

np$save('/home-4/sbhatt15@jhu.edu/data/all_motion_feature_data/train_test_sets/trainLabelsGCSe',trainLabelsGCSe)
np$save('/home-4/sbhatt15@jhu.edu/data/all_motion_feature_data/train_test_sets/trainLabelsGCSm',trainLabelsGCSm)
np$save('/home-4/sbhatt15@jhu.edu/data/all_motion_feature_data/train_test_sets/trainLabelsGCSv',trainLabelsGCSv)
np$save('/home-4/sbhatt15@jhu.edu/data/all_motion_feature_data/train_test_sets/trainLabelsGCST',trainLabelsGCST)

np$save('/home-4/sbhatt15@jhu.edu/data/all_motion_feature_data/train_test_sets/testLabelsGCSe',testLabelsGCSe)
np$save('/home-4/sbhatt15@jhu.edu/data/all_motion_feature_data/train_test_sets/testLabelsGCSm',testLabelsGCSm)
np$save('/home-4/sbhatt15@jhu.edu/data/all_motion_feature_data/train_test_sets/testLabelsGCSv',testLabelsGCSv)
np$save('/home-4/sbhatt15@jhu.edu/data/all_motion_feature_data/train_test_sets/testLabelsGCST',testLabelsGCST)

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

