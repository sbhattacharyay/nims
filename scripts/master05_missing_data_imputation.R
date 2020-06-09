#### Master Script 5: Multiple Imputation of Missing Accelerometery Values ####
# Decoding Quantitative Motor Features for Classification and Prediction
# in Severe Acquired Brain Injury
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

library(tidyverse)
library(forecast)
library(Amelia)
library(parallel)
library(R.matlab)
library(car)
library(bestNormalize)
library(tseries)
library(seastests)
if(.Platform$OS.type == "unix") {
  library(doMC)
} else {
  library(doParallel)
}

# Set the number of parallel cores for parallel tuning
no.parallel.cores <- floor(2 * detectCores() / 3)
if(.Platform$OS.type == "unix") {
  registerDoMC(cores = no.parallel.cores)
} else {
  registerDoParallel(cores = no.parallel.cores)
}

source('./functions/load_patient_clinical_data.R')

patient_clinical_data = load_patient_clinical_data('../clinical_data/patient_clinical_data.csv')

source('./functions/get_motion_features.R')

# Load Motion Features Organized by Time from Recording (TFR)
if (!exists("tfr_motion_features")) {
  tfr_sensors <-
    readMat('../tfr_motion_feature_data/complete_sensor_data.mat')$sensors
  tfr_motion_features <- do.call(rbind, tfr_sensors)
  featureLabels <- read.csv('../tfr_motion_feature_data/feature_names.csv',header = FALSE)
}

# Set number of features:
feature_count <- 7

temp_dims <- dim(tfr_motion_features[[1]])
sensor_count <- temp_dims[1]
n <- length(tfr_motion_features)/feature_count

# Recode all missing values to NA, find "totally missing" data-streams:
mfIdx <- c()
srIdx <- c()
corr_motion_features <- vector(mode = 'list')

for (i in 1:length(tfr_motion_features)) {
  curr_matrix <- ts(t(tfr_motion_features[[i]]),frequency = 720)
  colnames(curr_matrix) <- c("Bed","LA","LE","LW","RA","RE","RW")
  curr_matrix[is.nan(curr_matrix)] <- NA
  curr_matrix[is.infinite(curr_matrix)] <- NA
  for (j in 1:ncol(curr_matrix)) {
    if (all(curr_matrix[,j] == 0,na.rm = TRUE)){
      curr_matrix[,j] <- rep(NA, length(curr_matrix[,j]))
    }
    if (all(is.na(curr_matrix[,j]))) {
      mfIdx <- c(mfIdx,i)
      srIdx <- c(srIdx,j)
    }
  }
  corr_motion_features[[i]] <- curr_matrix
}

rm(tfr_motion_features, tfr_sensors)

# Convert mfIdx to ptIdx and ftIdx
ptIdx <- c()
ftIdx <- c()

for (i in 1:length(mfIdx)){
  if(mfIdx[i] %% n == 0){
    ptIdx <- c(ptIdx,n)
  } else {
    ptIdx <- c(ptIdx,mfIdx[i] %% n)
  }
  ftIdx <- c(ftIdx,ceiling(mfIdx[i]/69))
}
totallyMissingSet <- as.data.frame(cbind(ptIdx,ftIdx,srIdx,mfIdx))

# From this dataframe, we see that the only incident in which a patient had more than one sensor completely missing was patient 6, who
# had sensors 2 (LA) and 6 (RE) missing.

# Here we compile a complete set of all the recorded motion features across the set
completeFeatureSet <- data.frame(matrix(ncol = 9, nrow = 0))
names(completeFeatureSet) <- c("Bed","LA","LE","LW","RA","RE","RW","featureType","pNum")

for (i in 1:length(corr_motion_features)){
  tempDF <- as.data.frame(corr_motion_features[[i]]) %>% filter_all(any_vars(!is.na(.)))
  if(i %% n == 0){
    current_pNum <- n
  } else {
    current_pNum <- i %% n
  }
  current_feature <- featureLabels[,ceiling(i/69)]
  tempDF$featureType <- current_feature
  tempDF$pNum <- current_pNum
  completeFeatureSet <- rbind(completeFeatureSet,tempDF)
  print(paste('set no',i,'complete'))
}

# First separate out missing UE indices:
totallyMissingUE <- totallyMissingSet %>% filter(srIdx%in%c(3,4,6,7))
uniqUECombos <- unique(totallyMissingUE[,c('srIdx','ftIdx')])
UEmdls <- vector(mode = "list")
UEbxcx <- vector(mode = "list")

for (i in 1:nrow(uniqUECombos)){
  if(uniqUECombos$srIdx[i] %in% c(3,4)) {
    filteredSet <- completeFeatureSet %>% drop_na(LE,LW) %>% filter(featureType == featureLabels[[uniqUECombos$ftIdx[i]]])
    # Apply Box Cox
    LE_bxcx <- boxcox(filteredSet$LE,standardize = TRUE)
    LW_bxcx <- boxcox(filteredSet$LW,standardize = TRUE)
    UEbxcx[[i]]<- list(LE_bxcx,LW_bxcx)
    filteredSet$LE <- LE_bxcx$x.t
    filteredSet$LW <- LW_bxcx$x.t
    if (uniqUECombos$srIdx[i] == 3) {
      UEmdl <- lm(LE ~ LW, filteredSet)
    } else {
      UEmdl <- lm(LW ~ LE, filteredSet)
    }
  } else if (uniqUECombos$srIdx[i] %in% c(6,7)) {
    filteredSet <- completeFeatureSet %>% drop_na(RE,RW) %>% filter(featureType == featureLabels[[uniqUECombos$ftIdx[i]]])
    # Apply Box Cox
    RE_bxcx <- boxcox(filteredSet$RE,standardize = TRUE)
    RW_bxcx <- boxcox(filteredSet$RW,standardize = TRUE)
    UEbxcx[[i]]<- list(RE_bxcx,RW_bxcx)
    filteredSet$RE <- RE_bxcx$x.t
    filteredSet$RW <- RW_bxcx$x.t
    if (uniqUECombos$srIdx[i] == 6) {
      UEmdl <- lm(RE ~ RW, filteredSet)
    } else {
      UEmdl <- lm(RW ~ RE, filteredSet)
    }
  }
  UEmdls[[i]] <- UEmdl
  print(paste('combination no',i,'complete'))
}

# Impute all totally missing upper extremity values possible with regression models:
for (i in 1:nrow(totallyMissingUE)){
  curr_matrix <- corr_motion_features[[totallyMissingUE$mfIdx[i]]]
  uniqIdx <- which(uniqUECombos$srIdx == totallyMissingUE$srIdx[i] & uniqUECombos$ftIdx == totallyMissingUE$ftIdx[i])
  curr_UE_mdl <- UEmdls[[uniqIdx]]
  curr_E_bxcx <- UEbxcx[[uniqIdx]][[1]]
  curr_W_bxcx <- UEbxcx[[uniqIdx]][[2]]
  imputed_sensor<- curr_matrix[,totallyMissingUE$srIdx[i]] 
  if (totallyMissingUE$srIdx[i] == 3){
    tfPreds <- predict(curr_W_bxcx,newdata = curr_matrix[,4]) 
    predictor_DF <- data.frame(LE=imputed_sensor,LW=tfPreds)
    tfOuts<-predict(curr_UE_mdl,newdata = predictor_DF)
    Outs <- predict(curr_E_bxcx,newdata = tfOuts,inverse = TRUE)
  } else if (totallyMissingUE$srIdx[i] == 4) {
    tfPreds <- predict(curr_E_bxcx,newdata = curr_matrix[,3]) 
    predictor_DF <- data.frame(LE=tfPreds,LW=imputed_sensor)
    tfOuts<-predict(curr_UE_mdl,newdata = predictor_DF)
    Outs <- predict(curr_W_bxcx,newdata = tfOuts,inverse = TRUE)
  } else if (totallyMissingUE$srIdx[i] == 6) {
    tfPreds <- predict(curr_W_bxcx,newdata = curr_matrix[,7]) 
    predictor_DF <- data.frame(RE=imputed_sensor,RW=tfPreds)
    tfOuts<-predict(curr_UE_mdl,newdata = predictor_DF)
    Outs <- predict(curr_E_bxcx,newdata = tfOuts,inverse = TRUE)
  } else if (totallyMissingUE$srIdx[i] == 7) {
    tfPreds <- predict(curr_E_bxcx,newdata = curr_matrix[,6]) 
    predictor_DF <- data.frame(RE=tfPreds,LW=imputed_sensor)
    tfOuts<-predict(curr_UE_mdl,newdata = predictor_DF)
    Outs <- predict(curr_W_bxcx,new,newdata = tfOuts,inverse = TRUE)
  }
  curr_matrix[,totallyMissingUE$srIdx[i]] <- Outs
  corr_motion_features[[totallyMissingUE$mfIdx[i]]] <- curr_matrix
}
rm(UEmdls)

# Impute totally missing bed values by randomly sampling (with replacement) from available bed data of the same sensor:
totallyMissingBed <- totallyMissingSet %>% filter(srIdx == 1)
for (i in 1:nrow(totallyMissingBed)){
  curr_matrix <- corr_motion_features[[totallyMissingBed$mfIdx[i]]]
  filtered_bedSet <- completeFeatureSet %>% drop_na(Bed) %>% filter(featureType == featureLabels[[totallyMissingBed$ftIdx[i]]])
  curr_matrix[,1]<-sample(x = filtered_bedSet$Bed,size = nrow(curr_matrix),replace = TRUE)
  corr_motion_features[[totallyMissingBed$mfIdx[i]]] <- curr_matrix
}

# Focus on exploratory analysis of prediction of ankle regression:
# smaFilter <- completeFeatureSet %>% filter(featureType == "sma")
# LA_bxcx <- boxcox(smaFilter$LA,standardize = TRUE)
# LE_bxcx <- boxcox(smaFilter$LE,standardize = TRUE)
# LW_bxcx <- boxcox(smaFilter$LW,standardize = TRUE)
# RA_bxcx <- boxcox(smaFilter$RA,standardize = TRUE)
# RE_bxcx <- boxcox(smaFilter$RE,standardize = TRUE)
# RW_bxcx <- boxcox(smaFilter$RW,standardize = TRUE)
# 
# smaFilter$LA <- predict(LA_bxcx,newdata = smaFilter$LA) 
# smaFilter$LE <- predict(LE_bxcx,newdata = smaFilter$LE) 
# smaFilter$LW <- predict(LW_bxcx,newdata = smaFilter$LW) 
# smaFilter$RA <- predict(RA_bxcx,newdata = smaFilter$RA) 
# smaFilter$RE <- predict(RE_bxcx,newdata = smaFilter$RE) 
# smaFilter$RW <- predict(RW_bxcx,newdata = smaFilter$RW) 
# 
# tempMdl1 <- lm(LA ~ RA,smaFilter)
# tempMdl2 <- lm(LA ~ LE + LW + LE:LW,smaFilter)

# Based on this analysis, I will use the other available Ankle to impute missing ankle values
# First separate out missing Lower Extremity indices:
totallyMissingLE <- totallyMissingSet %>% filter(srIdx%in%c(2,5))
uniqLECombos <- unique(totallyMissingLE[,c('srIdx','ftIdx')])
LEmdls <- vector(mode = "list")
LEbxcx <- vector(mode = "list")

for (i in 1:nrow(uniqLECombos)){
  filteredSet <- completeFeatureSet %>% drop_na(LA,RA) %>% filter(featureType == featureLabels[[uniqLECombos$ftIdx[i]]]) 
  LA_bxcx <- boxcox(filteredSet$LA,standardize = TRUE)
  RA_bxcx <- boxcox(filteredSet$RA,standardize = TRUE)
  LEbxcx[[i]]<- list(LA_bxcx,RA_bxcx)
  filteredSet$LA <- LA_bxcx$x.t
  filteredSet$RA <- RA_bxcx$x.t
  if (uniqLECombos$srIdx[i] == 2){
    LEmdl <- lm(LA ~ RA, filteredSet)
  } else {
    LEmdl <- lm(RA ~ LA, filteredSet)
  }
  LEmdls[[i]] <- LEmdl
  print(paste('combination no',i,'complete'))
}

for (i in 1:nrow(totallyMissingLE)){
  curr_matrix <- corr_motion_features[[totallyMissingLE$mfIdx[i]]]
  uniqIdx <- which(uniqLECombos$srIdx == totallyMissingLE$srIdx[i] & uniqLECombos$ftIdx == totallyMissingLE$ftIdx[i])
  curr_mdl <- LEmdls[[uniqIdx]]
  curr_LA_bxcx <- LEbxcx[[uniqIdx]][[1]]
  curr_RA_bxcx <- LEbxcx[[uniqIdx]][[2]]
  imputed_sensor<- curr_matrix[,totallyMissingLE$srIdx[i]] 
  if (totallyMissingLE$srIdx[i] == 2){
    tfPreds <- predict(curr_RA_bxcx,newdata = curr_matrix[,5]) 
    predictor_DF <- data.frame(LA=imputed_sensor,RA=tfPreds)
    tfOuts<-predict(curr_mdl,newdata = predictor_DF)
    Outs <- predict(curr_LA_bxcx,newdata = tfOuts,inverse = TRUE)
  } else if (totallyMissingLE$srIdx[i] == 5){
    tfPreds <- predict(curr_LA_bxcx,newdata = curr_matrix[,2]) 
    predictor_DF <- data.frame(LA=tfPreds,RA=imputed_sensor)
    tfOuts<-predict(curr_mdl,newdata = predictor_DF)
    Outs <- predict(curr_RA_bxcx,newdata = tfOuts,inverse = TRUE)
  }
  curr_matrix[,totallyMissingLE$srIdx[i]] <- Outs
  corr_motion_features[[totallyMissingLE$mfIdx[i]]] <- curr_matrix
}
rm(LEmdls)

newCompleteFeatureSet <- data.frame(matrix(ncol = 10, nrow = 0))
names(newCompleteFeatureSet) <- c("Bed","LA","LE","LW","RA","RE","RW","featureType","pNum","timeCount")
new_bxcx <- vector(mode = "list")

for (i in 1:length(corr_motion_features)){
  tempDF <- as.data.frame(corr_motion_features[[i]])
  temp_bxcx <- vector(mode = "list")
  # transform each stream with univariate boxcox
  for (j in 1:ncol(tempDF)) {
    curr_bxcx <- boxcox(tempDF[,j],standardize = TRUE)
    tempDF[,j] <- curr_bxcx$x.t
    temp_bxcx[[j]] <- curr_bxcx
  }
  new_bxcx[[i]] <- temp_bxcx
  
  if(i %% n == 0){
    current_pNum <- n
  } else {
    current_pNum <- i %% n
  }
  current_feature <- featureLabels[,ceiling(i/69)]
  tempDF$featureType <- current_feature
  tempDF$pNum <- current_pNum
  tempDF$timeCount <- 1:nrow(corr_motion_features[[i]])
  newCompleteFeatureSet <- rbind(newCompleteFeatureSet,tempDF)
  print(paste('set no',i,'complete'))
}

newCompleteBPset <- newCompleteFeatureSet %>% filter(featureType=="band_power")%>%select(-featureType)
newCompleteFEset <- newCompleteFeatureSet %>% filter(featureType=="freq_entropy")%>%select(-featureType)
newCompleteFP1set <-newCompleteFeatureSet %>% filter(featureType=="freq_pairs1")%>%select(-featureType)
newCompleteFP2set <-newCompleteFeatureSet %>% filter(featureType=="freq_pairs2")%>%select(-featureType)
newCompleteMFset <- newCompleteFeatureSet %>% filter(featureType=="med_freq")%>%select(-featureType)
newCompleteSMset <- newCompleteFeatureSet %>% filter(featureType=="sma")%>%select(-featureType)
newCompleteWVset <- newCompleteFeatureSet %>% filter(featureType=="wavelets")%>%select(-featureType)

bp_amelia <- amelia(newCompleteBPset, m = 9, ts = "timeCount", cs ="pNum",polytime=2)
fe_amelia <- amelia(newCompleteFEset, m = 9, ts = "timeCount", cs ="pNum",polytime=2)
fp1_amelia<- amelia(newCompleteFP1set,m = 9, ts = "timeCount", cs ="pNum",polytime=2)
fp1_amelia<- amelia(newCompleteFP2set,m = 9, ts = "timeCount", cs ="pNum",polytime=2)
mf_amelia <- amelia(newCompleteMFset, m = 9, ts = "timeCount", cs ="pNum",polytime=2)
sm_amelia <- amelia(newCompleteSMset, m = 9, ts = "timeCount", cs ="pNum",polytime=2)
wv_amelia <- amelia(newCompleteWVset, m = 9, ts = "timeCount", cs ="pNum",polytime=2)

# Inverse-transform and save imputations in the same format as before:
# for (i in 1:length(corr_motion_features)){
#   tempDF <- as.data.frame(corr_motion_features[[i]])
#   temp_bxcx <- vector(mode = "list")
#   # transform each stream with univariate boxcox
#   for (j in 1:ncol(tempDF)) {
#     curr_bxcx <- boxcox(tempDF[,j],standardize = TRUE)
#     tempDF[,j] <- curr_bxcx$x.t
#     temp_bxcx[[j]] <- curr_bxcx
#   }
#   new_bxcx[[i]] <- temp_bxcx
#   
#   if(i %% n == 0){
#     current_pNum <- n
#   } else {
#     current_pNum <- i %% n
#   }
#   current_feature <- featureLabels[,ceiling(i/69)]
#   tempDF$featureType <- current_feature
#   tempDF$pNum <- current_pNum
#   tempDF$timeCount <- 1:nrow(corr_motion_features[[i]])
#   newCompleteFeatureSet <- rbind(newCompleteFeatureSet,tempDF)
#   print(paste('set no',i,'complete'))
# }
