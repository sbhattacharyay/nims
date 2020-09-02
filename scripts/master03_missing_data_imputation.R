#### Master Script 3: Multiple Imputation of Missing Accelerometery Values ####
# Decoding Quantitative Motor Features for Classification and Prediction
# in Severe Acquired Brain Injury
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

.libPaths(c("~/Rpackages" , .libPaths()))

library(tidyverse)
library(forecast)
library(Amelia)
library(imputeTS)
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

# Load Motion Features (all)
if (!exists("all_motion_features")) {
  all_sensors <-
    readMat('../all_motion_feature_data/complete_sensor_data.mat')$sensors
  all_motion_features <- do.call(rbind, all_sensors)
  featureLabels <- read.csv('../all_motion_feature_data/feature_names.csv',header = FALSE)
  all_mf_times <- read.csv('../all_motion_feature_data/indexed_times.csv') %>% mutate(times = as.POSIXct(times,tz="America/New_York",format ="%d-%b-%Y %H:%M:%S"))
}

n <- length(all_motion_features)/length(featureLabels)

source('./functions/mf_to_dataframe.R')
# Recode all missing values to NA, find "totally missing" data-streams, and store all feature values in single DF:

out.MF2DF <- mf_to_dataframe(all_motion_features,n,verbose = TRUE) 
completeFeatureSet <- out.MF2DF[[1]]
totallyMissingSet <- out.MF2DF[[2]]

featureLabels <- unlist(featureLabels[1,])

complete_timestamps <- as.data.frame(matrix(ncol = 3, nrow = 0))
for (i in 1:n){
  curr_timestamps <- all_mf_times %>% filter(ptIdx == i)
  curr_timestamps$timeCount <- seq(nrow(curr_timestamps))
  complete_timestamps <- rbind(complete_timestamps,curr_timestamps)
}
complete_timestamps <- rename(complete_timestamps,timeStamps = times)
completeFeatureSet <- inner_join(completeFeatureSet,complete_timestamps,by=c("ptIdx","timeCount"))

rm(all_motion_features,all_sensors,out.MF2DF,complete_timestamps,all_mf_times)
gc()
# From the totallyMissingSet dataframe, we see that the only incident in which a patient had more than one sensor completely missing was patient 6, who
# had sensors 2 (LA) and 6 (RE) missing.

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
  currRows <- which(completeFeatureSet$ptIdx == totallyMissingUE$ptIdx[i] & completeFeatureSet$featureType == featureLabels[[totallyMissingUE$ftIdx[i]]])
  currDF <- completeFeatureSet %>% filter(ptIdx == totallyMissingUE$ptIdx[i] & featureType == featureLabels[[totallyMissingUE$ftIdx[i]]])
  uniqIdx <- which(uniqUECombos$srIdx == totallyMissingUE$srIdx[i] & uniqUECombos$ftIdx == totallyMissingUE$ftIdx[i])
  curr_UE_mdl <- UEmdls[[uniqIdx]]
  curr_E_bxcx <- UEbxcx[[uniqIdx]][[1]]
  curr_W_bxcx <- UEbxcx[[uniqIdx]][[2]]
  if (totallyMissingUE$srIdx[i] == 3){
    currDF$LW <- predict(curr_W_bxcx,newdata = currDF$LW)
    tfOuts <- predict(curr_UE_mdl,newdata = currDF)
    Outs <- predict(curr_E_bxcx,newdata = tfOuts,inverse = TRUE)
  } else if (totallyMissingUE$srIdx[i] == 4) {
    currDF$LE <- predict(curr_E_bxcx,newdata = currDF$LE)
    tfOuts <- predict(curr_UE_mdl,newdata = currDF)
    Outs <- predict(curr_W_bxcx,newdata = tfOuts,inverse = TRUE)    
  } else if (totallyMissingUE$srIdx[i] == 6) {
    currDF$RW <- predict(curr_W_bxcx,newdata = currDF$RW)
    tfOuts <- predict(curr_UE_mdl,newdata = currDF)
    Outs <- predict(curr_E_bxcx,newdata = tfOuts,inverse = TRUE)
  } else if (totallyMissingUE$srIdx[i] == 7) {
    currDF$RE <- predict(curr_E_bxcx,newdata = currDF$RE)
    tfOuts <- predict(curr_UE_mdl,newdata = currDF)
    Outs <- predict(curr_W_bxcx,newdata = tfOuts,inverse = TRUE)   
  }
  completeFeatureSet[currRows,totallyMissingUE$srIdx[i]] <- Outs
}
rm(UEmdls,UEmdl,UEbxcx,curr_E_bxcx,curr_W_bxcx,curr_UE_mdl,LE_bxcx,LW_bxcx,RE_bxcx,RW_bxcx)

# Impute totally missing bed values by randomly sampling (with replacement) from available bed data of the same sensor:
totallyMissingBed <- totallyMissingSet %>% filter(srIdx == 1)
for (i in 1:nrow(totallyMissingBed)){
  currRows <- which(completeFeatureSet$ptIdx == totallyMissingBed$ptIdx[i] & completeFeatureSet$featureType == featureLabels[[totallyMissingBed$ftIdx[i]]])
  filtered_bedSet <- completeFeatureSet %>% drop_na(Bed) %>% filter(featureType == featureLabels[[totallyMissingBed$ftIdx[i]]])
  completeFeatureSet[currRows,1] <- sample(x =filtered_bedSet$Bed,size = length(currRows),replace = TRUE)
}
rm(filtered_bedSet)
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
  currRows <- which(completeFeatureSet$ptIdx == totallyMissingLE$ptIdx[i] & completeFeatureSet$featureType == featureLabels[[totallyMissingLE$ftIdx[i]]])
  currDF <- completeFeatureSet %>% filter(ptIdx == totallyMissingLE$ptIdx[i] & featureType == featureLabels[[totallyMissingLE$ftIdx[i]]])
  uniqIdx <- which(uniqLECombos$srIdx == totallyMissingLE$srIdx[i] & uniqLECombos$ftIdx == totallyMissingLE$ftIdx[i])
  curr_mdl <- LEmdls[[uniqIdx]]
  curr_LA_bxcx <- LEbxcx[[uniqIdx]][[1]]
  curr_RA_bxcx <- LEbxcx[[uniqIdx]][[2]]
  if (totallyMissingLE$srIdx[i] == 2){
    currDF$RA <- predict(curr_RA_bxcx,newdata = currDF$RA)
    tfOuts <- predict(curr_mdl,newdata = currDF)
    Outs <- predict(curr_LA_bxcx,newdata = tfOuts,inverse = TRUE) 
  } else if (totallyMissingLE$srIdx[i] == 5){
    currDF$LA <- predict(curr_LA_bxcx,newdata = currDF$LA)
    tfOuts <- predict(curr_mdl,newdata = currDF)
    Outs <- predict(curr_RA_bxcx,newdata = tfOuts,inverse = TRUE) 
  }
  completeFeatureSet[currRows,totallyMissingLE$srIdx[i]] <- Outs
}
rm(LEmdls,LEmdl,LEbxcx,curr_LA_bxcx,curr_RA_bxcx,curr_mdl,LA_bxcx,RA_bxcx,filteredSet)

stored_amelias <- vector(mode = "list")
stored_bxcx <- vector(mode = "list")

for (i in 1:length(featureLabels)){
  curr_amelia_DF <- data.frame(matrix(ncol = ncol(completeFeatureSet), nrow = 0))
  names(curr_amelia_DF) <- names(completeFeatureSet)
  amelia_bxcx <- vector(mode = "list")
  for (j in 1:n){
    currDF <- completeFeatureSet %>% filter(ptIdx == j & featureType == featureLabels[[i]])
    temp_bxcx <- vector(mode = "list")
    # transform each stream with univariate boxcox
    for (k in 1:7) {
      curr_bxcx <- boxcox(currDF[,k],standardize = TRUE)
      currDF[,k] <- curr_bxcx$x.t
      temp_bxcx[[k]] <- curr_bxcx
    }
    amelia_bxcx[[j]] <- temp_bxcx
    curr_amelia_DF <- rbind(curr_amelia_DF,currDF)
  }
  curr_amelia_DF <- curr_amelia_DF %>% dplyr::select(-featureType)
  curr_amelia <- amelia(curr_amelia_DF, m = 9, ts = "timeStamps", cs ="ptIdx",polytime=2,intercs = FALSE, p2s = 2)
  stored_amelias[[i]] <- curr_amelia
  stored_bxcx[[i]] <- amelia_bxcx
  print(paste('Feature no.',i,'complete'))
}

dir.create('../all_motion_feature_data/imputed_features',showWarnings = FALSE)

for (i in 1:length(stored_amelias)){
  curr_amelia <- stored_amelias[[i]]
  curr_bxcx <- stored_bxcx[[i]]
  for(l in 1:curr_amelia$m){
    curr_imp <- curr_amelia$imputations[[l]]
    for (j in 1:n){
      curr_imp_pt <- curr_imp %>% filter(ptIdx == j)
      rows_for_change <- which(curr_imp$ptIdx == j)
      temp_bxcx <- curr_bxcx[[j]]
      for (k in 1:length(curr_bxcx[[j]])){
        curr_vec <- predict(temp_bxcx[[k]],newdata = curr_imp_pt[,k],inverse = TRUE)
        if (sum(is.na(curr_vec)) > 0){
          curr_vec <- na_interpolation(curr_vec, option = "linear")
        }
        curr_imp_pt[,k] <- curr_vec
      }
      curr_imp[rows_for_change,] <- curr_imp_pt
    }
      fileName <- paste0(featureLabels[[i]],"_",l,".csv")
      write.csv(curr_imp,file.path("../all_motion_feature_data/imputed_features",fileName))
  }
  print(paste('Feature no.',i,'complete'))
}