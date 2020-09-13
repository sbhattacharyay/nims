#### Master Script 6: Multiple Imputation of Missing Accelerometery Values ####
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
library(readxl)
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

accel_data_files <- sort(list.files('~/data/pure_accel_data/',pattern = "*.csv"))


compiledData <- as.data.frame(matrix(ncol = 23, nrow = 0))

op <- options(digits.secs=3)

for (i in 1:length(accel_data_files)){
  
  curr_pt <- read_csv(file.path('~/data/pure_accel_data',accel_data_files[i]))
  substring(curr_pt$timeStamps,21,21) <- "."
  curr_pt$timeStamps <- strptime(curr_pt$timeStamps,format = "%d-%b-%Y %H:%M:%OS",tz = "America/New_York")
  compiledData <- rbind(compiledData,curr_pt)
  
  print(paste("Patient Idx No.",i,"Complete"))
}
n <- length(accel_data_files)

totallyMissingSet <- data.frame(matrix(ncol = 2, nrow = 0))
names(totallyMissingSet) <- c("ptIdx","axIdx")
for (ptIdx in 1:length(accel_data_files)){
  currPtIdx <- compiledData$ptIdx == ptIdx
  currPt <- compiledData[currPtIdx,]
  for (axIdx in 1:ncol(currPt)) {
    if (all(currPt[,axIdx] == 0,na.rm = TRUE)){
      currPt[,axIdx] <- rep(NA, nrow(currPt[,axIdx]))
    }
    if (all(is.na(currPt[,axIdx]))) {
      totallyMissingSet<-rbind(totallyMissingSet,data.frame(ptIdx,axIdx))
    } 
  }
  compiledData[currPtIdx,] <- currPt
  print(paste("Patient Idx No.",ptIdx,"Complete"))
}

save(compiledData, totallyMissingSet, file = "~/data/pure_accel_data/compiledData.RData")

load("~/data/pure_accel_data/compiledData.RData")

compiledData <- compiledData[,4:23]

totallyMissingSet$axIdx <- totallyMissingSet$axIdx - 3
totallyMissingSet <- totallyMissingSet[totallyMissingSet$axIdx>=1,]

totallyMissingUE <- totallyMissingSet %>% filter(axIdx%in%c(4:9,13:18))

uniqUECombos <- unique(totallyMissingUE[,'axIdx'])

UEmdls <- vector(mode = "list")
UEbxcx <- vector(mode = "list")

for (i in 1:length(uniqUECombos)){
  if(uniqUECombos[i] %in% c(4:9)) {
    
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


gcs_table <- read_xlsx('../clinical_data/GCS_table.xlsx')








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