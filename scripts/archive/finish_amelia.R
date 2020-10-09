.libPaths(c("~/Rpackages" , .libPaths()))
setwd("~/work/nims/scripts")

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
patient_clinical_data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv')

accel_data_files <- sort(list.files('~/scratch/pure_accel_data/',pattern = "*.csv"))
op <- options(digits.secs=3)
n <- length(accel_data_files)

load("~/scratch/pure_accel_data/checkpoint2.RData")

# Amelia II for multiple time-series normal missing data imputation
stored_amelias <- vector(mode = "list")
stored_bxcx <- vector(mode = "list")
curr_amelia_DF <- data.frame(matrix(ncol = ncol(compiledData), nrow = 0))
names(curr_amelia_DF) <- names(compiledData)
for (j in 1:n){
  currDF <- compiledData %>% filter(ptIdx == j)
  temp_bxcx <- vector(mode = "list")
  # transform each stream with univariate boxcox
  for (k in 1:18) {
    curr_bxcx <- boxcox(currDF[,k],standardize = TRUE)
    currDF[,k] <- curr_bxcx$x.t
    temp_bxcx[[k]] <- curr_bxcx
  }
  stored_bxcx[[j]] <- temp_bxcx
  curr_amelia_DF <- rbind(curr_amelia_DF,currDF)
  print(paste('Patient no.',j,'complete'))
}
curr_amelia_DF <- curr_amelia_DF[,1:20]
stored_amelias <- amelia(curr_amelia_DF, m = 5, ts = "timeStamps", cs ="ptIdx",polytime=2,intercs = FALSE, p2s = 2)
dir.create('~/scratch/pure_accel_data/imputed_features',showWarnings = FALSE)

# Invert boxcox transformations and save Amelia II imputations
for(l in 1:stored_amelias$m){
  curr_imp <- stored_amelias$imputations[[l]]
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
    print(paste('Patient no.',j,'complete'))
  }
  fileName <- paste0("imp_no","_",l,".csv")
  write.csv(curr_imp,file.path("~/scratch/pure_accel_data/imputed_features",fileName))
  print(paste('Imputation no.',l,'complete'))
}