#### Master Script 6: Multiple Imputation of Missing Accelerometery Values ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

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

# compiledData <- as.data.frame(matrix(ncol = 23, nrow = 0))

op <- options(digits.secs=3)

# for (i in 1:length(accel_data_files)){
#   curr_pt <- read_csv(file.path('~/scratch/pure_accel_data',accel_data_files[i]))
#   substring(curr_pt$timeStamps,21,21) <- "."
#   curr_pt$timeStamps <- strptime(curr_pt$timeStamps,format = "%d-%b-%Y %H:%M:%OS",tz = "America/New_York")
#   compiledData <- rbind(compiledData,curr_pt)
#   print(paste("Patient Idx No.",i,"Complete"))
# }
n <- length(accel_data_files)

# totallyMissingSet <- data.frame(matrix(ncol = 2, nrow = 0))
# names(totallyMissingSet) <- c("ptIdx","axIdx")
# for (ptIdx in 1:length(accel_data_files)){
#   currPtIdx <- compiledData$ptIdx == ptIdx
#   currPt <- compiledData[currPtIdx,]
#   for (axIdx in 1:ncol(currPt)) {
#     if (all(currPt[,axIdx] == 0,na.rm = TRUE)){
#       currPt[,axIdx] <- rep(NA, nrow(currPt[,axIdx]))
#     }
#     if (all(is.na(currPt[,axIdx]))) {
#       totallyMissingSet<-rbind(totallyMissingSet,data.frame(ptIdx,axIdx))
#     } 
#   }
#   compiledData[currPtIdx,] <- currPt
#   print(paste("Patient Idx No.",ptIdx,"Complete"))
# }
# save(compiledData, totallyMissingSet, file = "~/scratch/pure_accel_data/compiledData.RData")

## Checkpoint 1

load("~/scratch/pure_accel_data/compiledData.RData")

# Removal of Bed Data information
compiledData <- compiledData[,4:23]

# Take absolute value of data
compiledData[,1:18] <- abs(compiledData[,1:18])

# Correct sensor indices to account for removal of bed sensors
totallyMissingSet$axIdx <- totallyMissingSet$axIdx - 3
totallyMissingSet <- totallyMissingSet[totallyMissingSet$axIdx>=1,]

# Filter out upper extremity indices
totallyMissingUE <- totallyMissingSet %>% filter(axIdx %in% c(4:9,13:18))
uniqUECombos <- unique(totallyMissingUE[,'axIdx'])

 # Intialize empty vectors to store upper extremity prediction models and boxcox objects
UEmdls <- vector(mode = "list")
UEbxcx <- vector(mode = "list")

for (i in 1:length(uniqUECombos)){
  print(paste('combination no',i,'started'))
  if(uniqUECombos[i] %in% c(4:9)) {
    filteredSet <- compiledData %>% drop_na(LE.x,LE.y,LE.z,LW.x,LW.y,LW.z)
    
    # Apply Box Cox
    LE.x_bxcx <- boxcox(filteredSet$LE.x,standardize = TRUE)
    LE.y_bxcx <- boxcox(filteredSet$LE.y,standardize = TRUE)
    LE.z_bxcx <- boxcox(filteredSet$LE.z,standardize = TRUE)
    
    LW.x_bxcx <- boxcox(filteredSet$LW.x,standardize = TRUE)
    LW.y_bxcx <- boxcox(filteredSet$LW.y,standardize = TRUE)
    LW.z_bxcx <- boxcox(filteredSet$LW.z,standardize = TRUE)    
    
    UEbxcx[[i]]<- list(LE.x_bxcx,LE.y_bxcx,LE.z_bxcx,LW.x_bxcx,LW.y_bxcx,LW.z_bxcx)
    
    filteredSet$LE.x <- LE.x_bxcx$x.t
    filteredSet$LE.y <- LE.y_bxcx$x.t
    filteredSet$LE.z <- LE.z_bxcx$x.t
    
    filteredSet$LW.x <- LW.x_bxcx$x.t
    filteredSet$LW.y <- LW.y_bxcx$x.t
    filteredSet$LW.z <- LW.z_bxcx$x.t

    if (uniqUECombos[i] == 4) {
      UEmdl <- lm(LE.x ~ LW.x + LW.y + LW.z, filteredSet)
    } else if (uniqUECombos[i] == 5) { 
      UEmdl <- lm(LE.y ~ LW.x + LW.y + LW.z, filteredSet)
    } else if (uniqUECombos[i] == 6) {
      UEmdl <- lm(LE.z ~ LW.x + LW.y + LW.z, filteredSet)
    } else if (uniqUECombos[i] == 7) {
      UEmdl <- lm(LW.x ~ LE.x + LE.y + LE.z, filteredSet)
    } else if (uniqUECombos[i] == 8) {
      UEmdl <- lm(LW.y ~ LE.x + LE.y + LE.z, filteredSet)
    } else if (uniqUECombos[i] == 9) {
      UEmdl <- lm(LW.z ~ LE.x + LE.y + LE.z, filteredSet)
    }
    
  } else if (uniqUECombos[i] %in% c(13:18)) {
    filteredSet <- compiledData %>% drop_na(RE.x,RE.y,RE.z,RW.x,RW.y,RW.z)
    
    # Apply Box Cox
    RE.x_bxcx <- boxcox(filteredSet$RE.x,standardize = TRUE)
    RE.y_bxcx <- boxcox(filteredSet$RE.y,standardize = TRUE)
    RE.z_bxcx <- boxcox(filteredSet$RE.z,standardize = TRUE)
    
    RW.x_bxcx <- boxcox(filteredSet$RW.x,standardize = TRUE)
    RW.y_bxcx <- boxcox(filteredSet$RW.y,standardize = TRUE)
    RW.z_bxcx <- boxcox(filteredSet$RW.z,standardize = TRUE)    
    
    UEbxcx[[i]]<- list(RE.x_bxcx,RE.y_bxcx,RE.z_bxcx,RW.x_bxcx,RW.y_bxcx,RW.z_bxcx)
    
    filteredSet$RE.x <- RE.x_bxcx$x.t
    filteredSet$RE.y <- RE.y_bxcx$x.t
    filteredSet$RE.z <- RE.z_bxcx$x.t
    
    filteredSet$RW.x <- RW.x_bxcx$x.t
    filteredSet$RW.y <- RW.y_bxcx$x.t
    filteredSet$RW.z <- RW.z_bxcx$x.t
    
    if (uniqUECombos[i] == 13) {
      UEmdl <- lm(RE.x ~ RW.x + RW.y + RW.z, filteredSet)
    } else if (uniqUECombos[i] == 14) { 
      UEmdl <- lm(RE.y ~ RW.x + RW.y + RW.z, filteredSet)
    } else if (uniqUECombos[i] == 15) {
      UEmdl <- lm(RE.z ~ RW.x + RW.y + RW.z, filteredSet)
    } else if (uniqUECombos[i] == 16) {
      UEmdl <- lm(RW.x ~ RE.x + RE.y + RE.z, filteredSet)
    } else if (uniqUECombos[i] == 17) {
      UEmdl <- lm(RW.y ~ RE.x + RE.y + RE.z, filteredSet)
    } else if (uniqUECombos[i] == 18) {
      UEmdl <- lm(RW.z ~ RE.x + RE.y + RE.z, filteredSet)
    }
    
  }
  UEmdls[[i]] <- UEmdl
  print(paste('combination no',i,'complete'))
}

# Impute all totally missing upper extremity values possible with regression models:
for (i in 1:nrow(totallyMissingUE)){
  currRows <- which(compiledData$ptIdx == totallyMissingUE$ptIdx[i])
  currDF <- compiledData[currRows,]
  uniqIdx <- which(uniqUECombos == totallyMissingUE$axIdx[i])
  curr_UE_mdl <- UEmdls[[uniqIdx]]
  
  curr_E.x_bxcx <- UEbxcx[[uniqIdx]][[1]]
  curr_E.y_bxcx <- UEbxcx[[uniqIdx]][[2]]
  curr_E.z_bxcx <- UEbxcx[[uniqIdx]][[3]]
  curr_W.x_bxcx <- UEbxcx[[uniqIdx]][[4]]
  curr_W.y_bxcx <- UEbxcx[[uniqIdx]][[5]]
  curr_W.z_bxcx <- UEbxcx[[uniqIdx]][[6]]
  
  if (totallyMissingUE$axIdx[i] %in% 4:6) {
    currDF$LW.x <- predict(curr_W.x_bxcx,newdata = currDF$LW.x)
    currDF$LW.y <- predict(curr_W.y_bxcx,newdata = currDF$LW.y)
    currDF$LW.z <- predict(curr_W.z_bxcx,newdata = currDF$LW.z)
  } else if (totallyMissingUE$axIdx[i] %in% 7:9) {
    currDF$LE.x <- predict(curr_E.x_bxcx,newdata = currDF$LE.x)
    currDF$LE.y <- predict(curr_E.y_bxcx,newdata = currDF$LE.y)
    currDF$LE.z <- predict(curr_E.z_bxcx,newdata = currDF$LE.z)
  } else if (totallyMissingUE$axIdx[i] %in% 13:15) {
    currDF$RW.x <- predict(curr_W.x_bxcx,newdata = currDF$RW.x)
    currDF$RW.y <- predict(curr_W.y_bxcx,newdata = currDF$RW.y)
    currDF$RW.z <- predict(curr_W.z_bxcx,newdata = currDF$RW.z)
  } else if (totallyMissingUE$axIdx[i] %in% 16:18) {
    currDF$RE.x <- predict(curr_E.x_bxcx,newdata = currDF$RE.x)
    currDF$RE.y <- predict(curr_E.y_bxcx,newdata = currDF$RE.y)
    currDF$RE.z <- predict(curr_E.z_bxcx,newdata = currDF$RE.z)
  }
  tfOuts <- predict(curr_UE_mdl,newdata = currDF)
  if (totallyMissingUE$axIdx[i] %in% c(4,13)) {
    outs <- predict(curr_E.x_bxcx,newdata = tfOuts,inverse = TRUE)
  } else if (totallyMissingUE$axIdx[i] %in% c(5,14)){
    outs <- predict(curr_E.y_bxcx,newdata = tfOuts,inverse = TRUE)
  } else if (totallyMissingUE$axIdx[i] %in% c(6,15)){
    outs <- predict(curr_E.z_bxcx,newdata = tfOuts,inverse = TRUE)
  } else if (totallyMissingUE$axIdx[i] %in% c(7,16)){
    outs <- predict(curr_W.x_bxcx,newdata = tfOuts,inverse = TRUE)
  } else if (totallyMissingUE$axIdx[i] %in% c(8,17)){
    outs <- predict(curr_W.y_bxcx,newdata = tfOuts,inverse = TRUE)
  } else if (totallyMissingUE$axIdx[i] %in% c(9,18)){
    outs <- predict(curr_W.z_bxcx,newdata = tfOuts,inverse = TRUE)
  }
  compiledData[currRows,totallyMissingUE$axIdx[i]] <- outs
}
rm(UEmdls,UEmdl,UEbxcx,curr_UE_mdl)

# Filter out missing lower extremity indices and initialize empty vectors to store lower extremity prediction models and boxcox objects
totallyMissingLE <- totallyMissingSet %>% filter(axIdx %in% c(1:3,10:12))
uniqLECombos <- unique(totallyMissingLE[,c('axIdx')])
LEmdls <- vector(mode = "list")
LEbxcx <- vector(mode = "list")

for (i in 1:length(uniqLECombos)){
  print(paste('combination no',i,'started'))
  filteredSet <- compiledData %>% drop_na(LA.x,LA.y,LA.z,RA.x,RA.y,RA.z)
  
  LA.x_bxcx <- boxcox(filteredSet$LA.x,standardize = TRUE)
  LA.y_bxcx <- boxcox(filteredSet$LA.y,standardize = TRUE)
  LA.z_bxcx <- boxcox(filteredSet$LA.z,standardize = TRUE)
  
  RA.x_bxcx <- boxcox(filteredSet$RA.x,standardize = TRUE)
  RA.y_bxcx <- boxcox(filteredSet$RA.y,standardize = TRUE)
  RA.z_bxcx <- boxcox(filteredSet$RA.z,standardize = TRUE)
  
  LEbxcx[[i]]<- list(LA.x_bxcx,LA.y_bxcx,LA.z_bxcx,RA.x_bxcx,RA.y_bxcx,RA.z_bxcx)

  filteredSet$LA.x <- LA.x_bxcx$x.t
  filteredSet$LA.y <- LA.y_bxcx$x.t
  filteredSet$LA.z <- LA.z_bxcx$x.t
  
  filteredSet$RA.x <- RA.x_bxcx$x.t
  filteredSet$RA.y <- RA.y_bxcx$x.t
  filteredSet$RA.z <- RA.z_bxcx$x.t
  
  if (uniqLECombos[i] == 1){
    LEmdl <- lm(LA.x ~ RA.x + RA.y + RA.z, filteredSet)
  } else if (uniqLECombos[i] == 2) {
    LEmdl <- lm(LA.y ~ RA.x + RA.y + RA.z, filteredSet)
  } else if (uniqLECombos[i] == 3) {
    LEmdl <- lm(LA.z ~ RA.x + RA.y + RA.z, filteredSet)
  } else if (uniqLECombos[i] == 10) {
    LEmdl <- lm(RA.x ~ LA.x + LA.y + LA.z, filteredSet)
  } else if (uniqLECombos[i] == 11) {
    LEmdl <- lm(RA.y ~ LA.x + LA.y + LA.z, filteredSet)
  } else if (uniqLECombos[i] == 12) {
    LEmdl <- lm(RA.z ~ LA.x + LA.y + LA.z, filteredSet)
  }
  LEmdls[[i]] <- LEmdl
  print(paste('combination no',i,'complete'))
}

# Impute all totally missing lower extremity values possible with regression models:
for (i in 1:nrow(totallyMissingLE)){
  currRows <- which(compiledData$ptIdx == totallyMissingLE$ptIdx[i])
  currDF <- compiledData[currRows,]
  uniqIdx <- which(uniqLECombos == totallyMissingLE$axIdx[i])
  curr_mdl <- LEmdls[[uniqIdx]]

  curr_LA.x_bxcx <- LEbxcx[[uniqIdx]][[1]]
  curr_LA.y_bxcx <- LEbxcx[[uniqIdx]][[2]]
  curr_LA.z_bxcx <- LEbxcx[[uniqIdx]][[3]]
  curr_RA.x_bxcx <- LEbxcx[[uniqIdx]][[4]]
  curr_RA.y_bxcx <- LEbxcx[[uniqIdx]][[5]]
  curr_RA.z_bxcx <- LEbxcx[[uniqIdx]][[6]]
  
  if (totallyMissingLE$axIdx[i] %in% 1:3){
    currDF$RA.x <- predict(curr_RA.x_bxcx,newdata = currDF$RA.x)
    currDF$RA.y <- predict(curr_RA.y_bxcx,newdata = currDF$RA.y)
    currDF$RA.z <- predict(curr_RA.z_bxcx,newdata = currDF$RA.z)
  } else {
    currDF$LA.x <- predict(curr_LA.x_bxcx,newdata = currDF$LA.x)
    currDF$LA.x <- predict(curr_LA.x_bxcx,newdata = currDF$LA.x)
    currDF$LA.x <- predict(curr_LA.x_bxcx,newdata = currDF$LA.x)
  }
  tfOuts <- predict(curr_mdl,newdata = currDF)
  if (totallyMissingLE$axIdx[i] == 1){
    outs <- predict(curr_LA.x_bxcx,newdata = tfOuts,inverse = TRUE) 
  } else if (totallyMissingLE$axIdx[i] == 2){
    outs <- predict(curr_LA.y_bxcx,newdata = tfOuts,inverse = TRUE) 
  } else if (totallyMissingLE$axIdx[i] == 3){
    outs <- predict(curr_LA.z_bxcx,newdata = tfOuts,inverse = TRUE) 
  } else if (totallyMissingLE$axIdx[i] == 10){
    outs <- predict(curr_RA.x_bxcx,newdata = tfOuts,inverse = TRUE) 
  } else if (totallyMissingLE$axIdx[i] == 11){
    outs <- predict(curr_RA.y_bxcx,newdata = tfOuts,inverse = TRUE) 
  } else if (totallyMissingLE$axIdx[i] == 12){
    outs <- predict(curr_RA.z_bxcx,newdata = tfOuts,inverse = TRUE) 
  }
  compiledData[currRows,totallyMissingLE$axIdx[i]] <- outs
}
rm(LEmdls,LEmdl,LEbxcx,curr_mdl)

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