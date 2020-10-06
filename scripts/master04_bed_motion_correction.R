#### Master Script 4: Bed Motion Correction and Collection of Multiple Imputations ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

# Denote number of imputations
m <- 9

library(tidyverse)
library(readxl)
source('./functions/find_thresholds.R')

featureLabels <- read.csv('../all_motion_feature_data/feature_names.csv',header = FALSE)
featureLabels <- unlist(featureLabels[1,])
compiledImputations <- vector(mode = "list")

for (i in 1:m){
  print(paste('Imputation no.',i,'started'))
  currPattern <- paste0('*',i,'.csv')
  currFileList <- list.files('../all_motion_feature_data/imputed_features/',pattern = currPattern)
  imputationDF <- data.frame(matrix(ncol = 10, nrow = 0))
  for (j in 1:length(featureLabels)){
    currFilePath <- file.path('../all_motion_feature_data/imputed_features',currFileList[j]) 
    currDF <- read.csv(currFilePath) %>% dplyr::select(-X)
    currDF$featureType <- featureLabels[[j]]
    imputationDF <- rbind(imputationDF,currDF)
    print(paste('Feature no.',j,'complete'))
  }
  compiledImputations[[i]] <- imputationDF
  print(paste('Imputation no.',i,'complete'))
}

sma_thresh <- .135

# Based on the SMA threshold, find corresponding thresholds for the other feature spaces. Note we only use the first imputation to do so since, when testing, we found that the thresholds are the same for all imputations
feature_thresholds <- find_thresholds(compiledImputations[[1]],sma_thresh)
featRanges <-
  list(
    bp_nm_range = c(0, feature_thresholds[1]),
    fe_nm_range = c(feature_thresholds[2], 1.707),
    fp1_nm_range = c(0, feature_thresholds[3]),
    fp2_nm_range = c(0, feature_thresholds[4]),
    mf_nm_range = c(feature_thresholds[5], 3.2),
    sm_nm_range = c(0, feature_thresholds[6]),
    wv_nm_range = c(0, feature_thresholds[7])
  )
# Correct bed motion by identifying time points and patient indices during which SMA exceeds a literature-reviewed threshold:
for (i in 1:length(compiledImputations)){
  print(paste("Imputation No.",i,"Started"))
  bedActDF <- data.frame()
  currDF <- compiledImputations[[i]]
  currSMAdf <-currDF %>% filter(featureType == "sma")
  bedSMARows <- which(currSMAdf$Bed > sma_thresh)
  bedActRows <- which(currDF$Bed > sma_thresh & currDF$featureType == "sma")
  bedActDF <- data.frame(ptIdx=currDF$ptIdx[bedActRows],timeCount=currDF$timeCount[bedActRows])
  for (j in 1:length(featureLabels)){
    a <- featRanges[[j]][1]
    b <- featRanges[[j]][2]
    print(paste("Feature No.",j,"Started"))
    currFeatRows <- which(currDF$featureType == featureLabels[[j]])
    currFeatDF <- currDF %>% filter(featureType == featureLabels[[j]])
    currFeatChangeDF <- currFeatDF[bedSMARows,]
    if (featureLabels[[j]] %in% c("freq_entropy","med_freq")) {
      tempMat <- currFeatChangeDF[,2:7] + currFeatChangeDF[,1]
      tempMat[tempMat > b] <- runif(sum(tempMat > b),a+(b-a)/2,b)
    } else {
      tempMat <- currFeatChangeDF[,2:7] - currFeatChangeDF[,1]
      tempMat[tempMat < 0] <- runif(sum(tempMat < 0),0,b/2)
    }
    currFeatChangeDF[,2:7] <- tempMat
    currFeatDF[bedSMARows,] <- currFeatChangeDF
    currDF[currFeatRows,] <- currFeatDF
    print(paste("Feature No.",j,"Complete"))
  }
  compiledImputations[[i]] <- currDF
  print(paste("Imputation No.",i,"Complete"))
}

for (i in 1:length(compiledImputations)){
  currDF <- compiledImputations[[i]]
  currDF$ImputNum <- i
  compiledImputations[[i]] <- currDF
}

# Load timestamps:
indexedTimes <- read.csv('~/data/all_motion_feature_data/indexed_times.csv') %>% mutate(times = as.POSIXct(times, format = "%d-%b-%Y %H:%M:%S",tz = "America/New_York"))

# Load complete bed-corrected, imputed feature set and add corresponding timestamps to accelerometry data:
for (imp in 1:length(compiledImputations)){
  currImp <- compiledImputations[[imp]] %>% select(-Bed, -ImputNum)
  currImp$TakenInstant <- ""
  for (patIdx in unique(indexedTimes$ptIdx)){
    currPtIdx <- currImp$ptIdx == patIdx
    currTimes <- indexedTimes %>% filter(ptIdx == patIdx) %>% rename(TakenInstant = times) %>% mutate(timeCount = 1:(sum(currPtIdx)/7)) %>% slice(rep(1:n(), times = 7))
    currTimes$TakenInstant <- as.character(currTimes$TakenInstant)
    currImp[currPtIdx,c("TakenInstant","ptIdx","timeCount")] <- currTimes
  }
  currImp$TakenInstant <- as.POSIXct(currImp$TakenInstant, tz = "America/New_York")
  save(currImp,file = file.path('~/data/all_motion_feature_data/final_imputed_features',paste0('imp',imp,'.RData')))
}

rm(compiledImputations,currPtIdx,currTimes,indexedTimes)
gc()