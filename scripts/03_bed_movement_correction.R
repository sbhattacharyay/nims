#### Master Script 3: Bed Motion Correction and Collection of Multiple Imputations ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# University of Cambridge
# Johns Hopkins University
# email address: sb2406@cam.ac.uk
#
### Contents:
# I. Initialization
# II. Load imputed motion feature data and collect into one object
# III. Determine feature space thresholds to correct bed motion 
# IV. Correct bed motion by identifying time points and patient indices during which SMA exceeds a literature-reviewed threshold

### I. Initialization
# Denote number of imputations
m <- 9

# Call requisite libraries
library(tidyverse)
library(readxl)

# Define feature label names
feature.labels <- c("BPW","FDE","HLF_h","HLF_l","MFR","SMA","WVL")

### II. Load imputed motion feature data and collect into one object
compiledImputations <- vector(mode = "list")
for (i in 1:m){
  print(paste('Imputation no.',i,'started'))
  currPattern <- paste0('*',i,'.csv')
  currFileList <- list.files('../features/01_imputed_features/',pattern = currPattern)
  imputation.df <- data.frame(matrix(ncol = 12, nrow = 0))
  for (j in 1:length(feature.labels)){
    currFilePath <- file.path('../features/01_imputed_features',currFileList[j]) 
    curr.df <- read.csv(currFilePath) %>% 
      mutate(Feature = feature.labels[j]) %>%
      relocate(Feature, .after = TimeOfDay)
    imputation.df <- rbind(imputation.df,curr.df)
    print(paste('Feature no.',j,'complete'))
  }
  compiledImputations[[i]] <- imputation.df %>% arrange(UPI,RecordingIdx,Feature)
  print(paste('Imputation no.',i,'complete'))
}

### III. Determine feature space thresholds to correct bed motion 
# Based on an SMA threshold for dynamic vs. static activity (https://doi.org/10.1016/j.medengphy.2013.06.005), 
# find corresponding thresholds for the other feature spaces. 
SMA.thresh <- .135

# Note we only use the first imputation to do so since, when testing, we found that the thresholds are the same for all imputations
source('./functions/find_thresholds.R')
feature.thresholds <- find_thresholds(compiledImputations[[1]],SMA.thresh)
featRanges <-
  list(
    BPW.nm.range = c(0, feature.thresholds[1]),
    FDE.nm.range = c(feature.thresholds[2], 1.707),
    HLF_l.nm.range = c(0, feature.thresholds[3]),
    HLF_h.nm.range = c(0, feature.thresholds[4]),
    MFR.nm.range = c(feature.thresholds[5], 3.2),
    SMA.nm.range = c(0, feature.thresholds[6]),
    WVL.nm.range = c(0, feature.thresholds[7])
  )

### IV. Correct bed motion by identifying time points and patient indices during which SMA exceeds a literature-reviewed threshold:
for (i in 1:length(compiledImputations)){
  print(paste("Imputation No.",i,"Started"))
  curr.df <- compiledImputations[[i]]
  curr.SMA.df <-curr.df %>% filter(Feature == "SMA")
  bed.SMA.rows <- which(curr.SMA.df$Bed > SMA.thresh)
  for (j in 1:length(feature.labels)){
    a <- featRanges[[j]][1]
    b <- featRanges[[j]][2]
    print(paste("Feature No.",j,"Started"))
    curr.feat.rows <- which(curr.df$Feature == feature.labels[j])
    currFeat.df <- curr.df %>% filter(Feature == feature.labels[j])
    currFeatChange.df <- currFeat.df[bed.SMA.rows,]
    if (feature.labels[j] %in% c("FDE","MFR")) {
      temp.mat <- currFeatChange.df[,c('LA','LE','LW','RA','RE','RW')] + currFeatChange.df[,'Bed']
      temp.mat[temp.mat > b] <- runif(sum(temp.mat > b),a+(b-a)/2,b)
    } else {
      temp.mat <- currFeatChange.df[,c('LA','LE','LW','RA','RE','RW')] - currFeatChange.df[,'Bed']
      temp.mat[temp.mat < 0] <- runif(sum(temp.mat < 0),0,b/2)
    }
    currFeatChange.df[,c('LA','LE','LW','RA','RE','RW')] <- temp.mat
    currFeat.df[bed.SMA.rows,] <- currFeatChange.df
    curr.df[curr.feat.rows,] <- currFeat.df
    print(paste("Feature No.",j,"Complete"))
  }
  curr.df <- curr.df %>% 
    mutate(ImputationNo = i) %>%
    relocate(ImputationNo,
             UPI,
             RecordingIdx,
             HoursFromICUAdmission,
             TimeOfDay,
             Feature)
  write.csv(currImp,paste0('../features/02_bed_corrected_imputed_features/bed_corrected_imputation_',i,'.csv'),row.names = F)
  print(paste("Imputation No.",i,"Complete"))
}
