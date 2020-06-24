#### Master Script 5: Temporal Segmentation, Dimensionality Reduction, and Predictor Preparation ####
# Decoding Quantitative Motor Features for Classification and Prediction
# in Severe Acquired Brain Injury
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

.libPaths( c( "~/Rpackages" , .libPaths(),"/tmp/Rtmpcc66n8/downloaded_packages" ) )
setwd('~/work/nims/scripts/')

library(devtools)
if (!require(lolR))
  install_github('neurodata/lol',
                 build_vignettes = TRUE,
                 force = TRUE)
library(tidyverse)
library(lolR)
library(R.matlab)
library(ggplot2)
library(plotly)
library(naniar)
library(MASS)
library(glmnet)
library(caret)
library(keras)
library(caTools)
library(kernlab)
library(pROC)
library(gridExtra)
library(ggplotify)
library(grid)
library(cowplot)
library(nnet)
library(parallel)
library(fastAdaboost)
library(cutpointr)
library(precrec)
library(reshape2)
library(MLmetrics)
library(DMwR)
if(.Platform$OS.type == "unix") {
  library(doMC)
} else {
  library(doParallel)
}
# Set the number of parallel cores
no.parallel.cores <- floor(2 * detectCores() / 3)
if(.Platform$OS.type == "unix") {
  registerDoMC(cores = no.parallel.cores)
} else {
  registerDoParallel(cores = no.parallel.cores)
}

if (!exists("compiledImputations")) {
  compiledImputations <- readRDS('../all_motion_feature_data/bc_i_c_dataset.rds')
}

for (seg_window in c(5,10,30,60,180)) {
  print(paste("Segment Window Size",seg_window,"min Started"))
  
  new_dir_name <- file.path("..","all_motion_feature_data",paste0("final_processed_features_",seg_window,"min"))
  
  n <- length(unique(compiledImputations[[1]]$ptIdx))
  
  # Algorithm strategy: assign new label to data that corresponds to chunk Index and provide a timeIdx within chunk
  dir.create(new_dir_name,showWarnings = FALSE)
  for (l in 1:length(compiledImputations)){
    print(paste("Imputation No.",l,"Started"))
    currImp <- compiledImputations[[l]]
    newImp <- as.data.frame(matrix(nrow = 0,ncol = 12))
    for (i in 1:n) {
      curr_patient <- currImp %>% filter(ptIdx == i)
      currPtTimePoints <- unique(curr_patient$timeCount)
      no_splits <- floor(length(currPtTimePoints)/(seg_window*12))
      truncTime <- currPtTimePoints[1:((seg_window*12)*no_splits)]
      indexSplits <- split(truncTime, ceiling(seq_along(truncTime)/(seg_window*12)))
      indexDF <- as.data.frame(indexSplits) %>% gather() %>% rename(segIdx=1,timeCount=2) %>% mutate(segIdx = gsub("[^0-9.-]", "", segIdx)) %>% mutate(segIdx=as.numeric(segIdx))
      indexDF$segTimePoint <- rep(1:length(indexSplits[[1]]),length(indexSplits))
      curr_patient <- curr_patient %>% inner_join(.,indexDF,by="timeCount") %>% dplyr::select(-timeCount)
      newImp <- rbind(newImp,curr_patient)
      print(paste("Patient No.",i,"Complete"))
    }
    imp_file_path <- file.path(new_dir_name,paste0("imp",l,".csv"))
    write.csv(newImp,imp_file_path)
    print(paste("Imputation No.",l,"Finished and Saved"))
  }
  print(paste("Segment Window Size",seg_window,"min Complete"))
}

### Load clinical metadata ###
source('./functions/load_patient_clinical_data.R')
gose_thresh <- 5 # inclusive threshold
mrs_thresh <- 3 # inclusive threshold
patient_clinical_data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv',gose_thresh,mrs_thresh) %>% add_column(ptIdx = 1:nrow(.),.after = 2)

### Create and save outer folds for ML based on available outcome labels ###
labels.temp <- expand.grid(label=c("fav_mort","fav_GOSE","fav_mRS"),
                           temp = c("dis","12m")) %>% 
  mutate(label.name = paste(label,temp,sep = "_"))

outer_fold_count <- 5
inner_fold_count <- 5

# source('./functions/createOuterCVFolds.R')
# outerFolds <- createOuterCVFolds(patient_clinical_data,outer_fold_count)
# saveRDS(outerFolds,'../all_motion_feature_data/outerFolds.rds')
# rm(outerFolds)
outerFolds <- readRDS('../all_motion_feature_data/outerFolds.rds')

### Write LOL-embedded motion features for training and testing ###
# Choose sensors and motion features to test under a given temporal segment window:
source('./functions/motion_feature_preparation.R')

# 5, 10, 30, 60, or 180 minutes
for (seg_window in c(5,10,30,60,180)){
  mf_choice<-c("band_power","freq_entropy","freq_pairs1","freq_pairs2","med_freq","sma","wavelets")
  sr_choice<-c("LA","LE","LW","RA","RE","RW")
  motion_feature_preparation(
    patient_clinical_data = patient_clinical_data,
    seg_window = seg_window,
    mf_choice = mf_choice,
    sr_choice = sr_choice,
    outerFolds = outerFolds,
    impNo = 1,
    saveDir = file.path('../all_motion_feature_data/complete_LOL_set',paste0(seg_window,"_min"))
  )
}