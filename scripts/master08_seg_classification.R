#### Master Script 8: Classification of Temporal Segments ####
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
library(lolR)
library(R.matlab)
library(tidyverse)
library(ggplot2)
library(plotly)
library(naniar)
library(MASS)
library(glmnet)
library(caret)
library(kernlab)
#library(rlist)
library(cvAUC)
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
#library(h2o)
if(.Platform$OS.type == "unix") {
  library(doMC)
} else {
  library(doParallel)
}

# Initiate H20 Cluster:
#h2o.init()

# Set the number of parallel cores for parallel tuning
#no.parallel.cores <- floor(2 * detectCores() / 3)
no.parallel.cores <- floor(2 * detectCores() / 3)
if(.Platform$OS.type == "unix") {
  registerDoMC(cores = no.parallel.cores)
} else {
  registerDoParallel(cores = no.parallel.cores)
}

source('./functions/load_patient_clinical_data.R')
source('./functions/seg_classification_function.R')

patient_clinical_data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv') %>% add_column(ptIdx = 1:nrow(.),.after = 2)

gose_thresh <- 5 # inclusive threshold
mrs_thresh <- 3 # inclusive threshold
# Code outcomes as favorable vs. unfavorable
patient_clinical_data <- mutate(patient_clinical_data,
                                fav_mort_dis = factor(DiedDuringThisHospitalStay_ == "0",labels = c("Unfav","Fav")),
                                fav_mort_12m = factor(Death12Months == "0",labels = c("Unfav","Fav")),
                                fav_GOSE_dis = factor(GOS_EDischarge >= gose_thresh,labels = c("Unfav","Fav")),
                                fav_GOSE_12m = factor(GOS_E12Months >= gose_thresh,labels = c("Unfav","Fav")),
                                fav_mRS_dis = factor(mRSDischarge <= mrs_thresh,labels = c("Unfav","Fav")),
                                fav_mRS_12m = factor(mRS12Months <= mrs_thresh,labels = c("Unfav","Fav")))

mf_choice<-c("band_power","freq_entropy","freq_pairs1","freq_pairs2","med_freq","sma","wavelets")
sr_choice<-c("LA","LE","LW","RA","RE","RW")
#classifier_choice<-c("svmRadialWeights","knn","lda","glmnet","avNNet","parRF","adaboost")
classifier_choice<-c("lda","glmnet_h2o","avNNet","parRF")

outcome <- "fav_GOSE_dis"

outer_fold_count <- 5
inner_fold_count <- 5

for (seg_windows in c(5,10,30,60,180)){
  
  print(paste('Segment size:',seg_windows,"min started"))
  
  seg_dir <- list.dirs('../all_motion_feature_data/')[str_detect(list.dirs('../all_motion_feature_data/'),pattern=as.character(seg_windows))]
  seg_files <- list.files(seg_dir,pattern = "*.csv",full.names = TRUE)
  
#  if (!exists("compiled_imps")) {
#  compiled_imps <- lapply(seg_files,read_csv) %>% lapply(.,function(x){x %>% dplyr::select(-1,-2)})
#  }

  curr_imp <-
    read_csv(seg_files[1]) %>% dplyr::select(-1, -2) %>% 
    filter(featureType %in% mf_choice) %>% 
    dplyr::select(all_of(c(
      sr_choice, "ptIdx", "segTimePoint", "featureType", "segIdx"
    )))
  
  #clinical_vars <- c("NCCUAdmissionAge_y_","Temp","MAP","HR","RespRate","Na","K","")
  mf_col_idx <- seq(3,length(sr_choice)*length(mf_choice)*seg_windows*12+2)

  dir.create(file.path('..','all_motion_feature_data','ML_predictions',outcome),showWarnings = FALSE)
  dir.create(file.path('..','all_motion_feature_data','ML_predictions',outcome,paste0(seg_windows,'_min')),showWarnings = FALSE)
  
  # Preparing covariates for ML training
  #for (l in 1:length(compiled_imps)){
  l=1
  #curr_imp <- compiled_imps[[l]] %>% 
  curr_data_table<-dcast(melt(curr_imp,id.vars = c("ptIdx","featureType","segIdx","segTimePoint")),ptIdx+segIdx~featureType+variable+segTimePoint) %>% full_join(.,patient_clinical_data,by="ptIdx")
  if (exists("curr_imp")) {
    rm(curr_imp)
    gc()
  }
  outList<-seg_classification_function(curr_data_table,outcome,classifier_choice,patient_clinical_data,mf_col_idx,outer_fold_count,inner_fold_count)
  save(outList,file=file.path('..','all_motion_feature_data','ML_predictions',outcome,paste0(seg_windows,'_min'),paste0('imp',l,'.RData')))
  gc()
  print(paste('Segment size:',seg_windows,"min complete"))
}