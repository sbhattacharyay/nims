#### Master Script 10: Grid Search of ML Predictors ####
# Decoding Quantitative Motor Features for Classification and Prediction
# in Severe Acquired Brain Injury
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

if(.Platform$OS.type == "unix") {
  setwd('~/work/SB_bims/scripts/')  
  .libPaths( c( "~/Rpackages" , .libPaths(),"/tmp/Rtmpcc66n8/downloaded_packages" ) )
}

if (!require(lolR))
  devtools::install_github('neurodata/lol',
                 build_vignettes = FALSE,
                 force = TRUE,
                 out_dir='~/Rpackages/')
library(lolR)
library(R.matlab)
# library(dplyr)
# library(tidyr)
# library(readr)
# library(purrr)
# library(tibble)
# library(stringr)
# library(forcats)
# library(ggplot2)
library(readxl)
library(plotly)
library(naniar)
library(MASS)
library(tidyverse)
library(glmnet)
library(caret)
library(kernlab)
# library(rlist)
library(cvAUC)
library(gridExtra)
library(ggplotify)
library(grid)
library(cowplot)
library(nnet)
library(parallel)
library(fastAdaboost)
library(cutpointr)
if(.Platform$OS.type == "unix") {
  library(doMC)
} else {
  library(doParallel)
}

source('./functions/list_cbind.R')
source('./functions/get_motion_features.R')
source('./functions/lol_project_motion_features.R')
source("./functions/load_tf_patient_covariates.R")
source('./functions/cv_lol_project_motion_features.R')
source('./functions/prepare_training_covariates.R')
source('./functions/prepare_testing_covariates.R')
source('./functions/train_caret_models.R')
source('./functions/predict_caret_models.R')
source('./functions/classification_function.R')
source('./functions/get_sens_info.R')
source('./functions/get_spec_info.R')
source('./functions/get_PPV_info.R')
source('./functions/get_NPV_info.R')
source('./functions/get_acc_info.R')
source('./functions/get_auc_info.R')
source('./functions/get_auc_cp.R')
source('./functions/get_auc_plots.R')
source('./functions/get_auc_info_ci.R')
source('./functions/generateRootDir.R')
source('./functions/get_cutpoint.R')
source('./functions/plot_metric_stripcharts.R')

# Set the number of parallel cores for parallel tuning
no.parallel.cores <- floor(2 * detectCores() / 3)
if(.Platform$OS.type == "unix") {
  registerDoMC(cores = no.parallel.cores)
} else {
  registerDoParallel(cores = no.parallel.cores)
}

# Load Motion Features Organized by Time of Day (TOD)
if (!exists("tod_motion_features")) {
  tod_sensors <-
    readMat('~/data/tod_motion_feature_data/bed_corrected_imputed_complete_sensor_data.mat')$bed.corrected.sensors
  tod_motion_features <- get_motion_features(do.call(rbind, tod_sensors))
}

# Load Motion Features Organized by Time from Recording (TFR)
if (!exists("tfr_motion_features")) {
  tfr_sensors <-
    readMat('~/data/tfr_motion_feature_data/bed_corrected_imputed_complete_sensor_data.mat')$bed.corrected.sensors
  tfr_motion_features <- get_motion_features(do.call(rbind, tfr_sensors))
}

# Create all combinations of motion feature choices:
mf_names <- names(tod_motion_features)
mf_options <- list()

count <- 0
for (i in 1:length(mf_names)){
  curr_comb<- combn(mf_names,i)
  for (j in 1:ncol(curr_comb)){
    count <- count + 1
    mf_options[[count]] <- curr_comb[,j]
  }
}

# Create all combinations of sensor placement choices:
sp_options <- list()

count <- 0
for (i in 1:length(sensor_options)){
  curr_comb<- combn(sensor_options,i)
  for (j in 1:ncol(curr_comb)){
    count <- count + 1
    sp_options[[count]] <- curr_comb[,j]
  }
}

preds_list <- list(mf_options,sp_options,time_choice=c('tod','tfr'),outcomes = c('dis','12m'))
complete_predictor_set <- expand.grid(preds_list)
names(complete_predictor_set) <- c("mf_options","sp_options","time_choice","outcomes")

classifier_choice<-c("svmRadialWeights","knn","lda","glmnet","avNNet","parRF")

out_complete_pred_set <- list()
for (i in 1:nrow(complete_predictor_set)){
  curr_timeChoice <- as.character(complete_predictor_set$time_choice[i])
  if (curr_timeChoice == 'tod'){
    curr_timeSlide <- c(as.POSIXct("2020-05-05 18:00:00"),as.POSIXct("2020-05-06 12:00:00"))
  } else {
    curr_timeSlide <- c(0,8)
  }
  curr_outcome_code <- as.character(complete_predictor_set$outcomes[i])
  if (curr_outcome_code == 'dis'){
    curr_outcomes <- patient_clinical_data$DiedDuringThisHospitalStay_
  } else {
    curr_outcomes <- patient_clinical_data$Death12Months
  }
  curr_mfChoice <- complete_predictor_set$mf_options[[i]]
  curr_sensorLoc <- complete_predictor_set$sp_options[[i]]
  
  out_complete_pred_set[[i]]<-classification_function(curr_timeChoice, 
                                                      curr_timeSlide, 
                                                      classifier_choice, 
                                                      curr_outcomes, 
                                                      curr_mfChoice, 
                                                      curr_sensorLoc)
}

