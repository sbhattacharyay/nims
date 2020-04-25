#### Master Script 9: GLM Classification ####
# Decoding Quantitative Motor Features for Classification and Prediction
# in Severe Acquired Brain Injury
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

library(devtools)
if (!require(lolR)) install_github('neurodata/lol', build_vignettes=TRUE, force=TRUE)
library(lolR)
library(R.matlab)
library(tidyverse)
library(ggplot2)
library(plotly)
library(naniar)
library(MASS)

source('./functions/load_patient_clinical_data.R')
source('./functions/update_clinicalVariableList.R')
source('./functions/get_motion_features.R')
source('./functions/lol_project_motion_features.R')
source('./functions/viz_lol_2D.R')
source("./functions/generateRootDir.R")
source("./functions/cross_val_splits.R")
source("./functions/load_tf_patient_covariates.R")
source('./functions/cv_lol_project_motion_features.R')
source('./functions/prepare_training_covariates.R')
source('./functions/prepare_testing_covariates.R')

# Load patient clinical data (sorts by PY numbering and corrects variable types)
patient_clinical_data<-load_patient_clinical_data()

# Load and update clinical variable list
clinicalVariableList<-update_clinicalVariableList()

# Load Motion Features Organized by Time of Day (TOD)
if(!exists("tod_sensors")) {
  tod_sensors<-readMat('../motion_feature_data/bed_corrected_imputed_complete_sensor_data.mat')$bed.corrected.sensors
}

# Load Motion Features Organized by Time from Recording (TFR)
if(!exists("tfr_sensors")) {
  tfr_sensors<-readMat('../tfr_motion_feature_data/bed_corrected_imputed_complete_sensor_data.mat')$bed.corrected.sensors
}

# Organizing Features for both motion feature categorizations
tod_temp_var<-do.call(rbind,tod_sensors)
tfr_temp_var<-do.call(rbind,tfr_sensors)
# ranked_temp_var<-lapply(temp_var, function(x) t(apply(x,1,sort,decreasing = TRUE)))
# fft_temp_var<-lapply(temp_var, function(x) abs(t(mvfft(t(x)))))

tod_motion_features <- get_motion_features(tod_temp_var)
tfr_motion_features <- get_motion_features(tfr_temp_var)

# Load transformed covariates for both TOD and TFR
tod_tf_covariates<-load_tf_patient_covariates('../motion_feature_data/tf_patient_covariates.csv')
tfr_tf_covariates<-load_tf_patient_covariates('../tfr_motion_feature_data/tf_patient_covariates.csv')

# Split for k-fold cross validation for discharge predictions:
k<-5
cvIdx<-cross_val_splits(patient_clinical_data,k)

for (i in 1:length(cvIdx)){
  currTestIdx<-cvIdx[[as.character(i)]]
  currTrainIdx<-seq(nrow(patient_clinical_data))[-currTestIdx]
  
  # Convert numeric to logical indexing
  logicalTotal <- rep(FALSE,nrow(patient_clinical_data))
  logicalTest <- logicalTotal
  logicalTest[currTestIdx] <- TRUE
  logicalTrain <- !logicalTest

  # Extract training and testing covariates based on current split
  train_tfr_tf_covariates<-tfr_tf_covariates[currTrainIdx,]
  test_tfr_tf_covariates<-tfr_tf_covariates[currTestIdx,]
  
  train_tod_tf_covariates<-tod_tf_covariates[currTrainIdx,]
  test_tod_tf_covariates<-tod_tf_covariates[currTestIdx,]
  
  # Perform LOL on training data
  r<-3
  Y<-as.factor(patient_clinical_data$gose)
  tod_GOSE_LOL<-cv_lol_project_motion_features(tod_motion_features,Y,r,logicalTrain)
  tfr_GOSE_LOL<-cv_lol_project_motion_features(tfr_motion_features,Y,r,logicalTrain)

  Y<-as.factor(patient_clinical_data$favorable)
  tod_fav_LOL<-cv_lol_project_motion_features(tod_motion_features,Y,r,logicalTrain)
  tfr_fav_LOL<-cv_lol_project_motion_features(tfr_motion_features,Y,r,logicalTrain)
  
  Y<-as.factor(patient_clinical_data$death)
  tod_death_LOL<-cv_lol_project_motion_features(tod_motion_features,Y,r,logicalTrain)
  tfr_death_LOL<-cv_lol_project_motion_features(tfr_motion_features,Y,r,logicalTrain)
  
  # Prepare training covariates
  train_tod_GOSE_CV<- prepare_training_covariates(tod_GOSE_LOL,train_tod_tf_covariates,TRUE)
  train_tod_fav_CV<-  prepare_training_covariates(tod_fav_LOL,train_tod_tf_covariates,TRUE)
  train_tod_death_CV<-prepare_training_covariates(tod_death_LOL,train_tod_tf_covariates,FALSE)
  
  train_tfr_GOSE_CV<- prepare_training_covariates(tfr_GOSE_LOL,train_tfr_tf_covariates,TRUE)
  train_tfr_fav_CV<-  prepare_training_covariates(tfr_fav_LOL,train_tfr_tf_covariates,TRUE)
  train_tfr_death_CV<-prepare_training_covariates(tfr_death_LOL,train_tfr_tf_covariates,FALSE)
  
  # Prepare testing covariates
  test_tod_GOSE_CV<- prepare_testing_covariates(tod_motion_features,tod_GOSE_LOL,test_tod_tf_covariates,currTestIdx,TRUE)
  test_tod_fav_CV<-  prepare_testing_covariates(tod_motion_features,tod_fav_LOL,test_tod_tf_covariates,currTestIdx,TRUE)
  test_tod_death_CV<-prepare_testing_covariates(tod_motion_features,tod_death_LOL,test_tod_tf_covariates,currTestIdx,FALSE)
  
  test_tfr_GOSE_CV<- prepare_testing_covariates(tfr_motion_features,tfr_GOSE_LOL,test_tfr_tf_covariates,currTestIdx,TRUE)
  test_tfr_fav_CV<-  prepare_testing_covariates(tfr_motion_features,tfr_fav_LOL,test_tfr_tf_covariates,currTestIdx,TRUE)
  test_tfr_death_CV<-prepare_testing_covariates(tfr_motion_features,tfr_death_LOL,test_tfr_tf_covariates,currTestIdx,FALSE)
  
  # Train GLM on training data
  Y<-as.factor(patient_clinical_data$favorable)
  tod_GOSE_LOL<-cv_lol_project_motion_features(tod_motion_features,Y,r,logicalTrain)
  tfr_GOSE_LOL<-cv_lol_project_motion_features(tfr_motion_features,Y,r,logicalTrain)
  
  tod_fav_LOL<-cv_lol_project_motion_features(tod_motion_features,Y,r,logicalTrain)
  tfr_fav_LOL<-cv_lol_project_motion_features(tfr_motion_features,Y,r,logicalTrain)
  
  Y<-as.factor(patient_clinical_data$death)
  tod_death_LOL<-cv_lol_project_motion_features(tod_motion_features,Y,r,logicalTrain)
  tfr_death_LOL<-cv_lol_project_motion_features(tfr_motion_features,Y,r,logicalTrain)
}

# Split for k-fold cross validation for 12 month predictions:
k<-5
cvIdx<-cross_val_splits(patient_clinical_data,k)


Y<-as.factor(patient_clinical_data$death_12mo)
tod_death12mo_LOL<-cv_lol_project_motion_features(tod_motion_features,Y,r,logicalTrain)
tfr_death12mo_LOL<-cv_lol_project_motion_features(tfr_motion_features,Y,r,logicalTrain)


Y<-as.factor(patient_clinical_data$favorable_12mo)
tod_fav12mo_LOL<-cv_lol_project_motion_features(tod_motion_features,Y,r,logicalTrain)
tfr_fav12mo_LOL<-cv_lol_project_motion_features(tfr_motion_features,Y,r,logicalTrain)


Y<-as.factor(patient_clinical_data$gose_12mo)
tod_GOSE12mo_LOL<-cv_lol_project_motion_features(tod_motion_features,r,3,logicalTrain)
tfr_GOSE12mo_LOL<-cv_lol_project_motion_features(tfr_motion_features,r,3,logicalTrain)
