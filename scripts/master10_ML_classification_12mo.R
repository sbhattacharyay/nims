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
source('./functions/train_GLM.R')
source('./functions/predict_GLM.R')

# Load patient clinical data (sorts by PY numbering and corrects variable types)
patient_clinical_data <- load_patient_clinical_data()

# Load and update clinical variable list
clinicalVariableList <- update_clinicalVariableList()

# Load Motion Features Organized by Time of Day (TOD)
if (!exists("tod_sensors")) {
  tod_sensors <-
    readMat('../motion_feature_data/bed_corrected_imputed_complete_sensor_data.mat')$bed.corrected.sensors
}

# Load Motion Features Organized by Time from Recording (TFR)
if (!exists("tfr_sensors")) {
  tfr_sensors <-
    readMat('../tfr_motion_feature_data/bed_corrected_imputed_complete_sensor_data.mat')$bed.corrected.sensors
}

# Organizing Features for both motion feature categorizations
tod_temp_var <- do.call(rbind, tod_sensors)
tfr_temp_var <- do.call(rbind, tfr_sensors)
# ranked_temp_var<-lapply(temp_var, function(x) t(apply(x,1,sort,decreasing = TRUE)))
# fft_temp_var<-lapply(temp_var, function(x) abs(t(mvfft(t(x)))))

tod_motion_features <- get_motion_features(tod_temp_var)
tfr_motion_features <- get_motion_features(tfr_temp_var)

# Load transformed covariates for both TOD and TFR
tod_tf_covariates <-
  load_tf_patient_covariates('../motion_feature_data/tf_patient_covariates.csv')
tfr_tf_covariates <-
  load_tf_patient_covariates('../tfr_motion_feature_data/tf_patient_covariates.csv')

# Split for k-fold cross validation for discharge predictions:
k <- 5
cvIdx <- cross_val_splits(patient_clinical_data, k)

self_val_predictionsDis <- vector(mode = "list", length = k)
test_val_predictionsDis <- vector(mode = "list", length = k)

for (i in 1:length(cvIdx)) {
  currTestIdx <- cvIdx[[as.character(i)]]
  currTrainIdx <- seq(nrow(patient_clinical_data))[-currTestIdx]
  
  # Convert numeric to logical indexing
  logicalTotal <- rep(FALSE, nrow(patient_clinical_data))
  logicalTest <- logicalTotal
  logicalTest[currTestIdx] <- TRUE
  logicalTrain <- !logicalTest
  
  # Extract training and testing covariates based on current split
  train_tfr_tf_covariates <- tfr_tf_covariates[currTrainIdx,]
  test_tfr_tf_covariates <- tfr_tf_covariates[currTestIdx,]
  
  train_tod_tf_covariates <- tod_tf_covariates[currTrainIdx,]
  test_tod_tf_covariates <- tod_tf_covariates[currTestIdx,]
  
  # Perform LOL on training data
  r <- 3
  Y <- as.factor(patient_clinical_data$gose)
  tod_GOSE_LOL <-
    cv_lol_project_motion_features(tod_motion_features, Y, r, logicalTrain)
  tfr_GOSE_LOL <-
    cv_lol_project_motion_features(tfr_motion_features, Y, r, logicalTrain)
  
  Y <- as.factor(patient_clinical_data$favorable)
  tod_fav_LOL <-
    cv_lol_project_motion_features(tod_motion_features, Y, r, logicalTrain)
  tfr_fav_LOL <-
    cv_lol_project_motion_features(tfr_motion_features, Y, r, logicalTrain)
  
  Y <- as.factor(patient_clinical_data$death)
  tod_death_LOL <-
    cv_lol_project_motion_features(tod_motion_features, Y, r, logicalTrain)
  tfr_death_LOL <-
    cv_lol_project_motion_features(tfr_motion_features, Y, r, logicalTrain)
  
  # Prepare training covariates
  train_tod_GOSE_CV <-
    prepare_training_covariates(tod_GOSE_LOL, train_tod_tf_covariates, TRUE)
  train_tod_fav_CV <-
    prepare_training_covariates(tod_fav_LOL, train_tod_tf_covariates, TRUE)
  train_tod_death_CV <-
    prepare_training_covariates(tod_death_LOL, train_tod_tf_covariates, FALSE)
  
  train_tfr_GOSE_CV <-
    prepare_training_covariates(tfr_GOSE_LOL, train_tfr_tf_covariates, TRUE)
  train_tfr_fav_CV <-
    prepare_training_covariates(tfr_fav_LOL, train_tfr_tf_covariates, TRUE)
  train_tfr_death_CV <-
    prepare_training_covariates(tfr_death_LOL, train_tfr_tf_covariates, FALSE)
  
  # Prepare testing covariates
  test_tod_GOSE_CV <-
    prepare_testing_covariates(tod_motion_features,
                               tod_GOSE_LOL,
                               test_tod_tf_covariates,
                               currTestIdx,
                               TRUE)
  test_tod_fav_CV <-
    prepare_testing_covariates(tod_motion_features,
                               tod_fav_LOL,
                               test_tod_tf_covariates,
                               currTestIdx,
                               TRUE)
  test_tod_death_CV <-
    prepare_testing_covariates(tod_motion_features,
                               tod_death_LOL,
                               test_tod_tf_covariates,
                               currTestIdx,
                               FALSE)
  
  test_tfr_GOSE_CV <-
    prepare_testing_covariates(tfr_motion_features,
                               tfr_GOSE_LOL,
                               test_tfr_tf_covariates,
                               currTestIdx,
                               TRUE)
  test_tfr_fav_CV <-
    prepare_testing_covariates(tfr_motion_features,
                               tfr_fav_LOL,
                               test_tfr_tf_covariates,
                               currTestIdx,
                               TRUE)
  test_tfr_death_CV <-
    prepare_testing_covariates(tfr_motion_features,
                               tfr_death_LOL,
                               test_tfr_tf_covariates,
                               currTestIdx,
                               FALSE)
  
  # Train GLM on training data
  Y <- as.factor(patient_clinical_data$favorable)
  
  tod_GOSE_GLM <-
    train_GLM(train_tod_GOSE_CV, Y, currTrainIdx, 0, 4)
  tfr_GOSE_GLM <-
    train_GLM(train_tfr_GOSE_CV, Y, currTrainIdx, 0, 4)
  
  tod_fav_GLM <- train_GLM(train_tod_fav_CV, Y, currTrainIdx, 0, 4)
  tfr_fav_GLM <- train_GLM(train_tfr_fav_CV, Y, currTrainIdx, 0, 4)
  
  Y <- as.factor(patient_clinical_data$death)
  
  tod_death_GLM <-
    train_GLM(train_tod_death_CV, Y, currTrainIdx, 0, 4)
  tfr_death_GLM <-
    train_GLM(train_tfr_death_CV, Y, currTrainIdx, 0, 4)
  
  # Self-validate the GLM models on training data
  
  Y <- as.factor(patient_clinical_data$favorable)
  
  tod_GOSE_self_val <-
    predict_GLM(tod_GOSE_GLM, train_tod_GOSE_CV, Y, currTrainIdx, TRUE)
  tfr_GOSE_self_val <-
    predict_GLM(tfr_GOSE_GLM, train_tfr_GOSE_CV, Y, currTrainIdx, TRUE)
  
  tod_fav_self_val <-
    predict_GLM(tod_fav_GLM, train_tod_fav_CV, Y, currTrainIdx, TRUE)
  tfr_fav_self_val <-
    predict_GLM(tfr_fav_GLM, train_tfr_fav_CV, Y, currTrainIdx, TRUE)
  
  Y <- as.factor(patient_clinical_data$death)
  
  tod_death_self_val <-
    predict_GLM(tod_death_GLM, train_tod_death_CV, Y, currTrainIdx, TRUE)
  tfr_death_self_val <-
    predict_GLM(tfr_death_GLM, train_tfr_death_CV, Y, currTrainIdx, TRUE)
  
  temp_self_val_list <-
    list(
      tod_GOSE_self_val,
      tfr_GOSE_self_val,
      tod_fav_self_val,
      tfr_fav_self_val,
      tod_death_self_val,
      tfr_death_self_val
    )
  names(temp_self_val_list) <-
    c(
      "tod_GOSE_self_val",
      "tfr_GOSE_self_val",
      "tod_fav_self_val",
      "tfr_fav_self_val",
      "tod_death_self_val",
      "tfr_death_self_val"
    )
  self_val_predictionsDis[[i]] <- temp_self_val_list
  
  # Predict outcomes for testing data
  
  Y <- as.factor(patient_clinical_data$favorable)
  
  tod_GOSE_test_val <-
    predict_GLM(tod_GOSE_GLM, test_tod_GOSE_CV, Y, currTestIdx, FALSE)
  tfr_GOSE_test_val <-
    predict_GLM(tfr_GOSE_GLM, test_tfr_GOSE_CV, Y, currTestIdx, FALSE)
  
  tod_fav_test_val <-
    predict_GLM(tod_fav_GLM, test_tod_fav_CV, Y, currTestIdx, FALSE)
  tfr_fav_test_val <-
    predict_GLM(tfr_fav_GLM, test_tfr_fav_CV, Y, currTestIdx, FALSE)
  
  Y <- as.factor(patient_clinical_data$death)
  
  tod_death_test_val <-
    predict_GLM(tod_death_GLM, test_tod_death_CV, Y, currTestIdx, FALSE)
  tfr_death_test_val <-
    predict_GLM(tfr_death_GLM, test_tfr_death_CV, Y, currTestIdx, FALSE)
  
  temp_test_val_list <-
    list(
      tod_GOSE_test_val,
      tfr_GOSE_test_val,
      tod_fav_test_val,
      tfr_fav_test_val,
      tod_death_test_val,
      tfr_death_test_val
    )
  names(temp_test_val_list) <-
    c(
      "tod_GOSE_test_val",
      "tfr_GOSE_test_val",
      "tod_fav_test_val",
      "tfr_fav_test_val",
      "tod_death_test_val",
      "tfr_death_test_val"
    )
  test_val_predictionsDis[[i]] <- temp_test_val_list
}

# Split for k-fold cross validation for 12 month predictions:
k <- 5
idx_for_12mo <- !is.na(patient_clinical_data$favorable_12mo)
seq_for_12mo <- seq(62)[idx_for_12mo]
cvIdx <- cross_val_splits(patient_clinical_data[idx_for_12mo, ], k)

self_val_predictions12mo <- vector(mode = "list", length = k)
test_val_predictions12mo <- vector(mode = "list", length = k)

for (i in 1:length(cvIdx)) {
  currTestIdx <- cvIdx[[as.character(i)]]
  currTrainIdx <- seq_for_12mo[-currTestIdx]
  
  # Convert numeric to logical indexing
  logicalTotal <- rep(FALSE, nrow(patient_clinical_data))
  
  logicalTest <- logicalTotal
  logicalTest[currTestIdx] <- TRUE
  
  logicalTrain <- logicalTotal
  logicalTrain[currTrainIdx] <- TRUE
  
  # Extract training and testing covariates based on current split
  train_tfr_tf_covariates <- tfr_tf_covariates[currTrainIdx,]
  test_tfr_tf_covariates <- tfr_tf_covariates[currTestIdx,]
  
  train_tod_tf_covariates <- tod_tf_covariates[currTrainIdx,]
  test_tod_tf_covariates <- tod_tf_covariates[currTestIdx,]
  
  # Perform LOL on training data
  r <- 3
  Y <- as.factor(patient_clinical_data$gose_12mo)
  tod_GOSE_LOL <-
    cv_lol_project_motion_features(tod_motion_features, Y, r, logicalTrain)
  tfr_GOSE_LOL <-
    cv_lol_project_motion_features(tfr_motion_features, Y, r, logicalTrain)
  
  Y <- as.factor(patient_clinical_data$favorable_12mo)
  tod_fav_LOL <-
    cv_lol_project_motion_features(tod_motion_features, Y, r, logicalTrain)
  tfr_fav_LOL <-
    cv_lol_project_motion_features(tfr_motion_features, Y, r, logicalTrain)
  
  Y <- as.factor(patient_clinical_data$death_12mo)
  tod_death_LOL <-
    cv_lol_project_motion_features(tod_motion_features, Y, r, logicalTrain)
  tfr_death_LOL <-
    cv_lol_project_motion_features(tfr_motion_features, Y, r, logicalTrain)
  
  # Prepare training covariates
  train_tod_GOSE_CV <-
    prepare_training_covariates(tod_GOSE_LOL, train_tod_tf_covariates, TRUE)
  train_tod_fav_CV <-
    prepare_training_covariates(tod_fav_LOL, train_tod_tf_covariates, TRUE)
  train_tod_death_CV <-
    prepare_training_covariates(tod_death_LOL, train_tod_tf_covariates, FALSE)
  
  train_tfr_GOSE_CV <-
    prepare_training_covariates(tfr_GOSE_LOL, train_tfr_tf_covariates, TRUE)
  train_tfr_fav_CV <-
    prepare_training_covariates(tfr_fav_LOL, train_tfr_tf_covariates, TRUE)
  train_tfr_death_CV <-
    prepare_training_covariates(tfr_death_LOL, train_tfr_tf_covariates, FALSE)
  
  # Prepare testing covariates
  test_tod_GOSE_CV <-
    prepare_testing_covariates(tod_motion_features,
                               tod_GOSE_LOL,
                               test_tod_tf_covariates,
                               currTestIdx,
                               TRUE)
  test_tod_fav_CV <-
    prepare_testing_covariates(tod_motion_features,
                               tod_fav_LOL,
                               test_tod_tf_covariates,
                               currTestIdx,
                               TRUE)
  test_tod_death_CV <-
    prepare_testing_covariates(tod_motion_features,
                               tod_death_LOL,
                               test_tod_tf_covariates,
                               currTestIdx,
                               FALSE)
  
  test_tfr_GOSE_CV <-
    prepare_testing_covariates(tfr_motion_features,
                               tfr_GOSE_LOL,
                               test_tfr_tf_covariates,
                               currTestIdx,
                               TRUE)
  test_tfr_fav_CV <-
    prepare_testing_covariates(tfr_motion_features,
                               tfr_fav_LOL,
                               test_tfr_tf_covariates,
                               currTestIdx,
                               TRUE)
  test_tfr_death_CV <-
    prepare_testing_covariates(tfr_motion_features,
                               tfr_death_LOL,
                               test_tfr_tf_covariates,
                               currTestIdx,
                               FALSE)
  
  # Train GLM on training data
  Y <- as.factor(patient_clinical_data$favorable_12mo)
  
  tod_GOSE_GLM <-
    train_GLM(train_tod_GOSE_CV, Y, currTrainIdx, 0, 4)
  tfr_GOSE_GLM <-
    train_GLM(train_tfr_GOSE_CV, Y, currTrainIdx, 0, 4)
  
  tod_fav_GLM <- train_GLM(train_tod_fav_CV, Y, currTrainIdx, 0, 4)
  tfr_fav_GLM <- train_GLM(train_tfr_fav_CV, Y, currTrainIdx, 0, 4)
  
  Y <- as.factor(patient_clinical_data$death_12mo)
  
  tod_death_GLM <-
    train_GLM(train_tod_death_CV, Y, currTrainIdx, 0, 4)
  tfr_death_GLM <-
    train_GLM(train_tfr_death_CV, Y, currTrainIdx, 0, 4)
  
  # Self-validate the GLM models on training data
  
  Y <- as.factor(patient_clinical_data$favorable_12mo)
  
  tod_GOSE_self_val <-
    predict_GLM(tod_GOSE_GLM, train_tod_GOSE_CV, Y, currTrainIdx, TRUE)
  tfr_GOSE_self_val <-
    predict_GLM(tfr_GOSE_GLM, train_tfr_GOSE_CV, Y, currTrainIdx, TRUE)
  
  tod_fav_self_val <-
    predict_GLM(tod_fav_GLM, train_tod_fav_CV, Y, currTrainIdx, TRUE)
  tfr_fav_self_val <-
    predict_GLM(tfr_fav_GLM, train_tfr_fav_CV, Y, currTrainIdx, TRUE)
  
  Y <- as.factor(patient_clinical_data$death_12mo)
  
  tod_death_self_val <-
    predict_GLM(tod_death_GLM, train_tod_death_CV, Y, currTrainIdx, TRUE)
  tfr_death_self_val <-
    predict_GLM(tfr_death_GLM, train_tfr_death_CV, Y, currTrainIdx, TRUE)
  
  temp_self_val_list <-
    list(
      tod_GOSE_self_val,
      tfr_GOSE_self_val,
      tod_fav_self_val,
      tfr_fav_self_val,
      tod_death_self_val,
      tfr_death_self_val
    )
  names(temp_self_val_list) <-
    c(
      "tod_GOSE_self_val",
      "tfr_GOSE_self_val",
      "tod_fav_self_val",
      "tfr_fav_self_val",
      "tod_death_self_val",
      "tfr_death_self_val"
    )
  self_val_predictions12mo[[i]] <- temp_self_val_list
  
  # Predict outcomes for testing data
  
  Y <- as.factor(patient_clinical_data$favorable_12mo)
  
  tod_GOSE_test_val <-
    predict_GLM(tod_GOSE_GLM, test_tod_GOSE_CV, Y, currTestIdx, FALSE)
  tfr_GOSE_test_val <-
    predict_GLM(tfr_GOSE_GLM, test_tfr_GOSE_CV, Y, currTestIdx, FALSE)
  
  tod_fav_test_val <-
    predict_GLM(tod_fav_GLM, test_tod_fav_CV, Y, currTestIdx, FALSE)
  tfr_fav_test_val <-
    predict_GLM(tfr_fav_GLM, test_tfr_fav_CV, Y, currTestIdx, FALSE)
  
  Y <- as.factor(patient_clinical_data$death_12mo)
  
  tod_death_test_val <-
    predict_GLM(tod_death_GLM, test_tod_death_CV, Y, currTestIdx, FALSE)
  tfr_death_test_val <-
    predict_GLM(tfr_death_GLM, test_tfr_death_CV, Y, currTestIdx, FALSE)
  
  temp_test_val_list <-
    list(
      tod_GOSE_test_val,
      tfr_GOSE_test_val,
      tod_fav_test_val,
      tfr_fav_test_val,
      tod_death_test_val,
      tfr_death_test_val
    )
  names(temp_test_val_list) <-
    c(
      "tod_GOSE_test_val",
      "tfr_GOSE_test_val",
      "tod_fav_test_val",
      "tfr_fav_test_val",
      "tod_death_test_val",
      "tfr_death_test_val"
    )
  test_val_predictions12mo[[i]] <- temp_test_val_list
}

rm(list = setdiff(
  ls(),
  c(
    "self_val_predictionsDis",
    "test_val_predictionsDis",
    "self_val_predictions12mo",
    "test_val_predictions12mo",
    "patient_clinical_data",
    "clinicalVariableList",
    "tod_motion_features",
    "tfr_motion_features"
  )
))