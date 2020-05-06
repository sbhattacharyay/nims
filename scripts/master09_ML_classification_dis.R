#### Master Script 9: Classification of Discharge Outcomes ####
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
library(caret)
library(kernlab)
library(rlist)

source('./functions/load_patient_clinical_data.R')
source('./functions/update_clinicalVariableList.R')
source('./functions/get_motion_features.R')
source('./functions/lol_project_motion_features.R')
source("./functions/cross_val_splits.R")
source("./functions/load_tf_patient_covariates.R")
source('./functions/cv_lol_project_motion_features.R')
source('./functions/prepare_training_covariates.R')
source('./functions/prepare_testing_covariates.R')
source('./functions/train_GLM.R')
source('./functions/train_SVM.R')
source('./functions/train_LDA.R')
source('./functions/train_kNN.R')
source('./functions/predict_GLM.R')
source('./functions/predict_SVM.R')
source('./functions/predict_LDA.R')
source('./functions/predict_kNN.R')

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
tod_tf_covariates <-load_tf_patient_covariates('../motion_feature_data/tf_patient_covariates.csv')
tfr_tf_covariates <-load_tf_patient_covariates('../tfr_motion_feature_data/tf_patient_covariates.csv')

# Split for k-fold cross validation for discharge predictions:
k <- 5
cvIdx <- cross_val_splits(patient_clinical_data, k)

self_val_GLM_predictionsDis <- vector(mode = "list", length = k)
test_val_GLM_predictionsDis <- vector(mode = "list", length = k)

self_val_SVM_predictionsDis <- vector(mode = "list", length = k)
test_val_SVM_predictionsDis <- vector(mode = "list", length = k)

self_val_LDA_predictionsDis <- vector(mode = "list", length = k)
test_val_LDA_predictionsDis <- vector(mode = "list", length = k)

self_val_kNN_predictionsDis <- vector(mode = "list", length = k)
test_val_kNN_predictionsDis <- vector(mode = "list", length = k)

for (i in 1:length(cvIdx)) {
  currTestIdx <- cvIdx[[as.character(i)]]
  currTrainIdx <- seq(nrow(patient_clinical_data))[-currTestIdx]
  
  # Convert numeric to logical indexing
  logicalTest <- rep(FALSE, nrow(patient_clinical_data))
  logicalTest[currTestIdx] <- TRUE
  logicalTrain <- !logicalTest

  # Perform LOL on training data
  r <- 3
  
  # - GOSE LOL:
  Y <- as.factor(patient_clinical_data$gose)
  
  tod_GOSE_LOL <-
    cv_lol_project_motion_features(tod_motion_features, Y, r, logicalTrain)
  tfr_GOSE_LOL <-
    cv_lol_project_motion_features(tfr_motion_features, Y, r, logicalTrain)
  
  # - Favorable Outcome LOL:
  Y <- as.factor(patient_clinical_data$favorable)
  
  tod_fav_LOL <-
    cv_lol_project_motion_features(tod_motion_features, Y, r, logicalTrain)
  tfr_fav_LOL <-
    cv_lol_project_motion_features(tfr_motion_features, Y, r, logicalTrain)
  
  # - Mortality Outcome LOL:
  Y <- as.factor(patient_clinical_data$death)
  
  tod_death_LOL <-
    cv_lol_project_motion_features(tod_motion_features, Y, r, logicalTrain)
  tfr_death_LOL <-
    cv_lol_project_motion_features(tfr_motion_features, Y, r, logicalTrain)
  
  # Prepare training covariates
  
  clinicalVars <- c("Sex","APACHE","CVA","ICH","SAH","BT","SDH.TBI")
  
  train_tod_GOSE_CV <-
    prepare_training_covariates(clinicalVars,tod_GOSE_LOL, tod_tf_covariates[currTrainIdx,],r)
  train_tod_fav_CV <-
    prepare_training_covariates(clinicalVars,tod_fav_LOL, tod_tf_covariates[currTrainIdx,],r)
  train_tod_death_CV <-
    prepare_training_covariates(clinicalVars,tod_death_LOL, tod_tf_covariates[currTrainIdx,],r)
  
  train_tfr_GOSE_CV <-
    prepare_training_covariates(clinicalVars,tfr_GOSE_LOL,tfr_tf_covariates[currTrainIdx,],r)
  train_tfr_fav_CV <-
    prepare_training_covariates(clinicalVars,tfr_fav_LOL,tfr_tf_covariates[currTrainIdx,],r)
  train_tfr_death_CV <-
    prepare_training_covariates(clinicalVars,tfr_death_LOL, tfr_tf_covariates[currTrainIdx,],r)
  
  # Prepare testing covariates
  test_tod_GOSE_CV <-
    prepare_testing_covariates(tod_motion_features,
                               tod_GOSE_LOL,
                               tod_tf_covariates[currTestIdx,],
                               currTestIdx,
                               r)
  test_tod_fav_CV <-
    prepare_testing_covariates(tod_motion_features,
                               tod_fav_LOL,
                               tod_tf_covariates[currTestIdx,],
                               currTestIdx,
                               r)
  test_tod_death_CV <-
    prepare_testing_covariates(tod_motion_features,
                               tod_death_LOL,
                               tod_tf_covariates[currTestIdx,],
                               currTestIdx,
                               r)
  
  test_tfr_GOSE_CV <-
    prepare_testing_covariates(tfr_motion_features,
                               tfr_GOSE_LOL,
                               tfr_tf_covariates[currTestIdx,],
                               currTestIdx,
                               r)
  test_tfr_fav_CV <-
    prepare_testing_covariates(tfr_motion_features,
                               tfr_fav_LOL,
                               tfr_tf_covariates[currTestIdx,],
                               currTestIdx,
                               r)
  test_tfr_death_CV <-
    prepare_testing_covariates(tfr_motion_features,
                               tfr_death_LOL,
                               tfr_tf_covariates[currTestIdx,],
                               currTestIdx,
                               r)
  
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
  
  # Train SVM on training data
  Y <- as.factor(patient_clinical_data$favorable)
  
  tod_GOSE_SVM <-
    train_SVM(train_tod_GOSE_CV, Y, currTrainIdx, 4)
  tfr_GOSE_SVM <-
    train_SVM(train_tfr_GOSE_CV, Y, currTrainIdx, 4)
  
  tod_fav_SVM <- train_SVM(train_tod_fav_CV, Y, currTrainIdx, 4)
  tfr_fav_SVM <- train_SVM(train_tfr_fav_CV, Y, currTrainIdx, 4)
  
  Y <- as.factor(patient_clinical_data$death)
  
  tod_death_SVM <-
    train_SVM(train_tod_death_CV, Y, currTrainIdx, 4)
  tfr_death_SVM <-
    train_SVM(train_tfr_death_CV, Y, currTrainIdx, 4)
  
  # Train LDA on training data
  
  Y <- as.factor(patient_clinical_data$favorable)
  
  tod_GOSE_LDA <-
    train_LDA(train_tod_GOSE_CV, Y, currTrainIdx)
  tfr_GOSE_LDA <-
    train_LDA(train_tfr_GOSE_CV, Y, currTrainIdx)
  
  tod_fav_LDA <- train_LDA(train_tod_fav_CV, Y, currTrainIdx)
  tfr_fav_LDA <- train_LDA(train_tfr_fav_CV, Y, currTrainIdx)
  
  Y <- as.factor(patient_clinical_data$death)
  
  tod_death_LDA <-
    train_LDA(train_tod_death_CV, Y, currTrainIdx)
  tfr_death_LDA <-
    train_LDA(train_tfr_death_CV, Y, currTrainIdx)
  
  # Train kNN on training data
  Y <- as.factor(patient_clinical_data$favorable)
  
  tod_GOSE_kNN <-
    train_kNN(train_tod_GOSE_CV, Y, currTrainIdx, 4)
  tfr_GOSE_kNN <-
    train_kNN(train_tfr_GOSE_CV, Y, currTrainIdx, 4)
  
  tod_fav_kNN <- train_kNN(train_tod_fav_CV, Y, currTrainIdx, 4)
  tfr_fav_kNN <- train_kNN(train_tfr_fav_CV, Y, currTrainIdx, 4)
  
  Y <- as.factor(patient_clinical_data$death)
  
  tod_death_kNN <-
    train_kNN(train_tod_death_CV, Y, currTrainIdx, 4)
  tfr_death_kNN <-
    train_kNN(train_tfr_death_CV, Y, currTrainIdx, 4)
  
  # Self-validate the GLM models on training data
  
  Y <- as.factor(patient_clinical_data$favorable)
  
  tod_GOSE_self_val <-
    predict_GLM(tod_GOSE_GLM, train_tod_GOSE_CV, Y, currTrainIdx)
  tfr_GOSE_self_val <-
    predict_GLM(tfr_GOSE_GLM, train_tfr_GOSE_CV, Y, currTrainIdx)
  
  tod_fav_self_val <-
    predict_GLM(tod_fav_GLM, train_tod_fav_CV, Y, currTrainIdx)
  tfr_fav_self_val <-
    predict_GLM(tfr_fav_GLM, train_tfr_fav_CV, Y, currTrainIdx)
  
  Y <- as.factor(patient_clinical_data$death)
  
  tod_death_self_val <-
    predict_GLM(tod_death_GLM, train_tod_death_CV, Y, currTrainIdx)
  tfr_death_self_val <-
    predict_GLM(tfr_death_GLM, train_tfr_death_CV, Y, currTrainIdx)
  
  temp_GLM_self_val_list <-
    list(
      tod_GOSE_self_val,
      tfr_GOSE_self_val,
      tod_fav_self_val,
      tfr_fav_self_val,
      tod_death_self_val,
      tfr_death_self_val
    )
  names(temp_GLM_self_val_list) <-
    c(
      "tod_GOSE_self_val",
      "tfr_GOSE_self_val",
      "tod_fav_self_val",
      "tfr_fav_self_val",
      "tod_death_self_val",
      "tfr_death_self_val"
    )
  self_val_GLM_predictionsDis[[i]] <- temp_GLM_self_val_list
  
  # Self-validate the SVM models on training data
  
  Y <- as.factor(patient_clinical_data$favorable)
  
  tod_GOSE_self_val <-
    predict_SVM(tod_GOSE_SVM, train_tod_GOSE_CV, Y, currTrainIdx)
  tfr_GOSE_self_val <-
    predict_SVM(tfr_GOSE_SVM, train_tfr_GOSE_CV, Y, currTrainIdx)
  
  tod_fav_self_val <-
    predict_SVM(tod_fav_SVM, train_tod_fav_CV, Y, currTrainIdx)
  tfr_fav_self_val <-
    predict_SVM(tfr_fav_SVM, train_tfr_fav_CV, Y, currTrainIdx)
  
  Y <- as.factor(patient_clinical_data$death)
  
  tod_death_self_val <-
    predict_SVM(tod_death_SVM, train_tod_death_CV, Y, currTrainIdx)
  tfr_death_self_val <-
    predict_SVM(tfr_death_SVM, train_tfr_death_CV, Y, currTrainIdx)
  
  temp_SVM_self_val_list <-
    list(
      tod_GOSE_self_val,
      tfr_GOSE_self_val,
      tod_fav_self_val,
      tfr_fav_self_val,
      tod_death_self_val,
      tfr_death_self_val
    )
  names(temp_SVM_self_val_list) <-
    c(
      "tod_GOSE_self_val",
      "tfr_GOSE_self_val",
      "tod_fav_self_val",
      "tfr_fav_self_val",
      "tod_death_self_val",
      "tfr_death_self_val"
    )
  self_val_SVM_predictionsDis[[i]] <- temp_SVM_self_val_list
  
  # Self-validate the LDA models on training data
  
  Y <- as.factor(patient_clinical_data$favorable)
  
  tod_GOSE_self_val <-
    predict_LDA(tod_GOSE_LDA, train_tod_GOSE_CV, Y, currTrainIdx)
  tfr_GOSE_self_val <-
    predict_LDA(tfr_GOSE_LDA, train_tfr_GOSE_CV, Y, currTrainIdx)
  
  tod_fav_self_val <-
    predict_LDA(tod_fav_LDA, train_tod_fav_CV, Y, currTrainIdx)
  tfr_fav_self_val <-
    predict_LDA(tfr_fav_LDA, train_tfr_fav_CV, Y, currTrainIdx)
  
  Y <- as.factor(patient_clinical_data$death)
  
  tod_death_self_val <-
    predict_LDA(tod_death_LDA, train_tod_death_CV, Y, currTrainIdx)
  tfr_death_self_val <-
    predict_LDA(tfr_death_LDA, train_tfr_death_CV, Y, currTrainIdx)
  
  temp_LDA_self_val_list <-
    list(
      tod_GOSE_self_val,
      tfr_GOSE_self_val,
      tod_fav_self_val,
      tfr_fav_self_val,
      tod_death_self_val,
      tfr_death_self_val
    )
  names(temp_LDA_self_val_list) <-
    c(
      "tod_GOSE_self_val",
      "tfr_GOSE_self_val",
      "tod_fav_self_val",
      "tfr_fav_self_val",
      "tod_death_self_val",
      "tfr_death_self_val"
    )
  self_val_LDA_predictionsDis[[i]] <- temp_LDA_self_val_list
  
  # Self-validate the kNN models on training data
  
  Y <- as.factor(patient_clinical_data$favorable)
  
  tod_GOSE_self_val <-
    predict_kNN(tod_GOSE_kNN, train_tod_GOSE_CV, Y, currTrainIdx)
  tfr_GOSE_self_val <-
    predict_kNN(tfr_GOSE_kNN, train_tfr_GOSE_CV, Y, currTrainIdx)
  
  tod_fav_self_val <-
    predict_kNN(tod_fav_kNN, train_tod_fav_CV, Y, currTrainIdx)
  tfr_fav_self_val <-
    predict_kNN(tfr_fav_kNN, train_tfr_fav_CV, Y, currTrainIdx)
  
  Y <- as.factor(patient_clinical_data$death)
  
  tod_death_self_val <-
    predict_kNN(tod_death_kNN, train_tod_death_CV, Y, currTrainIdx)
  tfr_death_self_val <-
    predict_kNN(tfr_death_kNN, train_tfr_death_CV, Y, currTrainIdx)
  
  temp_kNN_self_val_list <-
    list(
      tod_GOSE_self_val,
      tfr_GOSE_self_val,
      tod_fav_self_val,
      tfr_fav_self_val,
      tod_death_self_val,
      tfr_death_self_val
    )
  names(temp_kNN_self_val_list) <-
    c(
      "tod_GOSE_self_val",
      "tfr_GOSE_self_val",
      "tod_fav_self_val",
      "tfr_fav_self_val",
      "tod_death_self_val",
      "tfr_death_self_val"
    )
  self_val_kNN_predictionsDis[[i]] <- temp_kNN_self_val_list
  
  # Predict outcomes for testing data using GLM
  
  Y <- as.factor(patient_clinical_data$favorable)
  
  tod_GOSE_test_val <-
    predict_GLM(tod_GOSE_GLM, test_tod_GOSE_CV, Y, currTestIdx)
  tfr_GOSE_test_val <-
    predict_GLM(tfr_GOSE_GLM, test_tfr_GOSE_CV, Y, currTestIdx)
  
  tod_fav_test_val <-
    predict_GLM(tod_fav_GLM, test_tod_fav_CV, Y, currTestIdx)
  tfr_fav_test_val <-
    predict_GLM(tfr_fav_GLM, test_tfr_fav_CV, Y, currTestIdx)
  
  Y <- as.factor(patient_clinical_data$death)
  
  tod_death_test_val <-
    predict_GLM(tod_death_GLM, test_tod_death_CV, Y, currTestIdx)
  tfr_death_test_val <-
    predict_GLM(tfr_death_GLM, test_tfr_death_CV, Y, currTestIdx)
  
  temp_GLM_test_val_list <-
    list(
      tod_GOSE_test_val,
      tfr_GOSE_test_val,
      tod_fav_test_val,
      tfr_fav_test_val,
      tod_death_test_val,
      tfr_death_test_val
    )
  names(temp_GLM_test_val_list) <-
    c(
      "tod_GOSE_test_val",
      "tfr_GOSE_test_val",
      "tod_fav_test_val",
      "tfr_fav_test_val",
      "tod_death_test_val",
      "tfr_death_test_val"
    )
  test_val_GLM_predictionsDis[[i]] <- temp_GLM_test_val_list
  
  # Predict outcomes for testing data using SVM
  
  Y <- as.factor(patient_clinical_data$favorable)
  
  tod_GOSE_test_val <-
    predict_SVM(tod_GOSE_SVM, test_tod_GOSE_CV, Y, currTestIdx)
  tfr_GOSE_test_val <-
    predict_SVM(tfr_GOSE_SVM, test_tfr_GOSE_CV, Y, currTestIdx)
  
  tod_fav_test_val <-
    predict_SVM(tod_fav_SVM, test_tod_fav_CV, Y, currTestIdx)
  tfr_fav_test_val <-
    predict_SVM(tfr_fav_SVM, test_tfr_fav_CV, Y, currTestIdx)
  
  Y <- as.factor(patient_clinical_data$death)
  
  tod_death_test_val <-
    predict_SVM(tod_death_SVM, test_tod_death_CV, Y, currTestIdx)
  tfr_death_test_val <-
    predict_SVM(tfr_death_SVM, test_tfr_death_CV, Y, currTestIdx)
  
  temp_SVM_test_val_list <-
    list(
      tod_GOSE_test_val,
      tfr_GOSE_test_val,
      tod_fav_test_val,
      tfr_fav_test_val,
      tod_death_test_val,
      tfr_death_test_val
    )
  names(temp_SVM_test_val_list) <-
    c(
      "tod_GOSE_test_val",
      "tfr_GOSE_test_val",
      "tod_fav_test_val",
      "tfr_fav_test_val",
      "tod_death_test_val",
      "tfr_death_test_val"
    )
  test_val_SVM_predictionsDis[[i]] <- temp_SVM_test_val_list
  
  # Predict outcomes for testing data using LDA
  
  Y <- as.factor(patient_clinical_data$favorable)
  
  tod_GOSE_test_val <-
    predict_LDA(tod_GOSE_LDA, test_tod_GOSE_CV, Y, currTestIdx)
  tfr_GOSE_test_val <-
    predict_LDA(tfr_GOSE_LDA, test_tfr_GOSE_CV, Y, currTestIdx)
  
  tod_fav_test_val <-
    predict_LDA(tod_fav_LDA, test_tod_fav_CV, Y, currTestIdx)
  tfr_fav_test_val <-
    predict_LDA(tfr_fav_LDA, test_tfr_fav_CV, Y, currTestIdx)
  
  Y <- as.factor(patient_clinical_data$death)
  
  tod_death_test_val <-
    predict_LDA(tod_death_LDA, test_tod_death_CV, Y, currTestIdx)
  tfr_death_test_val <-
    predict_LDA(tfr_death_LDA, test_tfr_death_CV, Y, currTestIdx)
  
  temp_LDA_test_val_list <-
    list(
      tod_GOSE_test_val,
      tfr_GOSE_test_val,
      tod_fav_test_val,
      tfr_fav_test_val,
      tod_death_test_val,
      tfr_death_test_val
    )
  names(temp_LDA_test_val_list) <-
    c(
      "tod_GOSE_test_val",
      "tfr_GOSE_test_val",
      "tod_fav_test_val",
      "tfr_fav_test_val",
      "tod_death_test_val",
      "tfr_death_test_val"
    )
  test_val_LDA_predictionsDis[[i]] <- temp_LDA_test_val_list
  
  # Predict outcomes for testing data using kNN
  
  Y <- as.factor(patient_clinical_data$favorable)
  
  tod_GOSE_test_val <-
    predict_kNN(tod_GOSE_kNN, test_tod_GOSE_CV, Y, currTestIdx)
  tfr_GOSE_test_val <-
    predict_kNN(tfr_GOSE_kNN, test_tfr_GOSE_CV, Y, currTestIdx)
  
  tod_fav_test_val <-
    predict_kNN(tod_fav_kNN, test_tod_fav_CV, Y, currTestIdx)
  tfr_fav_test_val <-
    predict_kNN(tfr_fav_kNN, test_tfr_fav_CV, Y, currTestIdx)
  
  Y <- as.factor(patient_clinical_data$death)
  
  tod_death_test_val <-
    predict_kNN(tod_death_kNN, test_tod_death_CV, Y, currTestIdx)
  tfr_death_test_val <-
    predict_kNN(tfr_death_kNN, test_tfr_death_CV, Y, currTestIdx)
  
  temp_kNN_test_val_list <-
    list(
      tod_GOSE_test_val,
      tfr_GOSE_test_val,
      tod_fav_test_val,
      tfr_fav_test_val,
      tod_death_test_val,
      tfr_death_test_val
    )
  names(temp_kNN_test_val_list) <-
    c(
      "tod_GOSE_test_val",
      "tfr_GOSE_test_val",
      "tod_fav_test_val",
      "tfr_fav_test_val",
      "tod_death_test_val",
      "tfr_death_test_val"
    )
  test_val_kNN_predictionsDis[[i]] <- temp_kNN_test_val_list
}

rm(list = setdiff(
  ls(),
  c(
    "self_val_GLM_predictionsDis",
    "test_val_GLM_predictionsDis",
    "self_val_SVM_predictionsDis",
    "test_val_SVM_predictionsDis",
    "self_val_LDA_predictionsDis",
    "test_val_LDA_predictionsDis",
    "self_val_kNN_predictionsDis",
    "test_val_kNN_predictionsDis",
    "patient_clinical_data",
    "clinicalVariableList",
    "tod_motion_features",
    "tfr_motion_features"
  )
))