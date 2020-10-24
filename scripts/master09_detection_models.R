#### Master Script 9: Train detection models ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

# Add library path for R packages and set working directory:
.libPaths(c("~/Rpackages" , .libPaths()))
setwd('~/work/nims/scripts')

# Load necessary packages
library(caret)
library(tidyverse)
library(readxl)
library(caret)
library(lolR)
library(UBL)
library(glmnet)
library(keras)
library(caTools)
library(kernlab)
library(nnet)
library(parallel)
library(fastAdaboost)
library(xgboost)
library(kknn)
library(cutpointr)
library(precrec)
library(reshape2)
library(MLmetrics)
library(gridExtra)
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

# Load clinical patient data:
source('./functions/load_patient_clinical_data.R')
patient_clinical_data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv') %>% arrange(AccelPatientNo_) %>% mutate(ptIdx = 1:nrow(.))

# Directories of imputations
impDirs <- list.files('~/scratch/all_motion_feature_data/formatted_matrices',include.dirs = TRUE, full.names = TRUE)

# Load detection labels and partitions
source('./functions/multinomial_to_ordinal.R')
load('~/scratch/all_motion_feature_data/gcs_labels/detection_labels.RData')
load('~/scratch/all_motion_feature_data/gcs_labels/detection_partitions.RData')

# Define classifiers to train
classifier_choice <- c("kknn","adaboost","glmnet", "parRF", "svmRadialWeights","lda")

# Define tuning grids for FIRST TUNING RUN:

tune.grid.adaboost <- expand.grid(method = "Adaboost.M1", nIter = c(10,30,100,300,1000))
tune.grid.glmnet <- expand.grid(alpha = c(0,.5,1),lambda = seq(0.0001, 1, length = 50))
tune.grid.parRF <- expand.grid(mtry = (1:15))
tune.grid.svmRadialWeights <- expand.grid(C = c(1,3,5,10,20), 
                                          Weight = c(0.1,0.5,1,2,3,5,10),
                                          sigma = c(0.0005,0.001,0.005,0.01,0.05))
tune.grid.kknn <- expand.grid(kmax = seq(3,15,by = 2),distance = 2, kernel = "optimal")

nTunes <- max(as.numeric(lapply(
  list(
    tune.grid.adaboost,
    tune.grid.glmnet,
    tune.grid.parRF,
    tune.grid.svmRadialWeights
  ),
  nrow
)))

# Define inner fold count for parameter tuning and validation
inner_fold_count <- 5

# Set seeds for tuning
source('./functions/setSeeds.R')
seed.list <- setSeeds(method = "cv",numbers=inner_fold_count,repeats=1,tunes = nTunes,seed=2020)

# Set destination directory
dir.create(path = "~/scratch/all_motion_feature_data/detection_results",showWarnings = FALSE)  

# Load ML tuning function
source("./functions/tuneDetectionMLModels.R")

for (i in rev(1:length(impDirs))){
  print(paste("Imputation no.",i,"out of",length(impDirs),"started."))
  currImpDir <- impDirs[i]
  detection_folders <- list.files(path = currImpDir,pattern = 'detection_*',include.dirs = TRUE, full.names = TRUE)
  
  # Create directory for imputation
  dir.create(path = file.path("~/scratch/all_motion_feature_data/detection_results",paste0("imp",i)),showWarnings = FALSE)  
  
  for (j in 1:length(detection_folders)){
    
    gc()
    
    curr_window_size <- as.numeric(sub(".*window_", "", detection_folders[j]))
    curr_window_idx <- which(det_parameters$obs_windows == curr_window_size)
    
    curr_label_set <- det_gcs_labels[[curr_window_idx]]
    
    # Create directory for current detection paradigm
    dir.create(paste0('~/scratch/all_motion_feature_data/detection_results/imp',i,'/detection_window_',curr_window_size),showWarnings = FALSE)
    
    print(paste("Detection folder no.",j,"out of",length(detection_folders),"started."))
    
    tryCatch({
      curr_motor_train_smote <- readRDS(file.path(detection_folders[j],'motor_train_smote.rds'))
      curr_motor_ordinal_labels <- multinomial_to_ordinal(curr_motor_train_smote$labels)
      
      cat("ML ON MOTOR, WINDOW SIZE", curr_window_size,"HOURS STARTED.","\n")
      
      # Create directory for motor GCS detection results
      dir.create(paste0('~/scratch/all_motion_feature_data/detection_results/imp',i,'/detection_window_',curr_window_size,"/motor"),showWarnings = FALSE)
      
      tuneDetectionMLModels(
        curr_SMOTE = curr_motor_train_smote,
        inner_fold_count = inner_fold_count,
        classifier_choice = classifier_choice,
        save.path = paste0('~/scratch/all_motion_feature_data/detection_results/imp',i,'/detection_window_',curr_window_size,'/motor'),
        seed.list = seed.list,
        ordinalLabels = curr_motor_ordinal_labels
      )
      
    }, error=function(e){cat("ERROR ON MOTOR, WINDOW SIZE", curr_window_size,"HOURS:",conditionMessage(e), "\n")})
    
    tryCatch({
      curr_eye_train_smote <- readRDS(file.path(detection_folders[j],'eye_train_smote.rds'))
      curr_eye_ordinal_labels <- multinomial_to_ordinal(curr_eye_train_smote$labels)
      
      cat("ML ON EYE, WINDOW SIZE", curr_window_size,"HOURS STARTED.","\n")
      
      # Create directory for eye GCS detection results
      dir.create(paste0('~/scratch/all_motion_feature_data/detection_results/imp',i,'/detection_window_',curr_window_size,"/eye"),showWarnings = FALSE)
      
      tuneDetectionMLModels(
        curr_SMOTE = curr_eye_train_smote,
        inner_fold_count = inner_fold_count,
        classifier_choice = classifier_choice,
        save.path = paste0('~/scratch/all_motion_feature_data/detection_results/imp',i,'/detection_window_',curr_window_size,'/eye'),
        seed.list = seed.list,
        ordinalLabels = curr_eye_ordinal_labels
      )
      
    }, error=function(e){cat("ERROR ON EYE, WINDOW SIZE", curr_window_size,"HOURS:",conditionMessage(e), "\n")})
    
    print(paste("Detection folder no.",j,"out of",length(detection_folders),"completed."))
  }
  print(paste("Imputation no.",i,"out of",length(impDirs),"completed."))
}

# Produce model predictions on training and testing sets
rm(list = ls())
gc()

# Directories of imputations
impDirs <- list.files('~/scratch/all_motion_feature_data/detection_results',include.dirs = TRUE, full.names = TRUE, pattern = "imp*")

source("./functions/test_detection_model.R")
dir.create('~/scratch/all_motion_feature_data/detection_results/test_results',showWarnings = FALSE)
dir.create('~/scratch/all_motion_feature_data/detection_results/train_results',showWarnings = FALSE)
classifier_choice <- c("kknn","adaboost","glmnet", "parRF", "svmRadialWeights","lda")

# Load detection model labels and partitions
load('~/scratch/all_motion_feature_data/gcs_labels/detection_labels.RData')
load('~/scratch/all_motion_feature_data/gcs_labels/detection_partitions.RData')

for (i in 1:length(det_parameters$obs_windows)){
  
  curr_obs_window <- det_parameters$obs_windows[i]
  
  curr_labels <- det_gcs_labels[[i]]
  
  motor_indices <- which(!is.na(curr_labels$Best.Motor.Response))
  eye_indices <- which(!is.na(curr_labels$Eye.Opening))
  
  curr_det_motor_train_idx <- det_motor_train_idx[[i]]
  curr_det_eye_train_idx <- det_eye_train_idx[[i]]
  
  curr_det_motor_test_idx <- setdiff(motor_indices,curr_det_motor_train_idx)
  curr_det_eye_test_idx <- setdiff(eye_indices,curr_det_eye_train_idx)
  
  curr_det_motor_train_labels <- curr_labels$Best.Motor.Response[curr_det_motor_train_idx]
  curr_det_eye_train_labels <-  curr_labels$Eye.Opening[curr_det_eye_train_idx]
  
  curr_det_motor_test_labels <- curr_labels$Best.Motor.Response[curr_det_motor_test_idx]
  curr_det_eye_test_labels <-  curr_labels$Eye.Opening[curr_det_eye_test_idx]
  
  dir.create(paste0('~/scratch/all_motion_feature_data/detection_results/test_results/detection_window_',curr_obs_window),showWarnings = FALSE)
  dir.create(paste0('~/scratch/all_motion_feature_data/detection_results/train_results/detection_window_',curr_obs_window),showWarnings = FALSE)
  
  dir.create(paste0('~/scratch/all_motion_feature_data/detection_results/test_results/detection_window_',curr_obs_window,'/motor'),showWarnings = FALSE)
  dir.create(paste0('~/scratch/all_motion_feature_data/detection_results/test_results/detection_window_',curr_obs_window,'/eye'),showWarnings = FALSE)
  
  dir.create(paste0('~/scratch/all_motion_feature_data/detection_results/train_results/detection_window_',curr_obs_window,'/motor'),showWarnings = FALSE)
  dir.create(paste0('~/scratch/all_motion_feature_data/detection_results/train_results/detection_window_',curr_obs_window,'/eye'),showWarnings = FALSE)
  
  for (modelType in classifier_choice){
    
    motor_test_results <- as.data.frame(matrix(ncol = 3,nrow = 0))
    motor_train_results <- as.data.frame(matrix(ncol = 3,nrow = 0))
    eye_test_results <- as.data.frame(matrix(ncol = 3,nrow = 0))
    eye_train_results <- as.data.frame(matrix(ncol = 3,nrow = 0))
    
    for (impNo in 1:length(impDirs)){
      tryCatch({
        curr_motor_mdl <- readRDS(paste0("~/scratch/all_motion_feature_data/detection_results/imp",impNo,"/detection_window_",curr_obs_window,"/motor/",modelType,".rds"))
        
        curr_motor_test_lol <- as.data.frame(readRDS(paste0("~/scratch/all_motion_feature_data/formatted_matrices/imp",impNo,"/detection_window_",curr_obs_window,"/","motor_test_matrix_lol.rds")))
        curr_motor_train_smote <- as.data.frame(readRDS(paste0("~/scratch/all_motion_feature_data/formatted_matrices/imp",impNo,"/detection_window_",curr_obs_window,"/","motor_train_smote.rds")))
        curr_motor_train_smote <- curr_motor_train_smote[,1:(length(curr_motor_train_smote)-1)]
        
        curr_motor_train_predictions <- test_detection_model(modelOutput = curr_motor_mdl, predictor_matrix = curr_motor_train_smote, trueLabels = curr_det_motor_train_labels)
        curr_motor_train_predictions$imp <- impNo
        
        curr_motor_test_predictions <- test_detection_model(modelOutput = curr_motor_mdl, predictor_matrix = curr_motor_test_lol, trueLabels = curr_det_motor_test_labels)
        curr_motor_test_predictions$imp <- impNo
        
        motor_train_results <- rbind(motor_train_results,curr_motor_train_predictions)
        motor_test_results <- rbind(motor_test_results,curr_motor_test_predictions)
        
      }, error=function(e){cat("ERROR ON MOTOR, WINDOW SIZE", curr_obs_window,"HOURS, IMPUTATION NO",impNo,conditionMessage(e), "\n")})
      
      tryCatch({
        curr_eye_mdl <- readRDS(paste0("~/scratch/all_motion_feature_data/detection_results/imp",impNo,"/detection_window_",curr_obs_window,"/eye/",modelType,".rds"))
        
        curr_eye_test_lol <- as.data.frame(readRDS(paste0("~/scratch/all_motion_feature_data/formatted_matrices/imp",impNo,"/detection_window_",curr_obs_window,"/","eye_test_matrix_lol.rds")))
        curr_eye_train_smote <- as.data.frame(readRDS(paste0("~/scratch/all_motion_feature_data/formatted_matrices/imp",impNo,"/detection_window_",curr_obs_window,"/","eye_train_smote.rds")))
        curr_eye_train_smote <- curr_eye_train_smote[,1:(length(curr_eye_train_smote)-1)]
        
        curr_eye_train_predictions <- test_detection_model(modelOutput = curr_eye_mdl, predictor_matrix = curr_eye_train_smote, trueLabels = curr_det_eye_train_labels)
        curr_eye_train_predictions$imp <- impNo
        
        curr_eye_test_predictions <- test_detection_model(modelOutput = curr_eye_mdl, predictor_matrix = curr_eye_test_lol, trueLabels = curr_det_eye_test_labels)
        curr_eye_test_predictions$imp <- impNo
        
        eye_train_results <- rbind(eye_train_results,curr_eye_train_predictions)
        eye_test_results <- rbind(eye_test_results,curr_eye_test_predictions)
      }, error=function(e){cat("ERROR ON EYE, WINDOW SIZE", curr_obs_window,"HOURS, IMPUTATION NO",impNo,conditionMessage(e), "\n")})
    }
    
    saveRDS(motor_test_results,paste0('~/scratch/all_motion_feature_data/detection_results/test_results/detection_window_',curr_obs_window,'/motor/',modelType,'_results.rds'))
    saveRDS(motor_train_results,paste0('~/scratch/all_motion_feature_data/detection_results/train_results/detection_window_',curr_obs_window,'/motor/',modelType,'_results.rds'))
    
    saveRDS(eye_test_results,paste0('~/scratch/all_motion_feature_data/detection_results/test_results/detection_window_',curr_obs_window,'/eye/',modelType,'_results.rds'))
    saveRDS(eye_train_results,paste0('~/scratch/all_motion_feature_data/detection_results/train_results/detection_window_',curr_obs_window,'/eye/',modelType,'_results.rds'))
  }
  
}

# Produce confusion matrices for model results:
rm(list = ls())
gc()

# Directories of imputations
impDirs <- list.files('~/scratch/all_motion_feature_data/detection_results',include.dirs = TRUE, full.names = TRUE, pattern = "imp*")

dir.create('~/scratch/all_motion_feature_data/detection_results/test_results',showWarnings = FALSE)
dir.create('~/scratch/all_motion_feature_data/detection_results/train_results',showWarnings = FALSE)
classifier_choice <- c("kknn","adaboost","glmnet", "parRF", "svmRadialWeights","lda")

for (i in 1:length(det_parameters$obs_windows)){
  
  curr_obs_window <- det_parameters$obs_windows[i]
  
  for (modelType in classifier_choice){
    
    for (impNo in 1:length(impDirs)){
      tryCatch({
        curr_motor_mdl <- readRDS(paste0("~/scratch/all_motion_feature_data/detection_results/imp",impNo,"/detection_window_",curr_obs_window,"/motor/",modelType,".rds"))
        
        curr_motor_test_lol <- as.data.frame(readRDS(paste0("~/scratch/all_motion_feature_data/formatted_matrices/imp",impNo,"/detection_window_",curr_obs_window,"/","motor_test_matrix_lol.rds")))
        curr_motor_train_smote <- as.data.frame(readRDS(paste0("~/scratch/all_motion_feature_data/formatted_matrices/imp",impNo,"/detection_window_",curr_obs_window,"/","motor_train_smote.rds")))
        curr_motor_train_smote <- curr_motor_train_smote[,1:(length(curr_motor_train_smote)-1)]
        
        curr_motor_train_predictions <- test_detection_model(modelOutput = curr_motor_mdl, predictor_matrix = curr_motor_train_smote, trueLabels = curr_det_motor_train_labels)
        curr_motor_train_predictions$imp <- impNo
        
        curr_motor_test_predictions <- test_detection_model(modelOutput = curr_motor_mdl, predictor_matrix = curr_motor_test_lol, trueLabels = curr_det_motor_test_labels)
        curr_motor_test_predictions$imp <- impNo
        
        motor_train_results <- rbind(motor_train_results,curr_motor_train_predictions)
        motor_test_results <- rbind(motor_test_results,curr_motor_test_predictions)
        
      }, error=function(e){cat("ERROR ON MOTOR, WINDOW SIZE", curr_obs_window,"HOURS, IMPUTATION NO",impNo,conditionMessage(e), "\n")})
      
      tryCatch({
        curr_eye_mdl <- readRDS(paste0("~/scratch/all_motion_feature_data/detection_results/imp",impNo,"/detection_window_",curr_obs_window,"/eye/",modelType,".rds"))
        
        curr_eye_test_lol <- as.data.frame(readRDS(paste0("~/scratch/all_motion_feature_data/formatted_matrices/imp",impNo,"/detection_window_",curr_obs_window,"/","eye_test_matrix_lol.rds")))
        curr_eye_train_smote <- as.data.frame(readRDS(paste0("~/scratch/all_motion_feature_data/formatted_matrices/imp",impNo,"/detection_window_",curr_obs_window,"/","eye_train_smote.rds")))
        curr_eye_train_smote <- curr_eye_train_smote[,1:(length(curr_eye_train_smote)-1)]
        
        curr_eye_train_predictions <- test_detection_model(modelOutput = curr_eye_mdl, predictor_matrix = curr_eye_train_smote, trueLabels = curr_det_eye_train_labels)
        curr_eye_train_predictions$imp <- impNo
        
        curr_eye_test_predictions <- test_detection_model(modelOutput = curr_eye_mdl, predictor_matrix = curr_eye_test_lol, trueLabels = curr_det_eye_test_labels)
        curr_eye_test_predictions$imp <- impNo
        
        eye_train_results <- rbind(eye_train_results,curr_eye_train_predictions)
        eye_test_results <- rbind(eye_test_results,curr_eye_test_predictions)
      }, error=function(e){cat("ERROR ON EYE, WINDOW SIZE", curr_obs_window,"HOURS, IMPUTATION NO",impNo,conditionMessage(e), "\n")})
    }
    
    saveRDS(motor_test_results,paste0('~/scratch/all_motion_feature_data/detection_results/test_results/detection_window_',curr_obs_window,'/motor/',modelType,'_results.rds'))
    saveRDS(motor_train_results,paste0('~/scratch/all_motion_feature_data/detection_results/train_results/detection_window_',curr_obs_window,'/motor/',modelType,'_results.rds'))
    
    saveRDS(eye_test_results,paste0('~/scratch/all_motion_feature_data/detection_results/test_results/detection_window_',curr_obs_window,'/eye/',modelType,'_results.rds'))
    saveRDS(eye_train_results,paste0('~/scratch/all_motion_feature_data/detection_results/train_results/detection_window_',curr_obs_window,'/eye/',modelType,'_results.rds'))
  }
  
}


Prediction: 24 hours
Lead Tme: 1 hour
Obs: test different possibilites (0.5 --- 1)







trialResults <- readRDS("../all_motion_feature_data/detection_results/test_results/detection_window_1/motor/glmnet_results.rds")
confusionMatrix(factor(trialResults$predLabels),factor(trialResults$trueLabels))
