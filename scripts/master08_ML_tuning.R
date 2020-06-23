#### Master Script 8: Tuning ML Classification of Temporal Segments ####
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
library(kernlab)
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

path.save <- "../all_motion_feature_data/ML_results"
dir.create(path.save,showWarnings = F)

### Load clinical metadata ###
source('./functions/load_patient_clinical_data.R')
patient_clinical_data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv') %>% add_column(ptIdx = 1:nrow(.),.after = 2)
gose_thresh <- 5 # inclusive threshold
mrs_thresh <- 3 # inclusive threshold
patient_clinical_data <- mutate(patient_clinical_data,
                                fav_mort_dis = factor(DiedDuringThisHospitalStay_ == "0",labels = c("Unfav","Fav")),
                                fav_mort_12m = factor(Death12Months == "0",labels = c("Unfav","Fav")),
                                fav_GOSE_dis = factor(GOS_EDischarge >= gose_thresh,labels = c("Unfav","Fav")),
                                fav_GOSE_12m = factor(GOS_E12Months >= gose_thresh,labels = c("Unfav","Fav")),
                                fav_mRS_dis = factor(mRSDischarge <= mrs_thresh,labels = c("Unfav","Fav")),
                                fav_mRS_12m = factor(mRS12Months <= mrs_thresh,labels = c("Unfav","Fav")))

### Create and save outer folds for ML based on available outcome labels ###
labels.temp <- expand.grid(label=c("fav_mort","fav_GOSE","fav_mRS"),
                           temp = c("dis","12m"))
outer_fold_count <- 5
inner_fold_count <- 5

# source('./functions/createOuterCVFolds.R')
# outerFolds <- createOuterCVFolds(patient_clinical_data,outer_fold_count)
# saveRDS(outerFolds,'../all_motion_feature_data/outerFolds.rds')
# rm(outerFolds)
outerFolds <- readRDS('../all_motion_feature_data/outerFolds.rds')

### Write LOL-embedded motion features for training and testing ###
# Choose sensors and motion features to test under a given temporal segment window:

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

# Define models to tune and train
classifier_choice <- c("adaboost", "avNNet", "DeepNN", "glmnet", "parRF", "svmRadialWeights","lda")

# Set tuning grids for FIRST TUNING RUN:

tune.grid.adaboost <- expand.grid(method = "Adaboost.M1", nIter = c(10,30,100,300,1000))
tune.grid.avNNet <- expand.grid(size = c(2,3,5,10,20), 
                                decay = c(0.5,1,3,4,5,10), 
                                bag = c(TRUE,FALSE))
tune.grid.DeepNN <- expand.grid(complexity.multiplier = c(0.3,0.5,0.8,1,1.2),
                                activation.layer_1 = c('tanh', 'relu'),
                                activation.layer_2 = c('tanh', 'relu'),
                                activation.layer_3 = c('sigmoid', 'tanh', 'relu'),
                                activation.layer_4 = 'tanh',
                                activation.layer_5 = 'relu',
                                activation.layer_6 = 'relu',
                                rate.dropout_1 = c(0.3,0.4),
                                rate.dropout_2 = c(0.2,0.3),
                                rate.dropout_3 = 0.2,
                                rate.dropout_4 = 0.2,
                                rate.dropout_5 = 0.2,
                                rate.dropout_6 = 0.2,
                                Epochs = c(30,50))
tune.grid.glmnet <- expand.grid(alpha = c(0,.5,1),lambda = seq(0.0001, 1, length = 50))
tune.grid.parRF <- expand.grid(mtry = (1:15))
tune.grid.svmRadialWeights <- expand.grid(C = c(1,3,5,10,20), 
                                         Weight = c(0.1,0.5,1,2,3,5,10),
                                         sigma = c(0.0005,0.001,0.005,0.01,0.05))

# Set seeds for tuning caret models
nTunes <- max(as.numeric(lapply(
  list(
    tune.grid.adaboost,
    tune.grid.avNNet,
    tune.grid.glmnet,
    tune.grid.parRF,
    tune.grid.svmRadialWeights
  ),
  nrow
)))

source('./functions/setSeeds.R')
seed.list <- setSeeds(method = "cv",numbers=inner_fold_count,repeats=1,tunes = nTunes,seed=2020)

# define the train control for all models
train.control <- trainControl(method="repeatedcv",
                              number=inner_fold_count,
                              repeats=1,
                              summaryFunction = twoClassSummary, 
                              classProbs = TRUE,
                              seeds = seed.list,  
                              savePredictions = "all",
                              returnResamp = "all")

# Prepare motion features for model building

# Run the model building function
source('./functions/tuneMachineLearningModels.R')

iter <- 1 # iteration 1 for caret models
deep.iter <- 1 # iteration 1 for deep learning models

tuneMachineLearningModels(Iter = iter, 
                          DeepIter = deep.iter, 
                          classifier_choice = classifier_choice, 
                          seed.list = seed.list,
                          path.D = path.save, 
                          labelsList = labels.temp, 
                          fullData.imput.in = fullData.imput1)



curr_data_table,outcome,classifier_choice,patient_clinical_data,mf_col_idx,outer_fold_count,inner_fold_count


classifier_choice<-c("lda","glmnet","avNNet","parRF")

outcome <- "fav_GOSE_dis"



tune.grid.avNNet <- expand.grid(size = c(2,3,5,10,20), 
                                decay = c(0.5,1,3,4,5,10), 
                                bag = c(TRUE,FALSE))
tune.grid.svmRdialWeights <- expand.grid(C = c(1,3,5,10,20), 
                                         Weight = c(0.1,0.5,1,2,3,5,10),
                                         sigma = c(0.0005,0.001,0.005,0.01,0.05))

tune.grid.adaboost <- expand.grid(method = "Adaboost.M1", nIter = c(10,30,100,300,1000))
tune.grid.DeepNN <- expand.grid(complexity.multiplier = c(0.3,0.5,0.8,1,1.2),
                                activation.layer_1 = c('tanh', 'relu'),
                                activation.layer_2 = c('tanh', 'relu'),
                                activation.layer_3 = c('sigmoid', 'tanh', 'relu'),
                                activation.layer_4 = 'tanh',
                                activation.layer_5 = 'relu',
                                activation.layer_6 = 'relu',
                                rate.dropout_1 = c(0.3,0.4),
                                rate.dropout_2 = c(0.2,0.3),
                                rate.dropout_3 = 0.2,
                                rate.dropout_4 = 0.2,
                                rate.dropout_5 = 0.2,
                                rate.dropout_6 = 0.2,
                                Epochs = c(30,50))

for (seg_windows in c(5,10,30,60,180)){
  
  print(paste('Segment size:',seg_windows,"min started"))
  
  
#  if (!exists("compiled_imps")) {
#  compiled_imps <- lapply(seg_files,read_csv) %>% lapply(.,function(x){x %>% dplyr::select(-1,-2)})
#  }


  dir.create(file.path('..','all_motion_feature_data','ML_predictions',outcome),showWarnings = FALSE)
  dir.create(file.path('..','all_motion_feature_data','ML_predictions',outcome,paste0(seg_windows,'_min')),showWarnings = FALSE)
  
  # Preparing covariates for ML training
  #for (l in 1:length(compiled_imps)){
  l=1
  #curr_imp <- compiled_imps[[l]] %>% 

  outList<-seg_classification_function(curr_data_table,outcome,classifier_choice,patient_clinical_data,mf_col_idx,outer_fold_count,inner_fold_count)
  save(outList,file=file.path('..','all_motion_feature_data','ML_predictions',outcome,paste0(seg_windows,'_min'),paste0('imp',l,'.RData')))
  gc()
  print(paste('Segment size:',seg_windows,"min complete"))
}