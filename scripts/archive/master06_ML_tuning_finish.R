#### Master Script 6: Tuning ML Classification of Temporal Segments ####
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

# library(devtools)
# if (!require(lolR))
#   install_github('neurodata/lol',
#                  build_vignettes = TRUE,
#                  force = TRUE)
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

### Load clinical metadata ###
source('./functions/load_patient_clinical_data.R')
gose_thresh <- 5 # inclusive threshold
mrs_thresh <- 3 # inclusive threshold
patient_clinical_data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv',gose_thresh,mrs_thresh) %>% add_column(ptIdx = 1:nrow(.),.after = 2)

### Create and save outer folds for ML based on available outcome labels ###
labels.temp <- expand.grid(label=c("fav_mort","fav_GOSE","fav_mRS"),
                           temp = c("dis","12m")) %>% mutate(label.name = paste(label,temp,sep = "_"))
outer_fold_count <- 5
inner_fold_count <- 5
outerFolds <- readRDS('../all_motion_feature_data/outerFolds.rds')

# Define models to tune and train
classifier_choice <- c("lda","svmRadialWeights","parRF")

# Define tuning grids for FIRST TUNING RUN:

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

# Define the train control for all models
train.control <- trainControl(method="repeatedcv",
                              number=inner_fold_count,
                              repeats=1,
                              summaryFunction = twoClassSummary, 
                              classProbs = TRUE,
                              seeds = seed.list,  
                              savePredictions = "all",
                              returnResamp = "all")

### Run the model building function ###
source('./functions/finishTuneMachineLearningModels.R')
source('./functions/tuneDeepLearningModel.R')
source('./functions/saveRDSFiles.R')

iter <- 1 # iteration 1 for caret models
deep.iter <- 1 # iteration 1 for deep learning models
seg_window <- 30

path.save <- file.path("../all_motion_feature_data/ML_results",paste0(seg_window,'_min'))
dir.create(path.save,showWarnings = F,recursive = T)
finishTuneMachineLearningModels(seg_window = seg_window,
                                Iter = iter, 
                                DeepIter = deep.iter, 
                                classifier_choice = classifier_choice, 
                                seed.list = seed.list,
                                path.D = path.save, 
                                labelsList = labels.temp)

### Examine first-pass tuning outcomes ###

# determine which models have been run 
source('./functions/generateModelDF.R')

model.df.all <- generateModelDF(path.D = "../all_motion_feature_data/ML_results",outer_fold_count)

# combine results into dataframes: one dataframe for each combination of label &
# segment window:
source('./functions/buildPredResDF.R')
buildPredResDF(modelsDF = model.df.all,patient_clinical_data = patient_clinical_data,outerFolds,impNo = 1)

# d. Plot graphs to inspect tuning results and determine any further tuning required
plotTuningGraphs()
