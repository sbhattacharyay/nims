#### Master Script 6: Train and evaluate prediction models and measure feature significance scores ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# University of Cambridge
# Johns Hopkins University
# email address: sb2406@cam.ac.uk
#
### Contents:
# I. Initialization
# II.  Train and validate threshold-level GCSm detection models
# III. Train and validate threshold-level GOSE at discharge prediction models
# IV. Train and validate threshold-level GOSE at 12 months prediction models
# V. Perform LOL coefficient analysis for weighted feature analysis

### I. Initialization
# Load necessary packages
library(caret)
library(rms)
library(tidyverse)
library(tidymodels)
library(nnet)
library(glmnet)
library(MASS)
library(UBL)
library(ROCR)
library(foreach)
library(doParallel)

# Set the number of parallel cores
no.parallel.cores <- 10
registerDoParallel(cores = no.parallel.cores)

### II. Train and validate threshold-level GCSm detection models

# Directories of predictor matrix observation windows
obs.window.dirs <- list.files('../features/03_formatted_predictor_matrices/predictor_matrices','*_h_obs_window',include.dirs = T,full.names = T)

# Iterate through observation window directories
foreach(curr.obs.window.dir = obs.window.dirs) %dopar% {
  # Extract current observation window length
  curr.obs.window <- as.numeric(str_match(curr.obs.window.dir, "/03_formatted_predictor_matrices/predictor_matrices/\\s*(.*?)\\s*_h_obs_window")[,2])
  # Status update on current observation window
  print(paste("Observation window",curr.obs.window,"hour(s) started."))
  # Extract repeat directories in current observation window directory
  repeat.dirs <- list.files(curr.obs.window.dir,'repeat*',include.dirs = T,full.names = T)
  repeat.count <- length(repeat.dirs)
  # Get number of folds in each repeat directory
  fold.count <- length(list.files(repeat.dirs[1],'fold*',include.dirs = T,full.names = T))
  # For each of the 5 GCSm thresholds we test, we create a table of hyperparameter combinations and assign a configuration index
  viable.model.configurations <- vector(mode='list')
  curr.obs.window.compiled.predictions <- vector(mode='list') 
  for (temp.gcs.thresh in c('GCSm.gt.1','GCSm.gt.2','GCSm.gt.3','GCSm.gt.4','GCSm.gt.5')){
    viable.model.configurations[[temp.gcs.thresh]] <- as.data.frame(expand.grid(d = 2:20)) %>% mutate(config.idx = 1:nrow(.))
    curr.obs.window.compiled.predictions[[temp.gcs.thresh]] <- as.data.frame(matrix(ncol = 12,nrow=0))
  }
  
  # Iterate through repeat directories
  for (curr.repeat.dir in repeat.dirs){
    # Extract current repeat number
    curr.repeat.no <- as.numeric(sub(".*repeat", "", curr.repeat.dir))
    # Status update on current repeat
    print(paste("Repeat no.",curr.repeat.no,"out of",length(repeat.dirs), "started."))
    # Extract fold directories in current repeat directory
    fold.dirs <- list.files(curr.repeat.dir,'fold*',include.dirs = T,full.names = T)
    
    # Iterate through fold directories
    for (curr.fold.dir in fold.dirs){
      # Extract current fold number
      curr.fold.no <- as.numeric(sub(".*fold", "", curr.fold.dir))
      # Status update on current fold
      print(paste("Fold no.",curr.fold.no,"out of",length(fold.dirs), "started."))
      # Load training and testing keys of current fold
      curr.GCSm.train.keys <- read.csv(file.path(curr.fold.dir,'GCSm_train_key.csv'))
      curr.GCSm.test.keys <- read.csv(file.path(curr.fold.dir,'GCSm_test_key.csv'))
      
      # Iterate through the 5 GCSm thresholds
      for (curr.thresh.idx in 1:length(viable.model.configurations)){
        # Extract current threshold name
        curr.thresh.name <- names(viable.model.configurations)[curr.thresh.idx]
        # Status update on current threshold based on availability of predictors
        if (paste0('LOL_thresh_',curr.thresh.name,'_train_matrix.rds') %in% list.files(curr.fold.dir)){
          print(paste("Model training for threshold",curr.thresh.name,"started."))
        } else {
          print(paste("Model training for threshold",curr.thresh.name,"skipped due to missing predictors."))
          next
        }
        # Get current training predictors and labels
        curr.train.matrix <- readRDS(file.path(curr.fold.dir,paste0('LOL_thresh_',curr.thresh.name,'_train_matrix.rds')))
        curr.train.labels <- curr.GCSm.train.keys[[curr.thresh.name]]
        # Get current testing predictors and labels
        curr.test.matrix <- readRDS(file.path(curr.fold.dir,paste0('LOL_thresh_',curr.thresh.name,'_test_matrix.rds')))
        curr.test.labels <- curr.GCSm.test.keys[[curr.thresh.name]]
        ## Name the training and testing matrices, convert to dataframe, and add labels
        # train
        curr.train.matrix <- as.data.frame(curr.train.matrix)
        names(curr.train.matrix) <- paste(paste0(rep("MF.",ncol(curr.train.matrix)),1:ncol(curr.train.matrix)))
        curr.train.matrix[[curr.thresh.name]] <- curr.train.labels
        # test
        curr.test.matrix <- as.data.frame(curr.test.matrix)
        names(curr.test.matrix) <- paste(paste0(rep("MF.",ncol(curr.test.matrix)),1:ncol(curr.test.matrix)))
        curr.test.matrix[[curr.thresh.name]] <- curr.test.labels
        ## Perform preprocessing with `recipe` package
        # Create `recipe` object for Yeo-Johnson pre-processing and normalization
        curr.rec <- recipe(as.formula(paste(curr.thresh.name,"~.")), data = curr.train.matrix) %>%
          step_YeoJohnson(all_numeric_predictors()) %>%
          step_normalize(all_numeric_predictors())
        # Train pre-processing steps on training set
        trained.curr.rec <- prep(curr.rec, training = curr.train.matrix)
        # Transform both training and testing sets with the trained transformation
        curr.tf.train.matrix <- bake(trained.curr.rec, new_data = curr.train.matrix)
        curr.tf.test.matrix  <- bake(trained.curr.rec, new_data = curr.test.matrix)
        # Get current dataframe of viable hyperparameters
        curr.viable.model.configurations <- viable.model.configurations[[curr.thresh.name]]
        
        # Iterate through viable model configurations and train/validate models
        for (curr.config.row in 1:nrow(curr.viable.model.configurations)){
          # Get current target dimensionality (d)
          curr.d <- curr.viable.model.configurations$d[curr.config.row]
          # Get current config.index
          curr.config.idx <- curr.viable.model.configurations$config.idx[curr.config.row]
          # Define current formula
          curr.formula <- as.formula(paste(curr.thresh.name,"~",paste(paste0(rep("MF.",curr.d),1:(curr.d)),collapse = " + ")))
          # Train model depending on SMOTE status
          curr.mdl <- glm(curr.formula, data = curr.tf.train.matrix, family = "binomial")
          # Evaluate trained model on testing dataset
          curr.predictions.df <- data.frame(TrueLabel = curr.test.labels,
                                            Prob = predict(curr.mdl, newdata = curr.tf.test.matrix,type='response'),
                                            ConfigIdx = curr.config.idx,
                                            TargetDim = curr.d, 
                                            Threshold = curr.thresh.name,
                                            UPI = curr.GCSm.test.keys$UPI,
                                            HoursFromICUAdmission = curr.GCSm.test.keys$HoursFromICUAdmission,
                                            TimeOfDay = curr.GCSm.test.keys$TimeOfDay,
                                            Repeat = curr.repeat.no,
                                            Fold = curr.fold.no,
                                            ObsWindow = curr.obs.window)
          # Append validation set predictions to running compiled list
          curr.obs.window.compiled.predictions[[curr.thresh.name]] <- rbind(
            curr.obs.window.compiled.predictions[[curr.thresh.name]],
            curr.predictions.df
          )
          # Create new nested directory to save model
          dir.create(
            file.path(
              '../results/GCSm_threshold_prediction',
              paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
              paste0('repeat', sprintf('%02.f', curr.repeat.no)),
              paste0('fold', sprintf('%02.f', curr.fold.no)),
              curr.thresh.name
            ),recursive = T, showWarnings = F)
          saveRDS(curr.mdl,
                  file.path(
                    '../results/GCSm_threshold_prediction',
                    paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
                    paste0('repeat', sprintf('%02.f', curr.repeat.no)),
                    paste0('fold', sprintf('%02.f', curr.fold.no)),
                    curr.thresh.name,
                    paste0('model_config_',sprintf('%02.f',curr.config.idx),'.rds')))
        }
        
        # Status update on current threshold
        print(paste("Model training for threshold",curr.thresh.name,"completed."))
        
      }
      
      # Status update on current fold
      print(paste("Fold no.",curr.fold.no,"out of",length(fold.dirs), "completed."))
      
    }
    # Save current status of compiled predictions
    saveRDS(curr.obs.window.compiled.predictions,
            file.path(
              '../results/GCSm_threshold_prediction',
              paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
              paste0('repeat', sprintf('%02.f', curr.repeat.no)),
              paste0('compiled_predictions.rds')))
    # Status update on current repeat
    print(paste("Repeat no.",curr.repeat.no,"out of",length(repeat.dirs), "completed."))
    
  }
  for (temp.gcs.thresh.2 in 1:length(curr.obs.window.compiled.predictions)){
    # Current threshold name
    curr.name.of.thresh <- names(curr.obs.window.compiled.predictions)[temp.gcs.thresh.2]
    # Extract current threshold prediction dataframe
    curr.thresh.pred.df <- curr.obs.window.compiled.predictions[[temp.gcs.thresh.2]]
    # Save current threshold predictions as CSV
    write.csv(curr.thresh.pred.df,
              file.path(
                '../results/GCSm_threshold_prediction',
                paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
                paste0(curr.name.of.thresh,'_compiled_predictions.csv')),
              row.names = F)
  } 
  # Status update on current observation window
  print(paste("Observation window",curr.obs.window,"completed."))
}

# Stop implicit cluster
stopImplicitCluster()

### III. Train and validate threshold-level GOSE at discharge prediction models

# Directories of predictor matrix observation windows
obs.window.dirs <- list.files('../features/03_formatted_predictor_matrices/predictor_matrices','*_h_obs_window',include.dirs = T,full.names = T)

# Iterate through observation window directories
foreach(curr.obs.window.dir = obs.window.dirs) %dopar% {
  # Extract current observation window length
  curr.obs.window <- as.numeric(str_match(curr.obs.window.dir, "/03_formatted_predictor_matrices/predictor_matrices/\\s*(.*?)\\s*_h_obs_window")[,2])
  # Status update on current observation window
  print(paste("Observation window",curr.obs.window,"hour(s) started."))
  # Extract repeat directories in current observation window directory
  repeat.dirs <- list.files(curr.obs.window.dir,'repeat*',include.dirs = T,full.names = T)
  repeat.count <- length(repeat.dirs)
  # Get number of folds in each repeat directory
  fold.count <- length(list.files(repeat.dirs[1],'fold*',include.dirs = T,full.names = T))
  # For each of the 5 GOSE thresholds we test, we create a table of hyperparameter combinations and assign a configuration index
  viable.model.configurations <- vector(mode='list')
  curr.obs.window.compiled.predictions <- vector(mode='list') 
  for (temp.GOSE.thresh in c('GOSE.gt.1','GOSE.gt.2','GOSE.gt.3','GOSE.gt.4','GOSE.gt.5')){
    viable.model.configurations[[temp.GOSE.thresh]] <- as.data.frame(expand.grid(d = 2:20)) %>% mutate(config.idx = 1:nrow(.))
    curr.obs.window.compiled.predictions[[temp.GOSE.thresh]] <- as.data.frame(matrix(ncol = 12,nrow=0))
  }
  
  # Iterate through repeat directories
  for (curr.repeat.dir in repeat.dirs){
    # Extract current repeat number
    curr.repeat.no <- as.numeric(sub(".*repeat", "", curr.repeat.dir))
    # Status update on current repeat
    print(paste("Repeat no.",curr.repeat.no,"out of",length(repeat.dirs), "started."))
    # Extract fold directories in current repeat directory
    fold.dirs <- list.files(curr.repeat.dir,'fold*',include.dirs = T,full.names = T)
    
    # Iterate through fold directories
    for (curr.fold.dir in fold.dirs){
      # Extract current fold number
      curr.fold.no <- as.numeric(sub(".*fold", "", curr.fold.dir))
      # Status update on current fold
      print(paste("Fold no.",curr.fold.no,"out of",length(fold.dirs), "started."))
      # Load training and testing keys of current fold
      curr.GOSE.train.keys <- read.csv(file.path(curr.fold.dir,'GOSE_train_key.csv'))
      curr.GOSE.test.keys <- read.csv(file.path(curr.fold.dir,'GOSE_test_key.csv'))
      
      # Iterate through the 5 GOSE thresholds
      for (curr.thresh.idx in 1:length(viable.model.configurations)){
        # Extract current threshold name
        curr.thresh.name <- names(viable.model.configurations)[curr.thresh.idx]
        # Status update on current threshold based on availability of predictors
        if (paste0('LOL_thresh_',curr.thresh.name,'_train_matrix.rds') %in% list.files(curr.fold.dir)){
          print(paste("Model training for threshold",curr.thresh.name,"started."))
        } else {
          print(paste("Model training for threshold",curr.thresh.name,"skipped due to missing predictors."))
          next
        }
        # Get current training predictors and labels
        curr.train.matrix <- readRDS(file.path(curr.fold.dir,paste0('LOL_thresh_',curr.thresh.name,'_train_matrix.rds')))
        curr.train.labels <- curr.GOSE.train.keys[[curr.thresh.name]]
        # Get current testing predictors and labels
        curr.test.matrix <- readRDS(file.path(curr.fold.dir,paste0('LOL_thresh_',curr.thresh.name,'_test_matrix.rds')))
        curr.test.labels <- curr.GOSE.test.keys[[curr.thresh.name]]
        ## Name the training and testing matrices, convert to dataframe, and add labels
        # train
        curr.train.matrix <- as.data.frame(curr.train.matrix)
        names(curr.train.matrix) <- paste(paste0(rep("MF.",ncol(curr.train.matrix)),1:ncol(curr.train.matrix)))
        curr.train.matrix[[curr.thresh.name]] <- curr.train.labels
        # test
        curr.test.matrix <- as.data.frame(curr.test.matrix)
        names(curr.test.matrix) <- paste(paste0(rep("MF.",ncol(curr.test.matrix)),1:ncol(curr.test.matrix)))
        curr.test.matrix[[curr.thresh.name]] <- curr.test.labels
        ## Perform preprocessing with `recipe` package
        # Create `recipe` object for Yeo-Johnson pre-processing and normalization
        curr.rec <- recipe(as.formula(paste(curr.thresh.name,"~.")), data = curr.train.matrix) %>%
          step_YeoJohnson(all_numeric_predictors()) %>%
          step_normalize(all_numeric_predictors())
        # Train pre-processing steps on training set
        trained.curr.rec <- prep(curr.rec, training = curr.train.matrix)
        # Transform both training and testing sets with the trained transformation
        curr.tf.train.matrix <- bake(trained.curr.rec, new_data = curr.train.matrix)
        curr.tf.test.matrix  <- bake(trained.curr.rec, new_data = curr.test.matrix)
        # Get current dataframe of viable hyperparameters
        curr.viable.model.configurations <- viable.model.configurations[[curr.thresh.name]]
        
        # Iterate through viable model configurations and train/validate models
        for (curr.config.row in 1:nrow(curr.viable.model.configurations)){
          # Get current target dimensionality (d)
          curr.d <- curr.viable.model.configurations$d[curr.config.row]
          # Get current config.index
          curr.config.idx <- curr.viable.model.configurations$config.idx[curr.config.row]
          # Define current formula
          curr.formula <- as.formula(paste(curr.thresh.name,"~",paste(paste0(rep("MF.",curr.d),1:(curr.d)),collapse = " + ")))
          # Train model depending on SMOTE status
          curr.mdl <- glm(curr.formula, data = curr.tf.train.matrix, family = "binomial")
          # Evaluate trained model on testing dataset
          curr.predictions.df <- data.frame(TrueLabel = curr.test.labels,
                                            Prob = predict(curr.mdl, newdata = curr.tf.test.matrix,type='response'),
                                            ConfigIdx = curr.config.idx,
                                            TargetDim = curr.d, 
                                            Threshold = curr.thresh.name,
                                            UPI = curr.GOSE.test.keys$UPI,
                                            HoursFromICUAdmission = curr.GOSE.test.keys$HoursFromICUAdmission,
                                            TimeOfDay = curr.GOSE.test.keys$TimeOfDay,
                                            Repeat = curr.repeat.no,
                                            Fold = curr.fold.no,
                                            ObsWindow = curr.obs.window)
          # Append validation set predictions to running compiled list
          curr.obs.window.compiled.predictions[[curr.thresh.name]] <- rbind(
            curr.obs.window.compiled.predictions[[curr.thresh.name]],
            curr.predictions.df
          )
          # Create new nested directory to save model
          dir.create(
            file.path(
              '../results/GOSE_threshold_prediction',
              paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
              paste0('repeat', sprintf('%02.f', curr.repeat.no)),
              paste0('fold', sprintf('%02.f', curr.fold.no)),
              curr.thresh.name
            ),recursive = T, showWarnings = F)
          saveRDS(curr.mdl,
                  file.path(
                    '../results/GOSE_threshold_prediction',
                    paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
                    paste0('repeat', sprintf('%02.f', curr.repeat.no)),
                    paste0('fold', sprintf('%02.f', curr.fold.no)),
                    curr.thresh.name,
                    paste0('model_config_',sprintf('%02.f',curr.config.idx),'.rds')))
        }
        
        # Status update on current threshold
        print(paste("Model training for threshold",curr.thresh.name,"completed."))
        
      }
      
      # Status update on current fold
      print(paste("Fold no.",curr.fold.no,"out of",length(fold.dirs), "completed."))
      
    }
    # Save current status of compiled predictions
    saveRDS(curr.obs.window.compiled.predictions,
            file.path(
              '../results/GOSE_threshold_prediction',
              paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
              paste0('repeat', sprintf('%02.f', curr.repeat.no)),
              paste0('compiled_predictions.rds')))
    # Status update on current repeat
    print(paste("Repeat no.",curr.repeat.no,"out of",length(repeat.dirs), "completed."))
    
  }
  for (temp.GOSE.thresh.2 in 1:length(curr.obs.window.compiled.predictions)){
    # Current threshold name
    curr.name.of.thresh <- names(curr.obs.window.compiled.predictions)[temp.GOSE.thresh.2]
    # Extract current threshold prediction dataframe
    curr.thresh.pred.df <- curr.obs.window.compiled.predictions[[temp.GOSE.thresh.2]]
    # Save current threshold predictions as CSV
    write.csv(curr.thresh.pred.df,
              file.path(
                '../results/GOSE_threshold_prediction',
                paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
                paste0(curr.name.of.thresh,'_compiled_predictions.csv')),
              row.names = F)
  } 
  # Status update on current observation window
  print(paste("Observation window",curr.obs.window,"completed."))
}

# Stop implicit cluster
stopImplicitCluster()

### IV. Train and validate threshold-level GOSE at 12 months prediction models

# Directories of predictor matrix observation windows
obs.window.dirs <- list.files('../features/03_formatted_predictor_matrices/predictor_matrices','*_h_obs_window',include.dirs = T,full.names = T)

# Iterate through observation window directories
foreach(curr.obs.window.dir = obs.window.dirs) %dopar% {
  # Extract current observation window length
  curr.obs.window <- as.numeric(str_match(curr.obs.window.dir, "/03_formatted_predictor_matrices/predictor_matrices/\\s*(.*?)\\s*_h_obs_window")[,2])
  # Status update on current observation window
  print(paste("Observation window",curr.obs.window,"hour(s) started."))
  # Extract repeat directories in current observation window directory
  repeat.dirs <- list.files(curr.obs.window.dir,'repeat*',include.dirs = T,full.names = T)
  repeat.count <- length(repeat.dirs)
  # Get number of folds in each repeat directory
  fold.count <- length(list.files(repeat.dirs[1],'fold*',include.dirs = T,full.names = T))
  # For each of the 7 GOSE12m thresholds we test, we create a table of hyperparameter combinations and assign a configuration index
  viable.model.configurations <- vector(mode='list')
  curr.obs.window.compiled.predictions <- vector(mode='list') 
  for (temp.GOSE12m.thresh in c('GOSE12m.gt.1','GOSE12m.gt.2','GOSE12m.gt.3','GOSE12m.gt.4','GOSE12m.gt.5','GOSE12m.gt.6','GOSE12m.gt.7')){
    viable.model.configurations[[temp.GOSE12m.thresh]] <- as.data.frame(expand.grid(d = 2:20)) %>% mutate(config.idx = 1:nrow(.))
    curr.obs.window.compiled.predictions[[temp.GOSE12m.thresh]] <- as.data.frame(matrix(ncol = 12,nrow=0))
  }
  
  # Iterate through repeat directories
  for (curr.repeat.dir in repeat.dirs){
    # Extract current repeat number
    curr.repeat.no <- as.numeric(sub(".*repeat", "", curr.repeat.dir))
    # Status update on current repeat
    print(paste("Repeat no.",curr.repeat.no,"out of",length(repeat.dirs), "started."))
    # Extract fold directories in current repeat directory
    fold.dirs <- list.files(curr.repeat.dir,'fold*',include.dirs = T,full.names = T)
    
    # Iterate through fold directories
    for (curr.fold.dir in fold.dirs){
      # Extract current fold number
      curr.fold.no <- as.numeric(sub(".*fold", "", curr.fold.dir))
      # Status update on current fold
      print(paste("Fold no.",curr.fold.no,"out of",length(fold.dirs), "started."))
      # Load training and testing keys of current fold
      curr.GOSE12m.train.keys <- read.csv(file.path(curr.fold.dir,'GOSE12m_train_key.csv'))
      curr.GOSE12m.test.keys <- read.csv(file.path(curr.fold.dir,'GOSE12m_test_key.csv'))
      
      # Iterate through the 7 GOSE12m thresholds
      for (curr.thresh.idx in 1:length(viable.model.configurations)){
        # Extract current threshold name
        curr.thresh.name <- names(viable.model.configurations)[curr.thresh.idx]
        # Status update on current threshold based on availability of predictors
        if (paste0('LOL_thresh_',curr.thresh.name,'_train_matrix.rds') %in% list.files(curr.fold.dir)){
          print(paste("Model training for threshold",curr.thresh.name,"started."))
        } else {
          print(paste("Model training for threshold",curr.thresh.name,"skipped due to missing predictors."))
          next
        }
        # Get current training predictors and labels
        curr.train.matrix <- readRDS(file.path(curr.fold.dir,paste0('LOL_thresh_',curr.thresh.name,'_train_matrix.rds')))
        curr.train.labels <- curr.GOSE12m.train.keys[[curr.thresh.name]]
        # Get current testing predictors and labels
        curr.test.matrix <- readRDS(file.path(curr.fold.dir,paste0('LOL_thresh_',curr.thresh.name,'_test_matrix.rds')))
        curr.test.labels <- curr.GOSE12m.test.keys[[curr.thresh.name]]
        ## Name the training and testing matrices, convert to dataframe, and add labels
        # train
        curr.train.matrix <- as.data.frame(curr.train.matrix)
        names(curr.train.matrix) <- paste(paste0(rep("MF.",ncol(curr.train.matrix)),1:ncol(curr.train.matrix)))
        curr.train.matrix[[curr.thresh.name]] <- curr.train.labels
        # test
        curr.test.matrix <- as.data.frame(curr.test.matrix)
        names(curr.test.matrix) <- paste(paste0(rep("MF.",ncol(curr.test.matrix)),1:ncol(curr.test.matrix)))
        curr.test.matrix[[curr.thresh.name]] <- curr.test.labels
        ## Perform preprocessing with `recipe` package
        # Create `recipe` object for Yeo-Johnson pre-processing and normalization
        curr.rec <- recipe(as.formula(paste(curr.thresh.name,"~.")), data = curr.train.matrix) %>%
          step_YeoJohnson(all_numeric_predictors()) %>%
          step_normalize(all_numeric_predictors())
        # Train pre-processing steps on training set
        trained.curr.rec <- prep(curr.rec, training = curr.train.matrix)
        # Transform both training and testing sets with the trained transformation
        curr.tf.train.matrix <- bake(trained.curr.rec, new_data = curr.train.matrix)
        curr.tf.test.matrix  <- bake(trained.curr.rec, new_data = curr.test.matrix)
        # Get current dataframe of viable hyperparameters
        curr.viable.model.configurations <- viable.model.configurations[[curr.thresh.name]]
        
        # Iterate through viable model configurations and train/validate models
        for (curr.config.row in 1:nrow(curr.viable.model.configurations)){
          # Get current target dimensionality (d)
          curr.d <- curr.viable.model.configurations$d[curr.config.row]
          # Get current config.index
          curr.config.idx <- curr.viable.model.configurations$config.idx[curr.config.row]
          # Define current formula
          curr.formula <- as.formula(paste(curr.thresh.name,"~",paste(paste0(rep("MF.",curr.d),1:(curr.d)),collapse = " + ")))
          # Train model depending on SMOTE status
          curr.mdl <- glm(curr.formula, data = curr.tf.train.matrix, family = "binomial")
          # Evaluate trained model on testing dataset
          curr.predictions.df <- data.frame(TrueLabel = curr.test.labels,
                                            Prob = predict(curr.mdl, newdata = curr.tf.test.matrix,type='response'),
                                            ConfigIdx = curr.config.idx,
                                            TargetDim = curr.d, 
                                            Threshold = curr.thresh.name,
                                            UPI = curr.GOSE12m.test.keys$UPI,
                                            HoursFromICUAdmission = curr.GOSE12m.test.keys$HoursFromICUAdmission,
                                            TimeOfDay = curr.GOSE12m.test.keys$TimeOfDay,
                                            Repeat = curr.repeat.no,
                                            Fold = curr.fold.no,
                                            ObsWindow = curr.obs.window)
          # Append validation set predictions to running compiled list
          curr.obs.window.compiled.predictions[[curr.thresh.name]] <- rbind(
            curr.obs.window.compiled.predictions[[curr.thresh.name]],
            curr.predictions.df
          )
          # Create new nested directory to save model
          dir.create(
            file.path(
              '../results/GOSE12m_threshold_prediction',
              paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
              paste0('repeat', sprintf('%02.f', curr.repeat.no)),
              paste0('fold', sprintf('%02.f', curr.fold.no)),
              curr.thresh.name
            ),recursive = T, showWarnings = F)
          saveRDS(curr.mdl,
                  file.path(
                    '../results/GOSE12m_threshold_prediction',
                    paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
                    paste0('repeat', sprintf('%02.f', curr.repeat.no)),
                    paste0('fold', sprintf('%02.f', curr.fold.no)),
                    curr.thresh.name,
                    paste0('model_config_',sprintf('%02.f',curr.config.idx),'.rds')))
        }
        
        # Status update on current threshold
        print(paste("Model training for threshold",curr.thresh.name,"completed."))
        
      }
      
      # Status update on current fold
      print(paste("Fold no.",curr.fold.no,"out of",length(fold.dirs), "completed."))
      
    }
    # Save current status of compiled predictions
    saveRDS(curr.obs.window.compiled.predictions,
            file.path(
              '../results/GOSE12m_threshold_prediction',
              paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
              paste0('repeat', sprintf('%02.f', curr.repeat.no)),
              paste0('compiled_predictions.rds')))
    # Status update on current repeat
    print(paste("Repeat no.",curr.repeat.no,"out of",length(repeat.dirs), "completed."))
    
  }
  for (temp.GOSE12m.thresh.2 in 1:length(curr.obs.window.compiled.predictions)){
    # Current threshold name
    curr.name.of.thresh <- names(curr.obs.window.compiled.predictions)[temp.GOSE12m.thresh.2]
    # Extract current threshold prediction dataframe
    curr.thresh.pred.df <- curr.obs.window.compiled.predictions[[temp.GOSE12m.thresh.2]]
    # Save current threshold predictions as CSV
    write.csv(curr.thresh.pred.df,
              file.path(
                '../results/GOSE12m_threshold_prediction',
                paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
                paste0(curr.name.of.thresh,'_compiled_predictions.csv')),
              row.names = F)
  } 
  # Status update on current observation window
  print(paste("Observation window",curr.obs.window,"completed."))
}

# Stop implicit cluster
stopImplicitCluster()

### V. Perform LOL coefficient analysis for feature analysis
# Create directory to store feature analysis files
dir.create('../features/03_formatted_predictor_matrices/feature_analysis',showWarnings = F)

# Acquire list of observation window directories
obs.window.dirs <- list.files('../features/03_formatted_predictor_matrices/predictor_matrices','*_h_obs_window',include.dirs = T,full.names = T)

# Load patient outcome information
patient.outcomes <- read.csv('../clinical_data/patient_outcomes.csv')

# Create dataframe combination of sensors and features
sensor.feature.combos <- as.data.frame(expand.grid(sensor = c('LA','LE','LW','RA','RE','RW'),
                                                   feature = c('BPW','FDE','HLF_h','HLF_l','MFR','SMA','WVL','PhysActivity'))) %>%
  mutate(formatted.feature.name = paste0(sensor,'/',feature,'/'))

# Iterate through unique observation windows
foreach(curr.obs.window.dir = obs.window.dirs) %dopar% {
  
  # Extract observation window length based on current file name
  curr.obs.window.hours <- as.numeric(str_match(curr.obs.window.dir, "predictor_matrices/\\s*(.*?)\\s*_h_obs_window")[,2])
  
  # Status update on current observation window
  print(paste("Observation window",curr.obs.window.hours,"h started."))
  
  # Create directory base for saving all feature analysis outputs in current observation window
  curr.combo.directory <- paste0('../features/03_formatted_predictor_matrices/feature_analysis/',
                                 sprintf('%05.2f', curr.obs.window.hours),
                                 '_h_obs_window')
  dir.create(curr.combo.directory,showWarnings = F,recursive = T)
  
  # Load full feature matrix column names for current observation window
  curr.full.matrix.names <-
    colnames(readRDS(
      file.path(
        '../features/03_formatted_predictor_matrices/full_matrices',
        paste0(
          sprintf('%05.2f', curr.obs.window.hours),
          '_h_imputation_1_full_matrix.rds'
        )
      )
    ))
  
  ## Thresholded GCSm (dichotomous response)
  threshold.GCSm.names <- c('GCSm.gt.1','GCSm.gt.2','GCSm.gt.3','GCSm.gt.4','GCSm.gt.5')
  # Iterate through threshold outcomes and perform LOL on each outcome
  for (curr.threshold.name in threshold.GCSm.names){
    # Initialize dataframe for compiled threshold feature analysis coefficients
    compiled.thresh.GCSm.LOL.feature.analysis <- as.data.frame(matrix(ncol = 8,nrow = 0))
    # Load compiled metrics dataframe of current observation window/threshold combination and filter AUC results
    curr.compiled.AUC <-
      read.csv(file.path(
        '../results/GCSm_threshold_prediction',
        paste0(sprintf('%05.2f', curr.obs.window.hours), '_h_obs_window'),
        paste0(curr.threshold.name, '_compiled_metrics.csv')
      )) %>%
      filter(Metrics == 'AUC')
    # Find optimal configuration index based on AUC
    curr.optimal.config.idx <- as.numeric(names(which.max(table(curr.compiled.AUC$ConfigIdx))))
    # Extract list of trained LOL files for current observation window/threshold combination
    curr.LOL.file.list <- list.files(file.path(
      '../features/03_formatted_predictor_matrices/predictor_matrices',
      paste0(sprintf('%05.2f', curr.obs.window.hours), '_h_obs_window')),
      paste0('LOL_thresh_',curr.threshold.name, '.rds'),
      recursive = T,
      full.names = T)
    # Iterate through each trained LOL file
    for (curr.LOL.file in curr.LOL.file.list){
      
      # Load current LOL file
      curr.LOL.projection <- readRDS(curr.LOL.file)
      
      # Identify current repeat and fold number
      curr.repeat.no <- as.numeric(str_match(curr.LOL.file, "/repeat\\s*(.*?)\\s*/fold")[,2])
      curr.fold.no <- as.numeric(str_match(curr.LOL.file, "/fold\\s*(.*?)\\s*/LOL_thresh_")[,2])
      
      # Load glm object of current observation window, threshold, repeat, and fold combination
      curr.glm.object <- readRDS(file.path(
        '../results/GCSm_threshold_prediction',
        paste0(sprintf('%05.2f', curr.obs.window.hours), '_h_obs_window'),
        paste0('repeat',sprintf('%02.f',curr.repeat.no)),
        paste0('fold',sprintf('%02.f',curr.repeat.no)),
        curr.threshold.name,
        paste0('model_config_',sprintf('%02.f',curr.optimal.config.idx),'.rds')
      ))
      
      # Load non-intercept coefficients from GLM object
      curr.coeff.vector <- abs(curr.glm.object$coefficients[2:(curr.optimal.config.idx+2)])
      
      # Extract absolute LOL projection coefficients for feature analysis and append to running dataframe
      curr.thresh.GCSm.LOL.feature.analysis <- data.frame(names = curr.full.matrix.names,
                                                          Significance = rowMeans(abs(curr.LOL.projection$A[,1:(curr.optimal.config.idx+1)]) %*% diag(curr.coeff.vector))) %>%
        separate(names,sep = '/',into=c('Sensor','Feature','StepsBeforeGCS')) %>%
        mutate(ObservationWindow = curr.obs.window.hours,
               RepeatNo = curr.repeat.no,
               FoldNo = curr.fold.no,
               Threshold = curr.threshold.name) %>%
        relocate(ObservationWindow,RepeatNo,FoldNo,Threshold)
      
      compiled.thresh.GCSm.LOL.feature.analysis <- rbind(compiled.thresh.GCSm.LOL.feature.analysis,
                                                         curr.thresh.GCSm.LOL.feature.analysis)
    }
    write.csv(compiled.thresh.GCSm.LOL.feature.analysis,file.path(curr.combo.directory,paste0(curr.threshold.name,'_feature_analysis_values.csv')),row.names = F)
  }
  
  ## Thresholded GOSE (dichotomous response)
  threshold.GOSE.names <- c('GOSE.gt.1','GOSE.gt.2','GOSE.gt.3','GOSE.gt.4','GOSE.gt.5')
  # Iterate through threshold outcomes and perform LOL on each outcome
  for (curr.threshold.name in threshold.GOSE.names){
    # Initialize dataframe for compiled threshold feature analysis coefficients
    compiled.thresh.GOSE.LOL.feature.analysis <- as.data.frame(matrix(ncol = 8,nrow = 0))
    # Load compiled metrics dataframe of current observation window/threshold combination and filter AUC results
    curr.compiled.AUC <-
      read.csv(file.path(
        '../results/GOSE_threshold_prediction',
        paste0(sprintf('%05.2f', curr.obs.window.hours), '_h_obs_window'),
        paste0(curr.threshold.name, '_compiled_metrics.csv')
      )) %>%
      filter(Metrics == 'AUC')
    # Find optimal configuration index based on AUC
    curr.optimal.config.idx <- as.numeric(names(which.max(table(curr.compiled.AUC$ConfigIdx))))
    # Extract list of trained LOL files for current observation window/threshold combination
    curr.LOL.file.list <- list.files(file.path(
      '../features/03_formatted_predictor_matrices/predictor_matrices',
      paste0(sprintf('%05.2f', curr.obs.window.hours), '_h_obs_window')),
      paste0('LOL_thresh_',curr.threshold.name, '.rds'),
      recursive = T,
      full.names = T)
    # Iterate through each trained LOL file
    for (curr.LOL.file in curr.LOL.file.list){
      
      # Load current LOL file
      curr.LOL.projection <- readRDS(curr.LOL.file)
      
      # Identify current repeat and fold number
      curr.repeat.no <- as.numeric(str_match(curr.LOL.file, "/repeat\\s*(.*?)\\s*/fold")[,2])
      curr.fold.no <- as.numeric(str_match(curr.LOL.file, "/fold\\s*(.*?)\\s*/LOL_thresh_")[,2])
      
      # Load glm object of current observation window, threshold, repeat, and fold combination
      curr.glm.object <- readRDS(file.path(
        '../results/GOSE_threshold_prediction',
        paste0(sprintf('%05.2f', curr.obs.window.hours), '_h_obs_window'),
        paste0('repeat',sprintf('%02.f',curr.repeat.no)),
        paste0('fold',sprintf('%02.f',curr.repeat.no)),
        curr.threshold.name,
        paste0('model_config_',sprintf('%02.f',curr.optimal.config.idx),'.rds')
      ))
      
      # Load non-intercept coefficients from GLM object
      curr.coeff.vector <- abs(curr.glm.object$coefficients[2:(curr.optimal.config.idx+2)])
      
      # Extract absolute LOL projection coefficients for feature analysis and append to running dataframe
      curr.thresh.GOSE.LOL.feature.analysis <- data.frame(names = curr.full.matrix.names,
                                                          Significance = rowMeans(abs(curr.LOL.projection$A[,1:(curr.optimal.config.idx+1)]) %*% diag(curr.coeff.vector))) %>%
        separate(names,sep = '/',into=c('Sensor','Feature','StepsBeforeGCS')) %>%
        mutate(ObservationWindow = curr.obs.window.hours,
               RepeatNo = curr.repeat.no,
               FoldNo = curr.fold.no,
               Threshold = curr.threshold.name) %>%
        relocate(ObservationWindow,RepeatNo,FoldNo,Threshold)
      
      compiled.thresh.GOSE.LOL.feature.analysis <- rbind(compiled.thresh.GOSE.LOL.feature.analysis,
                                                         curr.thresh.GOSE.LOL.feature.analysis)
    }
    write.csv(compiled.thresh.GOSE.LOL.feature.analysis,file.path(curr.combo.directory,paste0(curr.threshold.name,'_feature_analysis_values.csv')),row.names = F)
  }
  
  ## Thresholded GOSE12m (dichotomous response)
  threshold.GOSE12m.names <- c('GOSE12m.gt.1','GOSE12m.gt.2','GOSE12m.gt.3','GOSE12m.gt.4','GOSE12m.gt.5','GOSE12m.gt.6','GOSE12m.gt.7')
  # Iterate through threshold outcomes and perform LOL on each outcome
  for (curr.threshold.name in threshold.GOSE12m.names){
    # Initialize dataframe for compiled threshold feature analysis coefficients
    compiled.thresh.GOSE12m.LOL.feature.analysis <- as.data.frame(matrix(ncol = 8,nrow = 0))
    # Load compiled metrics dataframe of current observation window/threshold combination and filter AUC results
    curr.compiled.AUC <-
      read.csv(file.path(
        '../results/GOSE12m_threshold_prediction',
        paste0(sprintf('%05.2f', curr.obs.window.hours), '_h_obs_window'),
        paste0(curr.threshold.name, '_compiled_metrics.csv')
      )) %>%
      filter(Metrics == 'AUC')
    # Find optimal configuration index based on AUC
    curr.optimal.config.idx <- as.numeric(names(which.max(table(curr.compiled.AUC$ConfigIdx))))
    # Extract list of trained LOL files for current observation window/threshold combination
    curr.LOL.file.list <- list.files(file.path(
      '../features/03_formatted_predictor_matrices/predictor_matrices',
      paste0(sprintf('%05.2f', curr.obs.window.hours), '_h_obs_window')),
      paste0('LOL_thresh_',curr.threshold.name, '.rds'),
      recursive = T,
      full.names = T)
    # Iterate through each trained LOL file
    for (curr.LOL.file in curr.LOL.file.list){
      
      # Load current LOL file
      curr.LOL.projection <- readRDS(curr.LOL.file)
      
      # Identify current repeat and fold number
      curr.repeat.no <- as.numeric(str_match(curr.LOL.file, "/repeat\\s*(.*?)\\s*/fold")[,2])
      curr.fold.no <- as.numeric(str_match(curr.LOL.file, "/fold\\s*(.*?)\\s*/LOL_thresh_")[,2])
      
      # Load glm object of current observation window, threshold, repeat, and fold combination
      curr.glm.object <- readRDS(file.path(
        '../results/GOSE12m_threshold_prediction',
        paste0(sprintf('%05.2f', curr.obs.window.hours), '_h_obs_window'),
        paste0('repeat',sprintf('%02.f',curr.repeat.no)),
        paste0('fold',sprintf('%02.f',curr.repeat.no)),
        curr.threshold.name,
        paste0('model_config_',sprintf('%02.f',curr.optimal.config.idx),'.rds')
      ))
      
      # Load non-intercept coefficients from GLM object
      curr.coeff.vector <- abs(curr.glm.object$coefficients[2:(curr.optimal.config.idx+2)])
      
      # Extract absolute LOL projection coefficients for feature analysis and append to running dataframe
      curr.thresh.GOSE12m.LOL.feature.analysis <- data.frame(names = curr.full.matrix.names,
                                                             Significance = rowMeans(abs(curr.LOL.projection$A[,1:(curr.optimal.config.idx+1)]) %*% diag(curr.coeff.vector))) %>%
        separate(names,sep = '/',into=c('Sensor','Feature','StepsBeforeGCS')) %>%
        mutate(ObservationWindow = curr.obs.window.hours,
               RepeatNo = curr.repeat.no,
               FoldNo = curr.fold.no,
               Threshold = curr.threshold.name) %>%
        relocate(ObservationWindow,RepeatNo,FoldNo,Threshold)
      
      compiled.thresh.GOSE12m.LOL.feature.analysis <- rbind(compiled.thresh.GOSE12m.LOL.feature.analysis,
                                                            curr.thresh.GOSE12m.LOL.feature.analysis)
    }
    write.csv(compiled.thresh.GOSE12m.LOL.feature.analysis,file.path(curr.combo.directory,paste0(curr.threshold.name,'_feature_analysis_values.csv')),row.names = F)
  }
}

# Stop implicit cluster
stopImplicitCluster()