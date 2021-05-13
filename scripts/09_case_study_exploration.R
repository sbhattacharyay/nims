#### Master Script 9: Case study exploration ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# University of Cambridge
# Johns Hopkins University
# email address: sb2406@cam.ac.uk
#
### Contents:
# I. Initialization
# II. Extract training and testing predictor sets for individual-level case studies
# III. Train and evalate GCSm < 4 detection models
# IV. Observe prediction trajectories for each candidate and build 95% CIs across imputations

### I. Initialization
# Load necessary packages
library(tidyverse)
library(lolR)
library(doParallel)
library(foreach)
library(tidymodels)

# Set step size for prediction on candidate patients
STEP.SIZE <- 10 # in minutes

### II. Extract training and testing predictor sets for individual-level case studies

## Identify candidates for case study
# Select patients who, during the time of accelerometry recording, have both GCSm > 4 and GCSm <= 4
candidates <- read.csv('../features/03_formatted_predictor_matrices/full_matrices/06.00_h_imputation_1_keys.csv') %>%
  mutate(GCSm.gt.4 = as.integer(GCSm > 4)) %>%
  group_by(UPI) %>%
  summarise(neg.cases = n() - sum(GCSm.gt.4),
            pos.cases = sum(GCSm.gt.4)) %>%
  filter(neg.cases >= 1 & pos.cases >= 1)

# Create a directory to store features for case study analysis
dir.create('../features/04_case_study_matrices/predictor_matrices_06.00_h',showWarnings = F,recursive = T)
dir.create('../features/04_case_study_matrices/predictor_matrices_00.45_h',showWarnings = F,recursive = T)

## Extract and store testing sets for each case study patient
# Identify list of full matrix files of 6-hour observation windows
full.matrix.file.list.6hr <- list.files('../features/03_formatted_predictor_matrices/full_matrices',
                                    glob2rx('06.00_h_imputation_*_full_matrix.rds'),
                                    full.names = T)

# Identify list of full matrix files of 27-min observation windows
full.matrix.file.list.27min <- list.files('../features/03_formatted_predictor_matrices/full_matrices',
                                          glob2rx('00.45_h_imputation_*_full_matrix.rds'),
                                          full.names = T)

# Create dataframe of combinations of sensors and features for normalization
sensor.feature.combos <- as.data.frame(expand.grid(sensor = c('LA','LE','LW','RA','RE','RW'),
                                                   feature = c('BPW','FDE','HLF_h','HLF_l','MFR','SMA','WVL','PhysActivity'))) %>%
  mutate(formatted.feature.name = paste0(sensor,'/',feature,'/'))

## 6-hour observation window:
# Iterate through full matrix files
for (curr.full.matrix.file in full.matrix.file.list.6hr){
  
  # From file name, extract current imputation number and observation window hours
  curr.imputation.no <- as.numeric(str_match(curr.full.matrix.file, "_h_imputation_\\s*(.*?)\\s*_full_matrix")[,2])
  curr.obs.window.hours <- as.numeric(str_match(curr.full.matrix.file, "full_matrices/\\s*(.*?)\\s*_h_imputation_")[,2])
  
  # Status update on current imputation
  print(paste0('Imputation no. ',curr.imputation.no,' out of ',length(full.matrix.file.list.6hr),' started'))
  
  # Load the row keys corresponding to the current full matrix
  curr.imputation.key <- read.csv(file.path('../features/03_formatted_predictor_matrices/full_matrices',
                                            paste0(sprintf('%05.2f',curr.obs.window.hours),
                                                   '_h_imputation_',
                                                   sprintf('%01.f',curr.imputation.no),
                                                   '_keys.csv'))) %>%
    mutate(MatrixIdx = 1:nrow(.))
  
  # Load the current full matrix
  curr.full.matrix <- readRDS(curr.full.matrix.file)
  
  # Iterate through each of the candidate patients
  for (i in 1:nrow(candidates)){
    
    # Status update on the current candidate case
    print(paste0('Candidate case ',i,' out of ',nrow(candidates),' started'))
    
    # Extract current candidate UPI
    curr.UPI <- candidates$UPI[i]
    
    # Create base directory of current candidate to store subsequent feature outputs
    curr.base.dir <- file.path('../features/04_case_study_matrices/predictor_matrices_06.00_h',curr.UPI)
    dir.create(curr.base.dir,showWarnings = F,recursive = T)
    
    # Remove candidate patient from current row keys to create keys for the training set
    filt.curr.imputation.key <- curr.imputation.key %>%
      filter(UPI != curr.UPI) %>%
      mutate(GCSm.gt.4 = as.integer(GCSm > 4))
    
    # Save training set key for current imputation
    write.csv(filt.curr.imputation.key,
              file.path(curr.base.dir,
                        paste0('training_key_imputation_',
                               sprintf('%01.f',curr.imputation.no),
                               '.csv')),
              row.names = F)
    
    # Extract the training matrix for the current candidate
    curr.train.matrix <- abs(curr.full.matrix[filt.curr.imputation.key$MatrixIdx,])
    
    # Create a dataframe to store means and standard deviations of the training set
    # sensor-feature combinations for later normalization of the testing set
    curr.train.mean.sd.info <- sensor.feature.combos %>%
      mutate(ImputationNo = curr.imputation.no,
             ObsWindow = curr.obs.window.hours,
             CaseStudyUPI = curr.UPI) %>%
      relocate(ImputationNo,ObsWindow,CaseStudyUPI) %>%
      mutate(Mean = NA,
             SD = NA)
    
    # Iterate through unique sensor-feature combinations
    for (curr.sf.combo.idx in 1:nrow(sensor.feature.combos)){
      # Get current sensor/feature combo
      curr.sf.combo <- sensor.feature.combos$formatted.feature.name[curr.sf.combo.idx]
      # Find indices of column names with current combination
      curr.sf.col.idx <- which(grepl(curr.sf.combo, colnames(curr.train.matrix), fixed = TRUE))
      # Get mean and standard deviation information from training matrix
      curr.sf.mean <- mean(curr.train.matrix[,curr.sf.col.idx],na.rm = T)
      curr.sf.std <- sd(curr.train.matrix[,curr.sf.col.idx],na.rm = T)
      # Save training set mean and standard deviations in the information dataframe
      curr.train.mean.sd.info$Mean[curr.sf.combo.idx] <- curr.sf.mean
      curr.train.mean.sd.info$SD[curr.sf.combo.idx] <- curr.sf.std
      # Transform training and testing sets accordingly
      curr.train.matrix[,curr.sf.col.idx] <- (curr.train.matrix[,curr.sf.col.idx] - curr.sf.mean)/curr.sf.std
    }
    
    # Save mean/standard deviation information for current imputation and candidate
    write.csv(curr.train.mean.sd.info,
              file.path(curr.base.dir,
                        paste0('train_mean_sd_imputation_',
                               sprintf('%01.f',curr.imputation.no),
                               '.csv')),
              row.names = F)
    
    # Train and save LOL projection
    curr.UPI.LOL <- lol.project.lol(curr.train.matrix,filt.curr.imputation.key[['GCSm.gt.4']],r = 20)
    saveRDS(curr.UPI.LOL$A,
            file.path(curr.base.dir,
                      paste0('LOL_transform_matrix_imputation_',
                             sprintf('%01.f',curr.imputation.no),
                             '.rds')))
    
    # Save LOL-embedded training matrix into the directory
    saveRDS(curr.UPI.LOL$Xr,
            file.path(curr.base.dir,
                      paste0('LOL_train_matrix_imputation_',
                             sprintf('%01.f',curr.imputation.no),
                             '.rds')))
    
  }
  # Status update on the current imputation
  print(paste0('Imputation no. ',curr.imputation.no,' out of ',length(full.matrix.file.list.6hr),' completed'))
}

## 27-minute observation window:
# Iterate through full matrix files
for (curr.full.matrix.file in full.matrix.file.list.27min){
  
  # From file name, extract current imputation number and observation window hours
  curr.imputation.no <- as.numeric(str_match(curr.full.matrix.file, "_h_imputation_\\s*(.*?)\\s*_full_matrix")[,2])
  curr.obs.window.hours <- as.numeric(str_match(curr.full.matrix.file, "full_matrices/\\s*(.*?)\\s*_h_imputation_")[,2])
  
  # Status update on current imputation
  print(paste0('Imputation no. ',curr.imputation.no,' out of ',length(full.matrix.file.list.27min),' started'))
  
  # Load the row keys corresponding to the current full matrix
  curr.imputation.key <- read.csv(file.path('../features/03_formatted_predictor_matrices/full_matrices',
                                            paste0(sprintf('%05.2f',curr.obs.window.hours),
                                                   '_h_imputation_',
                                                   sprintf('%01.f',curr.imputation.no),
                                                   '_keys.csv'))) %>%
    mutate(MatrixIdx = 1:nrow(.))
  
  # Load the current full matrix
  curr.full.matrix <- readRDS(curr.full.matrix.file)
  
  # Iterate through each of the candidate patients
  for (i in 1:nrow(candidates)){
    
    # Status update on the current candidate case
    print(paste0('Candidate case ',i,' out of ',nrow(candidates),' started'))
    
    # Extract current candidate UPI
    curr.UPI <- candidates$UPI[i]
    
    # Create base directory of current candidate to store subsequent feature outputs
    curr.base.dir <- file.path('../features/04_case_study_matrices/predictor_matrices_00.45_h',curr.UPI)
    dir.create(curr.base.dir,showWarnings = F,recursive = T)
    
    # Remove candidate patient from current row keys to create keys for the training set
    filt.curr.imputation.key <- curr.imputation.key %>%
      filter(UPI != curr.UPI) %>%
      mutate(GCSm.gt.4 = as.integer(GCSm > 4))
    
    # Save training set key for current imputation
    write.csv(filt.curr.imputation.key,
              file.path(curr.base.dir,
                        paste0('training_key_imputation_',
                               sprintf('%01.f',curr.imputation.no),
                               '.csv')),
              row.names = F)
    
    # Extract the training matrix for the current candidate
    curr.train.matrix <- abs(curr.full.matrix[filt.curr.imputation.key$MatrixIdx,])
    
    # Create a dataframe to store means and standard deviations of the training set
    # sensor-feature combinations for later normalization of the testing set
    curr.train.mean.sd.info <- sensor.feature.combos %>%
      mutate(ImputationNo = curr.imputation.no,
             ObsWindow = curr.obs.window.hours,
             CaseStudyUPI = curr.UPI) %>%
      relocate(ImputationNo,ObsWindow,CaseStudyUPI) %>%
      mutate(Mean = NA,
             SD = NA)
    
    # Iterate through unique sensor-feature combinations
    for (curr.sf.combo.idx in 1:nrow(sensor.feature.combos)){
      # Get current sensor/feature combo
      curr.sf.combo <- sensor.feature.combos$formatted.feature.name[curr.sf.combo.idx]
      # Find indices of column names with current combination
      curr.sf.col.idx <- which(grepl(curr.sf.combo, colnames(curr.train.matrix), fixed = TRUE))
      # Get mean and standard deviation information from training matrix
      curr.sf.mean <- mean(curr.train.matrix[,curr.sf.col.idx],na.rm = T)
      curr.sf.std <- sd(curr.train.matrix[,curr.sf.col.idx],na.rm = T)
      # Save training set mean and standard deviations in the information dataframe
      curr.train.mean.sd.info$Mean[curr.sf.combo.idx] <- curr.sf.mean
      curr.train.mean.sd.info$SD[curr.sf.combo.idx] <- curr.sf.std
      # Transform training and testing sets accordingly
      curr.train.matrix[,curr.sf.col.idx] <- (curr.train.matrix[,curr.sf.col.idx] - curr.sf.mean)/curr.sf.std
    }
    
    # Save mean/standard deviation information for current imputation and candidate
    write.csv(curr.train.mean.sd.info,
              file.path(curr.base.dir,
                        paste0('train_mean_sd_imputation_',
                               sprintf('%01.f',curr.imputation.no),
                               '.csv')),
              row.names = F)
    
    # Train and save LOL projection
    curr.UPI.LOL <- lol.project.lol(curr.train.matrix,filt.curr.imputation.key[['GCSm.gt.4']],r = 20)
    saveRDS(curr.UPI.LOL$A,
            file.path(curr.base.dir,
                      paste0('LOL_transform_matrix_imputation_',
                             sprintf('%01.f',curr.imputation.no),
                             '.rds')))
    
    # Save LOL-embedded training matrix into the directory
    saveRDS(curr.UPI.LOL$Xr,
            file.path(curr.base.dir,
                      paste0('LOL_train_matrix_imputation_',
                             sprintf('%01.f',curr.imputation.no),
                             '.rds')))
    
  }
  # Status update on the current imputation
  print(paste0('Imputation no. ',curr.imputation.no,' out of ',length(full.matrix.file.list.27min),' completed'))
}

## Extract and store testing sets for each case study patient
# Extract list of case study candidates
case.study.UPIs <- list.dirs('../features/04_case_study_matrices/predictor_matrices_06.00_h',recursive = F,full.names = F)

# Extract list of complete, bed-corrected, imputed feature sets
bed.corrected.imputation.files <- list.files('../features/02_bed_corrected_imputed_features',
                                             glob2rx('bed_corrected_imputation_*.csv'),
                                             full.names = T)

## 6-hour observation windows:
# Define number of columns for the testing set matrix based on observation window
obs.window <- 6 #in hours
num.columns <- round(((obs.window*720)+1)*7*6 + 6)
num.recording.idx <- ((6*720)+1)

# Iterate through each full feature set file
for (curr.bc.imputation.file in bed.corrected.imputation.files){
  
  # From the file name, extract the current imputation number
  curr.imputation.no <- as.numeric(str_match(curr.bc.imputation.file, "/bed_corrected_imputation_\\s*(.*?)\\s*.csv")[,2])
  
  # Status update on current imputation 
  print(paste0('Imputation no. ',curr.imputation.no,' out of ',length(bed.corrected.imputation.files),' started'))
  
  # Read large full feature set file of current imputation 
  curr.bc.imputation.features <- read.csv(curr.bc.imputation.file) 
  
  # Status update on loading of current imputation (requires 1+ min)
  print(paste('Imputation no.',curr.imputation.no,'file loaded'))
  
  # Iterate through candidate patients for case study
  for (curr.UPI in case.study.UPIs){
    
    # Status update on the current candidate case
    print(paste0('Candidate case ',curr.UPI,' started'))
    
    # Create base directory variable for saving outputs
    curr.base.dir <- file.path('../features/04_case_study_matrices/predictor_matrices_06.00_h',curr.UPI)
    
    # Load the LOL transform matrix of current candidate/imputation combination
    curr.LOL.transform.matrix <- readRDS(file.path(curr.base.dir,
                                                   paste0('LOL_transform_matrix_imputation_',
                                                          sprintf('%01.f',curr.imputation.no),
                                                          '.rds')))
    
    # Load current training set mean/sd info of current candidate/imputation combination
    curr.train.mean.sd.info <- read.csv(file.path(curr.base.dir,
                                                  paste0('train_mean_sd_imputation_',
                                                         sprintf('%01.f',curr.imputation.no),
                                                         '.csv')))
    
    # Filter out features of the candidate patient
    curr.UPI.features <- curr.bc.imputation.features %>%
      filter(UPI == curr.UPI)
    
    # Identify list of unique recording indices
    curr.unique.RecordingIdx <- unique(curr.UPI.features$RecordingIdx)
    
    # Define recording index endpoints based on defined step size
    curr.UPI.endpoints <- round(seq(num.recording.idx,max(curr.unique.RecordingIdx),by = STEP.SIZE*12))
    
    # Extract temporal information of the candidate patient time chunks
    curr.test.matrix.key <- curr.UPI.features %>% 
      dplyr::select(ImputationNo,UPI,HoursFromICUAdmission,TimeOfDay,RecordingIdx) %>%
      filter(RecordingIdx %in% curr.UPI.endpoints) %>%
      distinct() %>%
      dplyr::select(-RecordingIdx)
    
    # Save testing set key for current imputation
    write.csv(curr.test.matrix.key,
              file.path(curr.base.dir,
                        paste0('testing_key_imputation_',
                               sprintf('%01.f',curr.imputation.no),
                               '.csv')),
              row.names = F)
    
    # Initialize empty dataframe to store candidate patient features
    curr.test.matrix <- matrix(nrow = length(curr.UPI.endpoints),
                               ncol = num.columns)
    
    # Iterate through each endpoint for the current candidate
    for(curr.endpoint.idx in 1:length(curr.UPI.endpoints)){
      
      # Extract endpoint, starting point, and corresponding recording index information
      curr.endpoint <- curr.UPI.endpoints[curr.endpoint.idx]
      curr.startpoint <- round(curr.endpoint - (num.recording.idx - 1))
      curr.RecordingIdx <- seq(curr.startpoint,curr.endpoint)
      
      # Create dataframe to map recording indices to row indices for the matrix
      row.index.key <- data.frame(RecordingIdx = curr.RecordingIdx,
                                  RowIdx = rev(seq_along(curr.RecordingIdx))-1)
      
      # Filter out features for current endpoint and arrange appropriately
      curr.filt.UPI.features <- curr.UPI.features %>%
        filter(RecordingIdx %in% curr.RecordingIdx) %>% 
        left_join(row.index.key,by = 'RecordingIdx') %>%
        select(RowIdx,Feature,LA,LE,LW,RA,RE,RW) %>%
        pivot_longer(cols = -c(RowIdx,Feature), names_to = 'sensor') %>%
        relocate(sensor, Feature, RowIdx) %>%
        arrange(sensor, Feature, desc(RowIdx)) %>%
        distinct(sensor, Feature,RowIdx,.keep_all = T)
      
      # Calculate overall physical activity scores (based on SMA) and append to feature dataframe
      curr.physical.activity.scores <- curr.filt.UPI.features %>%
        filter(Feature == 'SMA') %>%
        group_by(sensor) %>%
        summarise(value = sum(value >= 0.135)/n()) %>%
        arrange(sensor) %>%
        mutate(Feature = 'PhysActivity',
               RowIdx = 'NotApplicable') %>%
        relocate(sensor,Feature,RowIdx,value)
      curr.filt.UPI.features <- rbind(curr.filt.UPI.features,curr.physical.activity.scores)
      
      # Place features of current endpoint into corresponding row of the matrix
      curr.test.matrix[curr.endpoint.idx,] <- curr.filt.UPI.features$value
      
      # Set column names of the matrix
      colnames(curr.test.matrix) <- paste(curr.filt.UPI.features$sensor,curr.filt.UPI.features$Feature,curr.filt.UPI.features$RowIdx,sep = '/')
      
      # Status update on current matrix row
      if (curr.endpoint.idx %% 10 == 0){
        print(paste("Matrix row",curr.endpoint.idx,"out of",length(curr.UPI.endpoints),"completed."))
      }
    }
    
    # Iterate through sensor/feature combinations to normalize testing matrix
    for (curr.sf.combo.idx in 1:nrow(curr.train.mean.sd.info)){
      # Get current sensor/feature combo
      curr.sf.combo <- curr.train.mean.sd.info$formatted.feature.name[curr.sf.combo.idx]
      # Find indices of column names with current combination
      curr.sf.col.idx <- which(grepl(curr.sf.combo, colnames(curr.test.matrix), fixed = TRUE))
      # Get mean and standard deviation information
      curr.sf.mean <- curr.train.mean.sd.info$Mean[curr.sf.combo.idx]
      curr.sf.std <- curr.train.mean.sd.info$SD[curr.sf.combo.idx]
      # Transform testing matrix accordingly
      curr.test.matrix[,curr.sf.col.idx] <- (curr.test.matrix[,curr.sf.col.idx] - curr.sf.mean)/curr.sf.std
    }
    
    # Reduce dimensionality of testing matrix with current transform matrix
    curr.LOL.test.matrix <- curr.test.matrix %*% curr.LOL.transform.matrix
    saveRDS(curr.LOL.test.matrix,
            file.path(curr.base.dir,
                      paste0('LOL_test_matrix_imputation_',
                             sprintf('%01.f',curr.imputation.no),
                             '.rds')))
    
    # Status update on the current candidate case
    print(paste0('Candidate case ',curr.UPI,' completed'))
  }
  # Status update on current imputation 
  print(paste0('Imputation no. ',curr.imputation.no,' out of ',length(bed.corrected.imputation.files),' started'))
}

## 27-min observation window:
# Iterate through each full feature set file
for (curr.bc.imputation.file in bed.corrected.imputation.files){
  
  # From the file name, extract the current imputation number
  curr.imputation.no <- as.numeric(str_match(curr.bc.imputation.file, "/bed_corrected_imputation_\\s*(.*?)\\s*.csv")[,2])
  
  # Status update on current imputation 
  print(paste0('Imputation no. ',curr.imputation.no,' out of ',length(bed.corrected.imputation.files),' started'))
  
  # Read large full feature set file of current imputation 
  curr.bc.imputation.features <- read.csv(curr.bc.imputation.file) 
  
  # Status update on loading of current imputation (requires 1+ min)
  print(paste('Imputation no.',curr.imputation.no,'file loaded'))
  
  # Iterate through candidate patients for case study
  for (curr.UPI in case.study.UPIs){
    
    # Status update on the current candidate case
    print(paste0('Candidate case ',curr.UPI,' started'))
    
    # Create base directory variable for saving outputs
    curr.base.dir <- file.path('../features/04_case_study_matrices/predictor_matrices_00.45_h',curr.UPI)
    
    # Load the LOL transform matrix of current candidate/imputation combination
    curr.LOL.transform.matrix <- readRDS(file.path(curr.base.dir,
                                                   paste0('LOL_transform_matrix_imputation_',
                                                          sprintf('%01.f',curr.imputation.no),
                                                          '.rds')))
    
    # Load current training set mean/sd info of current candidate/imputation combination
    curr.train.mean.sd.info <- read.csv(file.path(curr.base.dir,
                                                  paste0('train_mean_sd_imputation_',
                                                         sprintf('%01.f',curr.imputation.no),
                                                         '.csv')))
    
    # Filter out features of the candidate patient
    curr.UPI.features <- curr.bc.imputation.features %>%
      filter(UPI == curr.UPI)
    
    # Identify list of unique recording indices
    curr.unique.RecordingIdx <- unique(curr.UPI.features$RecordingIdx)
    
    # Define recording index endpoints based on defined step size
    curr.UPI.endpoints <- round(seq(num.recording.idx,max(curr.unique.RecordingIdx),by = STEP.SIZE*12))
    
    # Extract temporal information of the candidate patient time chunks
    curr.test.matrix.key <- curr.UPI.features %>% 
      dplyr::select(ImputationNo,UPI,HoursFromICUAdmission,TimeOfDay,RecordingIdx) %>%
      filter(RecordingIdx %in% curr.UPI.endpoints) %>%
      distinct() %>%
      dplyr::select(-RecordingIdx)
    
    # Save testing set key for current imputation
    write.csv(curr.test.matrix.key,
              file.path(curr.base.dir,
                        paste0('testing_key_imputation_',
                               sprintf('%01.f',curr.imputation.no),
                               '.csv')),
              row.names = F)
    
    # Initialize empty dataframe to store candidate patient features
    curr.test.matrix <- matrix(nrow = length(curr.UPI.endpoints),
                               ncol = num.columns)
    
    # Iterate through each endpoint for the current candidate
    for(curr.endpoint.idx in 1:length(curr.UPI.endpoints)){
      
      # Extract endpoint, starting point, and corresponding recording index information
      curr.endpoint <- curr.UPI.endpoints[curr.endpoint.idx]
      curr.startpoint <- round(curr.endpoint - (num.recording.idx - 1))
      curr.RecordingIdx <- seq(curr.startpoint,curr.endpoint)
      
      # Create dataframe to map recording indices to row indices for the matrix
      row.index.key <- data.frame(RecordingIdx = curr.RecordingIdx,
                                  RowIdx = rev(seq_along(curr.RecordingIdx))-1)
      
      # Filter out features for current endpoint and arrange appropriately
      curr.filt.UPI.features <- curr.UPI.features %>%
        filter(RecordingIdx %in% curr.RecordingIdx) %>% 
        left_join(row.index.key,by = 'RecordingIdx') %>%
        select(RowIdx,Feature,LA,LE,LW,RA,RE,RW) %>%
        pivot_longer(cols = -c(RowIdx,Feature), names_to = 'sensor') %>%
        relocate(sensor, Feature, RowIdx) %>%
        arrange(sensor, Feature, desc(RowIdx)) %>%
        distinct(sensor, Feature,RowIdx,.keep_all = T)
      
      # Calculate overall physical activity scores (based on SMA) and append to feature dataframe
      curr.physical.activity.scores <- curr.filt.UPI.features %>%
        filter(Feature == 'SMA') %>%
        group_by(sensor) %>%
        summarise(value = sum(value >= 0.135)/n()) %>%
        arrange(sensor) %>%
        mutate(Feature = 'PhysActivity',
               RowIdx = 'NotApplicable') %>%
        relocate(sensor,Feature,RowIdx,value)
      curr.filt.UPI.features <- rbind(curr.filt.UPI.features,curr.physical.activity.scores)
      
      # Place features of current endpoint into corresponding row of the matrix
      curr.test.matrix[curr.endpoint.idx,] <- curr.filt.UPI.features$value
      
      # Set column names of the matrix
      colnames(curr.test.matrix) <- paste(curr.filt.UPI.features$sensor,curr.filt.UPI.features$Feature,curr.filt.UPI.features$RowIdx,sep = '/')
      
      # Status update on current matrix row
      if (curr.endpoint.idx %% 10 == 0){
        print(paste("Matrix row",curr.endpoint.idx,"out of",length(curr.UPI.endpoints),"completed."))
      }
    }
    
    # Iterate through sensor/feature combinations to normalize testing matrix
    for (curr.sf.combo.idx in 1:nrow(curr.train.mean.sd.info)){
      # Get current sensor/feature combo
      curr.sf.combo <- curr.train.mean.sd.info$formatted.feature.name[curr.sf.combo.idx]
      # Find indices of column names with current combination
      curr.sf.col.idx <- which(grepl(curr.sf.combo, colnames(curr.test.matrix), fixed = TRUE))
      # Get mean and standard deviation information
      curr.sf.mean <- curr.train.mean.sd.info$Mean[curr.sf.combo.idx]
      curr.sf.std <- curr.train.mean.sd.info$SD[curr.sf.combo.idx]
      # Transform testing matrix accordingly
      curr.test.matrix[,curr.sf.col.idx] <- (curr.test.matrix[,curr.sf.col.idx] - curr.sf.mean)/curr.sf.std
    }
    
    # Reduce dimensionality of testing matrix with current transform matrix
    curr.LOL.test.matrix <- curr.test.matrix %*% curr.LOL.transform.matrix
    saveRDS(curr.LOL.test.matrix,
            file.path(curr.base.dir,
                      paste0('LOL_test_matrix_imputation_',
                             sprintf('%01.f',curr.imputation.no),
                             '.rds')))
    
    # Status update on the current candidate case
    print(paste0('Candidate case ',curr.UPI,' completed'))
  }
  # Status update on current imputation 
  print(paste0('Imputation no. ',curr.imputation.no,' out of ',length(bed.corrected.imputation.files),' started'))
}

### III. Train and evaluate GCSm < 4 detection models
## 6-hour observation window:
# Create new directory to store case-study-specific results
dir.create('../results/case_study_analysis_06.00_h',showWarnings = F)

# Extract list of case study candidates
case.study.UPIs <- list.dirs('../features/04_case_study_matrices/predictor_matrices_06.00_h',recursive = F,full.names = F)

# Set the number of parallel cores for model training
no.parallel.cores <- 6
registerDoParallel(cores = no.parallel.cores)

# Iterate through case study candidates
foreach(curr.UPI = case.study.UPIs) %dopar% {
  # Status update on the current candidate case
  print(paste0('Candidate case ',curr.UPI,' started'))
  
  # Create placeholder variable for current predictor matrix directory
  curr.base.pred.dir <- file.path('../features/04_case_study_matrices/predictor_matrices_06.00_h',curr.UPI)
  
  # Extract list of training matrices for current candidate
  train.matrix.files <- list.files(curr.base.pred.dir,
                                   pattern = glob2rx('LOL_train_matrix_imputation_*.rds'),
                                   full.names = T)
  
  # Create directory to store results of current candidate
  curr.base.results.dir <- file.path('../results/case_study_analysis_06.00_h',curr.UPI)
  dir.create(curr.base.results.dir, showWarnings = F)
  
  # Initialize dataframe to store predictions across imputations and model configurations
  curr.UPI.predictions <- data.frame(matrix(ncol=8,nrow=0))
  
  # Iterate through training matrix files
  for (curr.train.matrix.file in train.matrix.files){
    
    # From the file name, extract the current imputation number
    curr.imputation.no <- as.numeric(str_match(curr.train.matrix.file, "LOL_train_matrix_imputation_\\s*(.*?)\\s*.rds")[,2])
    
    # Status update on current imputation 
    print(paste0('Imputation no. ',curr.imputation.no,' out of ',length(train.matrix.files),' started'))
    
    # Load current training and testing matrices
    curr.train.matrix <- readRDS(file.path(curr.base.pred.dir,
                                           paste0('LOL_train_matrix_imputation_',
                                                  sprintf('%01.f',curr.imputation.no),
                                                  '.rds')))
    curr.test.matrix <- readRDS(file.path(curr.base.pred.dir,
                                          paste0('LOL_test_matrix_imputation_',
                                                 sprintf('%01.f',curr.imputation.no),
                                                 '.rds')))
    
    # Load current training and testing keys
    curr.train.key <- read.csv(file.path(curr.base.pred.dir,
                                         paste0('training_key_imputation_',
                                                sprintf('%01.f',curr.imputation.no),
                                                '.csv')))
    curr.test.key <- read.csv(file.path(curr.base.pred.dir,
                                        paste0('testing_key_imputation_',
                                               sprintf('%01.f',curr.imputation.no),
                                               '.csv')))
    
    ## Name the training and testing matrices, convert to dataframe, and add labels to training set
    # train
    curr.train.matrix <- as.data.frame(curr.train.matrix)
    names(curr.train.matrix) <- paste(paste0(rep("MF.",ncol(curr.train.matrix)),1:ncol(curr.train.matrix)))
    curr.train.matrix$GCSm.gt.4 <- curr.train.key$GCSm.gt.4
    # test
    curr.test.matrix <- as.data.frame(curr.test.matrix)
    names(curr.test.matrix) <- paste(paste0(rep("MF.",ncol(curr.test.matrix)),1:ncol(curr.test.matrix)))
    
    ## Perform preprocessing with `recipe` package
    # Create `recipe` object for Yeo-Johnson pre-processing and normalization
    curr.rec <- recipe(as.formula("GCSm.gt.4 ~."), data = curr.train.matrix) %>%
      step_YeoJohnson(all_numeric_predictors()) %>%
      step_normalize(all_numeric_predictors())
    # Train pre-processing steps on training set
    trained.curr.rec <- prep(curr.rec, training = curr.train.matrix)
    # Transform both training and testing sets with the trained transformation
    curr.tf.train.matrix <- bake(trained.curr.rec, new_data = curr.train.matrix)
    curr.tf.test.matrix  <- bake(trained.curr.rec, new_data = curr.test.matrix)
    
    # Iterate through possible target dimensions of LOL
    for (curr.d in 1:20){
      # Define current formula
      curr.formula <- as.formula(paste("GCSm.gt.4 ~",paste(paste0(rep("MF.",curr.d),1:(curr.d)),collapse = " + ")))
      
      # Train model depending on SMOTE status
      curr.mdl <- glm(curr.formula, data = curr.tf.train.matrix, family = "binomial")
      
      # Save current model
      saveRDS(curr.mdl,
              file.path(
                curr.base.results.dir,
                paste0('model_d',sprintf('%02.f',curr.d),'_imputation_',curr.imputation.no,'.rds')))
      
      # Evaluate trained model on testing dataset
      curr.predictions.df <- data.frame(Prob = predict(curr.mdl, newdata = curr.tf.test.matrix,type='response'),
                                        TargetDim = curr.d, 
                                        Threshold = 'GCSm.gt.4',
                                        UPI = curr.UPI,
                                        HoursFromICUAdmission = curr.test.key$HoursFromICUAdmission,
                                        TimeOfDay = curr.test.key$TimeOfDay,
                                        ImputationNo = curr.imputation.no,
                                        ObsWindow = 6)
      
      # Append current predictions dataframe to running compiled dataframe of current candidate
      curr.UPI.predictions <- rbind(curr.UPI.predictions, curr.predictions.df)
    }
    # Status update on current imputation 
    print(paste0('Imputation no. ',curr.imputation.no,' out of ',length(train.matrix.files),' started'))
  }
  # Save compiled predictions for current candidate
  write.csv(curr.UPI.predictions,
            file.path(curr.base.results.dir,'compiled_predictions.csv'),
            row.names = F)
}

## 27-min observation window:
# Create new directory to store case-study-specific results
dir.create('../results/case_study_analysis_00.45_h',showWarnings = F)

# Extract list of case study candidates
case.study.UPIs <- list.dirs('../features/04_case_study_matrices/predictor_matrices_00.45_h',recursive = F,full.names = F)

# Set the number of parallel cores for model training
no.parallel.cores <- length(case.study.UPIs)
registerDoParallel(cores = no.parallel.cores)

# Iterate through case study candidates
foreach(curr.UPI = case.study.UPIs) %dopar% {
  # Status update on the current candidate case
  print(paste0('Candidate case ',curr.UPI,' started'))
  
  # Create placeholder variable for current predictor matrix directory
  curr.base.pred.dir <- file.path('../features/04_case_study_matrices/predictor_matrices_00.45_h',curr.UPI)
  
  # Extract list of training matrices for current candidate
  train.matrix.files <- list.files(curr.base.pred.dir,
                                   pattern = glob2rx('LOL_train_matrix_imputation_*.rds'),
                                   full.names = T)
  
  # Create directory to store results of current candidate
  curr.base.results.dir <- file.path('../results/case_study_analysis_00.45_h',curr.UPI)
  dir.create(curr.base.results.dir, showWarnings = F)
  
  # Initialize dataframe to store predictions across imputations and model configurations
  curr.UPI.predictions <- data.frame(matrix(ncol=8,nrow=0))
  
  # Iterate through training matrix files
  for (curr.train.matrix.file in train.matrix.files){
    
    # From the file name, extract the current imputation number
    curr.imputation.no <- as.numeric(str_match(curr.train.matrix.file, "LOL_train_matrix_imputation_\\s*(.*?)\\s*.rds")[,2])
    
    # Status update on current imputation 
    print(paste0('Imputation no. ',curr.imputation.no,' out of ',length(train.matrix.files),' started'))
    
    # Load current training and testing matrices
    curr.train.matrix <- readRDS(file.path(curr.base.pred.dir,
                                           paste0('LOL_train_matrix_imputation_',
                                                  sprintf('%01.f',curr.imputation.no),
                                                  '.rds')))
    curr.test.matrix <- readRDS(file.path(curr.base.pred.dir,
                                          paste0('LOL_test_matrix_imputation_',
                                                 sprintf('%01.f',curr.imputation.no),
                                                 '.rds')))
    
    # Load current training and testing keys
    curr.train.key <- read.csv(file.path(curr.base.pred.dir,
                                         paste0('training_key_imputation_',
                                                sprintf('%01.f',curr.imputation.no),
                                                '.csv')))
    curr.test.key <- read.csv(file.path(curr.base.pred.dir,
                                        paste0('testing_key_imputation_',
                                               sprintf('%01.f',curr.imputation.no),
                                               '.csv')))
    
    ## Name the training and testing matrices, convert to dataframe, and add labels to training set
    # train
    curr.train.matrix <- as.data.frame(curr.train.matrix)
    names(curr.train.matrix) <- paste(paste0(rep("MF.",ncol(curr.train.matrix)),1:ncol(curr.train.matrix)))
    curr.train.matrix$GCSm.gt.4 <- curr.train.key$GCSm.gt.4
    # test
    curr.test.matrix <- as.data.frame(curr.test.matrix)
    names(curr.test.matrix) <- paste(paste0(rep("MF.",ncol(curr.test.matrix)),1:ncol(curr.test.matrix)))
    
    ## Perform preprocessing with `recipe` package
    # Create `recipe` object for Yeo-Johnson pre-processing and normalization
    curr.rec <- recipe(as.formula("GCSm.gt.4 ~."), data = curr.train.matrix) %>%
      step_YeoJohnson(all_numeric_predictors()) %>%
      step_normalize(all_numeric_predictors())
    # Train pre-processing steps on training set
    trained.curr.rec <- prep(curr.rec, training = curr.train.matrix)
    # Transform both training and testing sets with the trained transformation
    curr.tf.train.matrix <- bake(trained.curr.rec, new_data = curr.train.matrix)
    curr.tf.test.matrix  <- bake(trained.curr.rec, new_data = curr.test.matrix)
    
    # Iterate through possible target dimensions of LOL
    for (curr.d in 1:20){
      # Define current formula
      curr.formula <- as.formula(paste("GCSm.gt.4 ~",paste(paste0(rep("MF.",curr.d),1:(curr.d)),collapse = " + ")))
      
      # Train model depending on SMOTE status
      curr.mdl <- glm(curr.formula, data = curr.tf.train.matrix, family = "binomial")
      
      # Save current model
      saveRDS(curr.mdl,
              file.path(
                curr.base.results.dir,
                paste0('model_d',sprintf('%02.f',curr.d),'_imputation_',curr.imputation.no,'.rds')))
      
      # Evaluate trained model on testing dataset
      curr.predictions.df <- data.frame(Prob = predict(curr.mdl, newdata = curr.tf.test.matrix,type='response'),
                                        TargetDim = curr.d, 
                                        Threshold = 'GCSm.gt.4',
                                        UPI = curr.UPI,
                                        HoursFromICUAdmission = curr.test.key$HoursFromICUAdmission,
                                        TimeOfDay = curr.test.key$TimeOfDay,
                                        ImputationNo = curr.imputation.no,
                                        ObsWindow = obs.window)
      
      # Append current predictions dataframe to running compiled dataframe of current candidate
      curr.UPI.predictions <- rbind(curr.UPI.predictions, curr.predictions.df)
    }
    # Status update on current imputation 
    print(paste0('Imputation no. ',curr.imputation.no,' out of ',length(train.matrix.files),' started'))
  }
  # Save compiled predictions for current candidate
  write.csv(curr.UPI.predictions,
            file.path(curr.base.results.dir,'compiled_predictions.csv'),
            row.names = F)
}

# Stop implicit cluster once complete with parallel processing
stopImplicitCluster()

### IV. Observe prediction trajectories for each candidate and build 95% CIs across imputations
## 6-hour observation windows:
# Set the number of parallel cores for model training
no.parallel.cores <- 10
registerDoParallel(cores = no.parallel.cores)

# Set number of bootstraps
NUM.BOOTSTRAPS <- 10000

# Extract list of case study candidates
case.study.UPIs <- list.dirs('../features/04_case_study_matrices/predictor_matrices_06.00_h',recursive = F,full.names = F)
case.study.UPI.TimeCutoffs <- data.frame(UPI = sort(case.study.UPIs),
                                         TimeCutoff = c(272,308,49,43,312,48),
                                         Case = c('Case No. 1',
                                                  'Case No. 2',
                                                  'Case No. 3',
                                                  'Case No. 4',
                                                  'Case No. 5',
                                                  'Case No. 6'))

# Retrieve GCSm scores for candidate patients
case.study.gcs.data <- read.csv('../clinical_data/neurological_assessments.csv') %>%
  filter(UPI %in% case.study.UPIs, CoincidesWithAccelRecording == T) %>%
  drop_na(GCSm) %>%
  mutate(GCSm.gt.4 = as.integer(GCSm > 4)) %>% 
  dplyr::select(UPI,CoincidesWithAccelRecording,HoursFromICUAdmission,TimeOfDay,GCSm,GCSm.gt.4)

all.cases.compiled.CI <- data.frame(matrix(ncol = 8,nrow=0))
all.cases.GCSm.data <- data.frame(matrix(ncol = 7,nrow=0))

# Create directory to store plots of the current date
dir.create(file.path('../plots',Sys.Date()),showWarnings = F)

for (curr.UPI in case.study.UPIs){
  
  curr.base.dir <- file.path('../results/case_study_analysis_06.00_h',curr.UPI)
  
  curr.TimeCutoff <- case.study.UPI.TimeCutoffs$TimeCutoff[case.study.UPI.TimeCutoffs$UPI == curr.UPI]
  
  curr.case.label <- case.study.UPI.TimeCutoffs$Case[case.study.UPI.TimeCutoffs$UPI == curr.UPI]
  
  # Filter out the GCSm information of the current patient
  curr.UPI.GCSm.data <- case.study.gcs.data %>% 
    filter(UPI == curr.UPI) %>%
    filter(HoursFromICUAdmission >= curr.TimeCutoff) %>%
    mutate(Case = curr.case.label)
  
  # Extract compiled predictions of the current case study UPI
  curr.UPI.compiled.predictions <- read.csv(file.path(curr.base.dir,
                                                      'compiled_predictions.csv')) %>%
    filter(TargetDim == 2,HoursFromICUAdmission >= curr.TimeCutoff) %>%
    filter(HoursFromICUAdmission %in% sort(unique(.$HoursFromICUAdmission))[seq(1,length(unique(.$HoursFromICUAdmission)),2)])
  
  curr.UPI.compiled.bootstrap.predictions <- foreach(icount(NUM.BOOTSTRAPS), .combine=rbind) %dopar%{
    
    # Draw current imputation samples for current bootstrap resamples
    curr.imputation.resample <- sample(unique(curr.UPI.compiled.predictions$ImputationNo),length(unique(curr.UPI.compiled.predictions$ImputationNo)),replace = T)
    
    # Isolate samples of current resample imputations and calculate mean
    curr.resample.predictions <- curr.UPI.compiled.predictions %>%
      filter(ImputationNo %in% unique(curr.imputation.resample)) %>%
      group_by(UPI,HoursFromICUAdmission,TimeOfDay) %>%
      summarise(Prob = mean(Prob,na.rm = T))
    
    # Return currrent bootstrap resample means for overall combination
    curr.resample.predictions
  }
  
  # Calculate confidence intervals across the bootstrap resamples
  curr.UPI.compiled.CI <- curr.UPI.compiled.bootstrap.predictions %>%
    group_by(UPI,HoursFromICUAdmission,TimeOfDay) %>%
    summarise(meanProb = mean(Prob,na.rm = T),
              medianProb = median(Prob,na.rm = T),
              lowerProb = quantile(Prob,0.025,na.rm = T),
              upperProb = quantile(Prob,0.975,na.rm = T)) %>%
    mutate(Case = curr.case.label)
  
  curr.UPI.trajectories.plot <- ggplot(data = NULL,mapping = aes(x = HoursFromICUAdmission)) +
    geom_line(data = curr.UPI.compiled.CI, mapping = aes(y=meanProb)) +
    geom_ribbon(data = curr.UPI.compiled.CI, mapping = aes(ymin = lowerProb,ymax = upperProb),alpha = 0.2,fill = 'black') + 
    geom_hline(yintercept = .5,color='gray') +
    xlab('Hours from ICU Admission') +
    ylab('Pr(GCSm > 4)') +
    guides(linetype = guide_legend(title = ''),
           color = guide_legend(title = 'GCSm')) +
    scale_linetype_manual(breaks = c(0,1),labels = c('GCSm <= 4','GCSm > 4'),values = c('dashed','solid')) +
    geom_vline(data = curr.UPI.GCSm.data, mapping = aes(xintercept = HoursFromICUAdmission, color = factor(GCSm),linetype=factor(GCSm.gt.4))) +
    coord_cartesian(xlim = c(curr.TimeCutoff,max(curr.UPI.compiled.CI$HoursFromICUAdmission)),
                    ylim=c(0,1),
                    expand = F)+
    ggtitle(curr.UPI)+
    theme_bw() + 
    theme(panel.grid = element_blank(),
          legend.position = 'bottom',
          plot.title = element_text(hjust = .5, face = 'bold'))
  
  
  ggsave(file.path('../plots',Sys.Date(),paste0('06.00_h_trajectory_plot_',curr.UPI,'.pdf')),curr.UPI.trajectories.plot,device = 'pdf',units = 'in',dpi = 300,height = 8.25,width = 11.75)
  
  all.cases.compiled.CI <- rbind(all.cases.compiled.CI,curr.UPI.compiled.CI)
  all.cases.GCSm.data <- rbind(all.cases.GCSm.data,curr.UPI.GCSm.data)
}

# Save dataframes into case study analysis results folder
write.csv(all.cases.compiled.CI,'../results/case_study_analysis_06.00_h/compiled_case_study_predictions.csv',row.names = F)
write.csv(all.cases.GCSm.data,'../results/case_study_analysis_06.00_h/compiled_GCSm_evaluations.csv',row.names = F)

# Stop implicit cluster from parallel processing
stopImplicitCluster()

case.study.UPI.time.info <- case.study.UPI.TimeCutoffs %>%
  rename(startHoursFromICUAdmission = TimeCutoff) %>%
  full_join(all.cases.compiled.CI %>% 
              group_by(UPI,Case) %>%
              summarise(endHoursFromICUAdmission = max(HoursFromICUAdmission,na.rm = T)),
            by = c('UPI','Case')) %>%
  relocate(Case,UPI)

write.csv(case.study.UPI.time.info,'../results/case_study_analysis_06.00_h/case_study_time_limits.csv',row.names = F)

## 27-min observation window:
# Set the number of parallel cores for model training
no.parallel.cores <- 10
registerDoParallel(cores = no.parallel.cores)

# Set number of bootstraps
NUM.BOOTSTRAPS <- 10000

# Extract list of case study candidates
case.study.UPIs <- list.dirs('../features/04_case_study_matrices/predictor_matrices_00.45_h',recursive = F,full.names = F)
case.study.UPI.TimeCutoffs <- data.frame(UPI = sort(case.study.UPIs),
                                         TimeCutoff = c(272,308,49,43,312,48),
                                         Case = c('Case No. 1',
                                                  'Case No. 2',
                                                  'Case No. 3',
                                                  'Case No. 4',
                                                  'Case No. 5',
                                                  'Case No. 6'))

# Retrieve GCSm scores for candidate patients
case.study.gcs.data <- read.csv('../clinical_data/neurological_assessments.csv') %>%
  filter(UPI %in% case.study.UPIs, CoincidesWithAccelRecording == T) %>%
  drop_na(GCSm) %>%
  mutate(GCSm.gt.4 = as.integer(GCSm > 4)) %>% 
  dplyr::select(UPI,CoincidesWithAccelRecording,HoursFromICUAdmission,TimeOfDay,GCSm,GCSm.gt.4)

all.cases.compiled.CI <- data.frame(matrix(ncol = 8,nrow=0))
all.cases.GCSm.data <- data.frame(matrix(ncol = 7,nrow=0))

# Create directory to store plots of the current date
dir.create(file.path('../plots',Sys.Date()),showWarnings = F)

for (curr.UPI in case.study.UPIs){
  
  curr.base.dir <- file.path('../results/case_study_analysis_00.45_h',curr.UPI)
  
  curr.TimeCutoff <- case.study.UPI.TimeCutoffs$TimeCutoff[case.study.UPI.TimeCutoffs$UPI == curr.UPI]
  
  curr.case.label <- case.study.UPI.TimeCutoffs$Case[case.study.UPI.TimeCutoffs$UPI == curr.UPI]
  
  # Filter out the GCSm information of the current patient
  curr.UPI.GCSm.data <- case.study.gcs.data %>% 
    filter(UPI == curr.UPI) %>%
    filter(HoursFromICUAdmission >= curr.TimeCutoff) %>%
    mutate(Case = curr.case.label)
  
  # Extract compiled predictions of the current case study UPI
  curr.UPI.compiled.predictions <- read.csv(file.path(curr.base.dir,
                                                      'compiled_predictions.csv')) %>%
    filter(TargetDim == 2,HoursFromICUAdmission >= curr.TimeCutoff) %>%
    filter(HoursFromICUAdmission %in% sort(unique(.$HoursFromICUAdmission))[seq(1,length(unique(.$HoursFromICUAdmission)),2)])
  
  curr.UPI.compiled.bootstrap.predictions <- foreach(icount(NUM.BOOTSTRAPS), .combine=rbind) %dopar%{
    
    # Draw current imputation samples for current bootstrap resamples
    curr.imputation.resample <- sample(unique(curr.UPI.compiled.predictions$ImputationNo),length(unique(curr.UPI.compiled.predictions$ImputationNo)),replace = T)
    
    # Isolate samples of current resample imputations and calculate mean
    curr.resample.predictions <- curr.UPI.compiled.predictions %>%
      filter(ImputationNo %in% unique(curr.imputation.resample)) %>%
      group_by(UPI,HoursFromICUAdmission,TimeOfDay) %>%
      summarise(Prob = mean(Prob,na.rm = T))
    
    # Return currrent bootstrap resample means for overall combination
    curr.resample.predictions
  }
  
  # Calculate confidence intervals across the bootstrap resamples
  curr.UPI.compiled.CI <- curr.UPI.compiled.bootstrap.predictions %>%
    group_by(UPI,HoursFromICUAdmission,TimeOfDay) %>%
    summarise(meanProb = mean(Prob,na.rm = T),
              medianProb = median(Prob,na.rm = T),
              lowerProb = quantile(Prob,0.025,na.rm = T),
              upperProb = quantile(Prob,0.975,na.rm = T)) %>%
    mutate(Case = curr.case.label)
  
  curr.UPI.trajectories.plot <- ggplot(data = NULL,mapping = aes(x = HoursFromICUAdmission)) +
    geom_line(data = curr.UPI.compiled.CI, mapping = aes(y=meanProb)) +
    geom_ribbon(data = curr.UPI.compiled.CI, mapping = aes(ymin = lowerProb,ymax = upperProb),alpha = 0.2,fill = 'black') + 
    geom_hline(yintercept = .5,color='gray') +
    xlab('Hours from ICU Admission') +
    ylab('Pr(GCSm > 4)') +
    guides(linetype = guide_legend(title = ''),
           color = guide_legend(title = 'GCSm')) +
    scale_linetype_manual(breaks = c(0,1),labels = c('GCSm <= 4','GCSm > 4'),values = c('dashed','solid')) +
    geom_vline(data = curr.UPI.GCSm.data, mapping = aes(xintercept = HoursFromICUAdmission, color = factor(GCSm),linetype=factor(GCSm.gt.4))) +
    coord_cartesian(xlim = c(curr.TimeCutoff,max(curr.UPI.compiled.CI$HoursFromICUAdmission)),
                    ylim=c(0,1),
                    expand = F)+
    ggtitle(curr.UPI)+
    theme_bw() + 
    theme(panel.grid = element_blank(),
          legend.position = 'bottom',
          plot.title = element_text(hjust = .5, face = 'bold'))
  
  
  ggsave(file.path('../plots',Sys.Date(),paste0('00.45_h_trajectory_plot_',curr.UPI,'.pdf')),curr.UPI.trajectories.plot,device = 'pdf',units = 'in',dpi = 300,height = 8.25,width = 11.75)
  
  all.cases.compiled.CI <- rbind(all.cases.compiled.CI,curr.UPI.compiled.CI)
  all.cases.GCSm.data <- rbind(all.cases.GCSm.data,curr.UPI.GCSm.data)
}

# Save dataframes into case study analysis results folder
write.csv(all.cases.compiled.CI,'../results/case_study_analysis_00.45_h/compiled_case_study_predictions.csv',row.names = F)
write.csv(all.cases.GCSm.data,'../results/case_study_analysis_00.45_h/compiled_GCSm_evaluations.csv',row.names = F)

# Stop implicit cluster from parallel processing
stopImplicitCluster()

case.study.UPI.time.info <- case.study.UPI.TimeCutoffs %>%
  rename(startHoursFromICUAdmission = TimeCutoff) %>%
  full_join(all.cases.compiled.CI %>% 
              group_by(UPI,Case) %>%
              summarise(endHoursFromICUAdmission = max(HoursFromICUAdmission,na.rm = T)),
            by = c('UPI','Case')) %>%
  relocate(Case,UPI)

write.csv(case.study.UPI.time.info,'../results/case_study_analysis_00.45_h/case_study_time_limits.csv',row.names = F)

