#### Master Script 5: LOL embedding for dimensionality reduction ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# University of Cambridge
# Johns Hopkins University
# email address: sb2406@cam.ac.uk
#
### Contents:
# I. Initialization
# II. Compile features into formatted matrices for LOL
# III. Create predictor matrices for GCSm detection and perform LOL projection to reduce dimensionality
# IV. Create predictor matrices for GOSE prediction at discharge and perform LOL projection to reduce dimensionality
# V. Create predictor matrices for GOSE prediction at 12 months and perform LOL projection to reduce dimensionality

### I. Initialization
# Load necessary packages
library(tidyverse)
library(lolR)
library(foreach)
library(doParallel)

# Set the number of parallel cores
no.parallel.cores <- 10
registerDoParallel(cores = no.parallel.cores)

# Create directory to store formatted predictor matrices
dir.create("../features/03_formatted_predictor_matrices",showWarnings = FALSE)

### II. Compile features into formatted matrices for LOL
# Acquire list of full bed-corrected imputation files
imp.files <- list.files("../features/02_bed_corrected_imputed_features",glob2rx("bed_corrected_imputation_*.csv"),full.names = T)

# Acquire list of repeated cross-validation splits for each observation window for GCSm and GOSE (discharge and 12m)
GCSm.validation.split.files <- list.files('../validation_sampling','*_h_GCSm_folds.csv',full.names = T)
GOSE.validation.split.files <- list.files('../validation_sampling','*_h_GOSE_folds.csv',full.names = T)
GOSE.12m.validation.split.files <- list.files('../validation_sampling','*_h_GOSE12m_folds.csv',full.names = T)

# Create directory to store full compiled matrices
dir.create("../features/03_formatted_predictor_matrices/full_matrices",showWarnings = FALSE)

# Iterate through each large imputation file
for (curr.imp.file in imp.files){
  curr.imp.no <- as.integer(str_match(curr.imp.file, "bed_corrected_imputation_\\s*(.*?)\\s*.csv")[,2])
  # Status update on current imputation
  print(paste("Imputation no.",curr.imp.no,"out of",length(imp.files),"started."))
  # Load current imputation's bed-corrected features
  curr.imp.features <- read.csv(curr.imp.file)
  # Iterate through different observation window lengths
  for (curr.validation.split.file in GCSm.validation.split.files){
    curr.obs.window.hours <- as.numeric(str_match(curr.validation.split.file, "validation_sampling/\\s*(.*?)\\s*_h_GCSm_folds.csv")[,2])
    # Status update on current observation window length
    print(paste("Observation window",sprintf('%05.2f',curr.obs.window.hours),"started."))
    # Load current repeated cross-validation splits based on observation window length and extract unique GCS observations
    curr.validation.splits <- read.csv(curr.validation.split.file) %>% select(-c(Fold,Repeat,Split)) %>% unique()
    # Save unique observations as a 'key' for the corresponding matrix
    write.csv(curr.validation.splits,paste0('../features/03_formatted_predictor_matrices/full_matrices/',sprintf('%05.2f',curr.obs.window.hours),'_h_imputation_',curr.imp.no,'_keys.csv'),row.names = F)
    # Initialize matrix for storing compiled features based on dimensions dictated by observation window
    num.columns <- round(((curr.obs.window.hours*720)+1)*7*6 + 6)
    num.rows <- nrow(curr.validation.splits)
    curr.validation.split.matrix <- matrix(nrow = num.rows,
                                           ncol = num.columns)
    # In the rare case (base R bug), that the number of columns was not properly initialized:
    if (ncol(curr.validation.split.matrix) != num.columns){
      curr.validation.split.matrix <- cbind(curr.validation.split.matrix,rep(NA,num.rows))
    }
    matrix.row.idx <- 0 # Dummy variable to iteratively count matrix row
    # Iterate through patients and extract patient-specific motion features
    for (curr.UPI in unique(curr.validation.splits$UPI)){
      curr.UPI.validation.split <- curr.validation.splits %>% filter(UPI == curr.UPI)
      curr.UPI.imp.features <- curr.imp.features %>% filter(UPI == curr.UPI)
      # Iterate through observations of chosen patient and reshape patient features into matrix
      for (curr.row.of.UPI in 1:nrow(curr.UPI.validation.split)){
        matrix.row.idx <- matrix.row.idx + 1
        # Filter out features within current observation window
        curr.row.UPI.imp.features <- curr.UPI.imp.features %>% filter(HoursFromICUAdmission <= curr.UPI.validation.split$HoursFromICUAdmission[curr.row.of.UPI],
                                                                      HoursFromICUAdmission >= curr.UPI.validation.split$HoursFromICUAdmission[curr.row.of.UPI] - curr.obs.window.hours) %>%
          mutate(RowIdx = round(720*(curr.UPI.validation.split$HoursFromICUAdmission[curr.row.of.UPI] - HoursFromICUAdmission)))
        # In the (rare) case that the number of filtered features does not equal the matrix dimensionality:
        if (sum(!(1:(curr.obs.window.hours*720) %in% unique(curr.row.UPI.imp.features$RowIdx))) >= 1){
          missing.row.indices <- (1:(curr.obs.window.hours*720))[!(1:(curr.obs.window.hours*720) %in% unique(curr.row.UPI.imp.features$RowIdx))]
          # If an endpoint is missing, add the next point to the recording
          curr.min.recording.idx <- min(curr.row.UPI.imp.features$RecordingIdx)
          features.to.append <- curr.UPI.imp.features %>% 
            filter(RecordingIdx <= (curr.min.recording.idx-1)) %>%
            filter(RecordingIdx == max(RecordingIdx)) %>%
            mutate(RowIdx = (curr.obs.window.hours*720))
          curr.row.UPI.imp.features <- rbind(curr.row.UPI.imp.features,features.to.append)
        }
        # Arrange filtered feature set for vectorization into matrix row
        curr.row.UPI.imp.features <- curr.row.UPI.imp.features %>% 
          select(RowIdx,Feature,LA,LE,LW,RA,RE,RW) %>%
          pivot_longer(cols = -c(RowIdx,Feature), names_to = 'sensor') %>%
          relocate(sensor, Feature, RowIdx) %>%
          arrange(sensor, Feature, desc(RowIdx)) %>%
          distinct(sensor, Feature,RowIdx,.keep_all = T)
        
        # Calculate physical activity score per sensor
        curr.physical.activity.scores <- curr.row.UPI.imp.features %>%
          filter(Feature == 'SMA') %>%
          group_by(sensor) %>%
          summarise(value = sum(value >= 0.135)/n()) %>%
          arrange(sensor) %>%
          mutate(Feature = 'PhysActivity',
                 RowIdx = 'NotApplicable') %>%
          relocate(sensor,Feature,RowIdx,value)
        
        # Append physical activity scores to the filtered feature set dataframe
        curr.row.UPI.imp.features <- rbind(curr.row.UPI.imp.features,curr.physical.activity.scores)
        
        # Assign reshaped variables into corresponding matrix row and name matrix column based on sensor placement, feature type, and time before observation
        curr.validation.split.matrix[matrix.row.idx,] <- curr.row.UPI.imp.features$value
        colnames(curr.validation.split.matrix) <- paste(curr.row.UPI.imp.features$sensor,curr.row.UPI.imp.features$Feature,curr.row.UPI.imp.features$RowIdx,sep = '/')
        # Status update on current matrix row
        if (matrix.row.idx %% 10 == 0){
          print(paste("Matrix row",matrix.row.idx,"out of",nrow(curr.validation.splits),"completed."))
        }
      }
    }
    # Save complete feature matrix of current observation window and imputation combination
    saveRDS(curr.validation.split.matrix,paste0('../features/03_formatted_predictor_matrices/full_matrices/',sprintf('%05.2f',curr.obs.window.hours),'_h_imputation_',curr.imp.no,'_full_matrix.rds'))
    print(paste("Observation window",sprintf('%05.2f',curr.obs.window.hours),"completed."))
  }
  # Status update on current imputation
  print(paste("Imputation no.",curr.imp.no,"out of",length(imp.files),"completed."))
}

### III. Create predictor matrices for GCSm detection and perform LOL projection to reduce dimensionality
# Create directory to store GCSm predictor matrices
dir.create('../features/03_formatted_predictor_matrices/predictor_matrices',showWarnings = F)

# Load patient outcome information
patient.outcomes <- read.csv('../clinical_data/patient_outcomes.csv')

# Acquire list of repeated cross-validation splits for each observation window for GCSm and GOSE
GCSm.validation.split.files <- list.files('../validation_sampling','*_h_GCSm_folds.csv',full.names = T)

# Acquire list of full bed-corrected imputation files to calculate number of imputations
imp.files <- list.files("../features/02_bed_corrected_imputed_features",glob2rx("bed_corrected_imputation_*.csv"),full.names = T)
m.imputations <- length(imp.files)

# Create dataframe combination of sensors and features
sensor.feature.combos <- as.data.frame(expand.grid(sensor = c('LA','LE','LW','RA','RE','RW'),
                                                   feature = c('BPW','FDE','HLF_h','HLF_l','MFR','SMA','WVL','PhysActivity'))) %>%
  mutate(formatted.feature.name = paste0(sensor,'/',feature,'/'))

# Iterate through through different observation window lengths
for (curr.validation.split.file in GCSm.validation.split.files[length(GCSm.validation.split.files)]){
  curr.obs.window.hours <- as.numeric(str_match(curr.validation.split.file, "validation_sampling/\\s*(.*?)\\s*_h_GCSm_folds.csv")[,2])
  # Status update on current observation window
  print(paste("Observation window",sprintf('%05.2f',curr.obs.window.hours),"started."))
  # Load current repeated cross-validation splits based on observation window length
  curr.validation.splits <- read.csv(curr.validation.split.file)
  # Extract unique repeat/fold combinations and randomly assign (with replacement) an imputation to each one
  unique.repeat.fold.combos <- curr.validation.splits %>%
    select(Repeat,Fold) %>% 
    unique() %>%
    arrange(Repeat,Fold) %>%
    mutate(Imputation = sample(1:m.imputations, nrow(.), replace = T))
  # Save imputation key for repeat/fold combinations in current observation window
  write.csv(unique.repeat.fold.combos,paste0('../validation_sampling/',sprintf('%05.2f',curr.obs.window.hours),'_h_GCSm_imputation_key.csv'),row.names = F)
  # Based on random imputation assignments, find imputations that must be loaded
  imps.to.load <- unique(unique.repeat.fold.combos$Imputation)
  # Load and save appropriate imputations (and their keys) into a list
  curr.obs.window.matrices <- vector(mode = "list")
  curr.obs.window.keys <- vector(mode = "list")
  for (curr.imp in imps.to.load){
    curr.obs.window.matrices[[curr.imp]] <- readRDS(paste0('../features/03_formatted_predictor_matrices/full_matrices/',
                                                           sprintf('%05.2f',curr.obs.window.hours),
                                                           '_h_imputation_',
                                                           curr.imp,
                                                           '_full_matrix.rds'))
    curr.obs.window.keys[[curr.imp]] <- read.csv(paste0('../features/03_formatted_predictor_matrices/full_matrices/',
                                                        sprintf('%05.2f',curr.obs.window.hours),
                                                        '_h_imputation_',
                                                        curr.imp,
                                                        '_keys.csv'))
  }
  # Iterate through unique combinations of repeats and folds
  for (curr.combo.idx in 1:nrow(unique.repeat.fold.combos)){
    # Status update on current repeat/fold combination
    print(paste("Repeat/Fold combination no.",curr.combo.idx,"out of",nrow(unique.repeat.fold.combos), "started."))
    # Extract repeat, fold, and imputation information based on current combination
    curr.repeat <- unique.repeat.fold.combos$Repeat[curr.combo.idx]
    curr.fold <- unique.repeat.fold.combos$Fold[curr.combo.idx]
    curr.imputation <- unique.repeat.fold.combos$Imputation[curr.combo.idx]
    # Create directory base for saving all outputs in current repeat/fold combination
    curr.combo.directory <- paste0('../features/03_formatted_predictor_matrices/predictor_matrices/',
                                   sprintf('%05.2f',curr.obs.window.hours),
                                   '_h_obs_window/repeat',
                                   sprintf('%02.f',curr.repeat),
                                   '/fold',
                                   sprintf('%02.f',curr.fold))
    # Create directory to store current repeat/fold combo information
    dir.create(curr.combo.directory,showWarnings = F,recursive = T)
    # Separate current cross-validation observation information for training and testing 
    curr.train.splits <- curr.validation.splits %>%
      filter(Fold == curr.fold,
             Repeat == curr.repeat,
             Split == 'Train')
    curr.test.splits <- curr.validation.splits %>%
      filter(Fold == curr.fold,
             Repeat == curr.repeat,
             Split == 'Test')
    
    # Extract full matrix based on the drawn imputation of the current repeat/fold combination
    curr.full.matrix <- curr.obs.window.matrices[[curr.imputation]]
    curr.full.matrix.key <- curr.obs.window.keys[[curr.imputation]] %>%
      mutate(MatrixRowIdx = 1:nrow(.))
    
    # Separate training and testing keys of current repeat/fold combination and encode outcome (GCS and GOSE) information
    curr.train.filt.matrix.key <- inner_join(curr.full.matrix.key,curr.train.splits,
                                             by = c("UPI", 
                                                    "HoursFromICUAdmission", 
                                                    "TimeOfDay", 
                                                    "GCST", 
                                                    "GCSm", 
                                                    "GCSe", 
                                                    "GCSv", 
                                                    "CoincidesWithAccelRecording")) %>%
      mutate(GCSm.gt.1 = as.integer(GCSm > 1),
             GCSm.gt.2 = as.integer(GCSm > 2),
             GCSm.gt.3 = as.integer(GCSm > 3),
             GCSm.gt.4 = as.integer(GCSm > 4),
             GCSm.gt.5 = as.integer(GCSm > 5)) %>%
      left_join(patient.outcomes %>% select(UPI,GOSEDischarge),by='UPI') %>%
      mutate(GOSE.gt.1 = as.integer(GOSEDischarge > 1),
             GOSE.gt.2 = as.integer(GOSEDischarge > 2),
             GOSE.gt.3 = as.integer(GOSEDischarge > 3),
             GOSE.gt.4 = as.integer(GOSEDischarge > 4),
             GOSE.gt.5 = as.integer(GOSEDischarge > 5))
    
    curr.test.filt.matrix.key <- inner_join(curr.full.matrix.key,curr.test.splits,
                                            by = c("UPI", 
                                                   "HoursFromICUAdmission", 
                                                   "TimeOfDay", 
                                                   "GCST", 
                                                   "GCSm", 
                                                   "GCSe", 
                                                   "GCSv", 
                                                   "CoincidesWithAccelRecording")) %>%
      mutate(GCSm.gt.1 = as.integer(GCSm > 1),
             GCSm.gt.2 = as.integer(GCSm > 2),
             GCSm.gt.3 = as.integer(GCSm > 3),
             GCSm.gt.4 = as.integer(GCSm > 4),
             GCSm.gt.5 = as.integer(GCSm > 5)) %>%
      left_join(patient.outcomes %>% select(UPI,GOSEDischarge),by='UPI') %>%
      mutate(GOSE.gt.1 = as.integer(GOSEDischarge > 1),
             GOSE.gt.2 = as.integer(GOSEDischarge > 2),
             GOSE.gt.3 = as.integer(GOSEDischarge > 3),
             GOSE.gt.4 = as.integer(GOSEDischarge > 4),
             GOSE.gt.5 = as.integer(GOSEDischarge > 5))
    
    # Save current training and testing keys
    write.csv(curr.train.filt.matrix.key,file.path(curr.combo.directory,'GCSm_train_key.csv'),row.names = F)
    write.csv(curr.test.filt.matrix.key,file.path(curr.combo.directory,'GCSm_test_key.csv'),row.names = F)
    
    # Separate full training and testing matrices and scale columns (each unique combination of feature and sensor type) based on training set
    curr.train.matrix <- abs(curr.full.matrix[curr.train.filt.matrix.key$MatrixRowIdx,])
    curr.test.matrix <- abs(curr.full.matrix[curr.test.filt.matrix.key$MatrixRowIdx,])
    
    for (curr.sf.combo.idx in 1:nrow(sensor.feature.combos)){
      # Get current sensor/feature combo
      curr.sf.combo <- sensor.feature.combos$formatted.feature.name[curr.sf.combo.idx]
      # Find indices of column names with current combination
      curr.sf.col.idx <- which(grepl(curr.sf.combo, colnames(curr.train.matrix), fixed = TRUE))
      # Get mean and standard deviation information from training matrix
      curr.sf.mean <- mean(curr.train.matrix[,curr.sf.col.idx],na.rm = T)
      curr.sf.std <- sd(curr.train.matrix[,curr.sf.col.idx],na.rm = T)
      # Transform training and testing sets accordingly
      curr.train.matrix[,curr.sf.col.idx] <- (curr.train.matrix[,curr.sf.col.idx] - curr.sf.mean)/curr.sf.std
      curr.test.matrix[,curr.sf.col.idx] <- (curr.test.matrix[,curr.sf.col.idx] - curr.sf.mean)/curr.sf.std
    }
    
    # Save full testing and training matrices to appropriate directory
    saveRDS(curr.train.matrix,file.path(curr.combo.directory,'GCSm_full_train_matrix.rds'))
    saveRDS(curr.test.matrix,file.path(curr.combo.directory,'GCSm_full_test_matrix.rds'))
    
    ## Full GCSm (ordinal response)
    # Train and save LOL projection
    curr.full.GCSm.LOL <- lol.project.lol(curr.train.matrix,curr.train.filt.matrix.key$GCSm,r = 20)
    saveRDS(curr.full.GCSm.LOL,file.path(curr.combo.directory,'LOL_full_GCSm.rds'))
    # Separate LOL-embedded training and testing matrices and save into the directory
    curr.full.GCSm.LOL.train.matrix <- curr.full.GCSm.LOL$Xr
    saveRDS(curr.full.GCSm.LOL.train.matrix,file.path(curr.combo.directory,'LOL_full_GCSm_train_matrix.rds'))
    curr.full.GCSm.LOL.test.matrix <- curr.test.matrix %*% curr.full.GCSm.LOL$A
    saveRDS(curr.full.GCSm.LOL.test.matrix,file.path(curr.combo.directory,'LOL_full_GCSm_test_matrix.rds'))
    
    ## Thresholded GCSm (dichotomous response)
    # Extract GCSm thresholded response variable names
    threshold.GCSm.names <- names(curr.train.filt.matrix.key)[startsWith(names(curr.train.filt.matrix.key),'GCSm.gt.')]
    # Iterate through threshold outcomes and perform LOL on each outcome
    for (curr.threshold.name in threshold.GCSm.names){
      # Skip current threshold if only one outcome case available 
      if (length(unique(curr.train.filt.matrix.key[[curr.threshold.name]])) == 1){
        print(paste0('Only one outcome for ',curr.threshold.name,' in repeat/fold combination no. ',curr.combo.idx))
        next
      }
      # Train and save LOL projection
      curr.thresh.GCSm.LOL <- lol.project.lol(curr.train.matrix,curr.train.filt.matrix.key[[curr.threshold.name]],r = 20)
      saveRDS(curr.thresh.GCSm.LOL,file.path(curr.combo.directory,paste0('LOL_thresh_',curr.threshold.name,'.rds')))
      # Separate LOL-embedded training and testing matrices and save into the directory
      curr.thresh.GCSm.LOL.train.matrix <- curr.thresh.GCSm.LOL$Xr
      saveRDS(curr.thresh.GCSm.LOL.train.matrix,file.path(curr.combo.directory,paste0('LOL_thresh_',curr.threshold.name,'_train_matrix.rds')))
      curr.thresh.GCSm.LOL.test.matrix <- curr.test.matrix %*% curr.thresh.GCSm.LOL$A
      saveRDS(curr.thresh.GCSm.LOL.test.matrix,file.path(curr.combo.directory,paste0('LOL_thresh_',curr.threshold.name,'_test_matrix.rds')))
    }
    
    # Status update on current repeat/fold combination
    print(paste("Repeat/Fold combination no.",curr.combo.idx,"out of",nrow(unique.repeat.fold.combos), "completed."))
  }
  # Status update on current observation window
  print(paste("Observation window",sprintf('%05.2f',curr.obs.window.hours),"completed."))
}

### IV. Create predictor matrices for GOSE prediction at discharge and perform LOL projection to reduce dimensionality
# Create directory to store GOSE predictor matrices
dir.create('../features/03_formatted_predictor_matrices/predictor_matrices',showWarnings = F)

# Load patient outcome information
patient.outcomes <- read.csv('../clinical_data/patient_outcomes.csv')

# Acquire list of repeated cross-validation splits for each observation window for GOSE
GOSE.validation.split.files <- list.files('../validation_sampling','*_h_GOSE_folds.csv',full.names = T)

# Acquire list of full bed-corrected imputation files to calculate number of imputations
imp.files <- list.files("../features/02_bed_corrected_imputed_features",glob2rx("bed_corrected_imputation_*.csv"),full.names = T)
m.imputations <- length(imp.files)

# Create dataframe combination of sensors and features
sensor.feature.combos <- as.data.frame(expand.grid(sensor = c('LA','LE','LW','RA','RE','RW'),
                                                   feature = c('BPW','FDE','HLF_h','HLF_l','MFR','SMA','WVL','PhysActivity'))) %>%
  mutate(formatted.feature.name = paste0(sensor,'/',feature,'/'))

# Iterate through through different observation window lengths
for (curr.validation.split.file in GOSE.validation.split.files){
  curr.obs.window.hours <- as.numeric(str_match(curr.validation.split.file, "validation_sampling/\\s*(.*?)\\s*_h_GOSE_folds.csv")[,2])
  # Status update on current observation window
  print(paste("Observation window",sprintf('%05.2f',curr.obs.window.hours),"started."))
  # Load current repeated cross-validation splits based on observation window length
  curr.validation.splits <- read.csv(curr.validation.split.file)
  # Extract unique repeat/fold combinations and randomly assign (with replacement) an imputation to each one
  unique.repeat.fold.combos <- curr.validation.splits %>%
    select(Repeat,Fold) %>% 
    unique() %>%
    arrange(Repeat,Fold) %>%
    mutate(Imputation = sample(1:m.imputations, nrow(.), replace = T))
  # Save imputation key for repeat/fold combinations in current observation window
  write.csv(unique.repeat.fold.combos,paste0('../validation_sampling/',sprintf('%05.2f',curr.obs.window.hours),'_h_GOSE_imputation_key.csv'),row.names = F)
  # Based on random imputation assignments, find imputations that must be loaded
  imps.to.load <- unique(unique.repeat.fold.combos$Imputation)
  # Load and save appropriate imputations (and their keys) into a list
  curr.obs.window.matrices <- vector(mode = "list")
  curr.obs.window.keys <- vector(mode = "list")
  for (curr.imp in imps.to.load){
    curr.obs.window.matrices[[curr.imp]] <- readRDS(paste0('../features/03_formatted_predictor_matrices/full_matrices/',
                                                           sprintf('%05.2f',curr.obs.window.hours),
                                                           '_h_imputation_',
                                                           curr.imp,
                                                           '_full_matrix.rds'))
    curr.obs.window.keys[[curr.imp]] <- read.csv(paste0('../features/03_formatted_predictor_matrices/full_matrices/',
                                                        sprintf('%05.2f',curr.obs.window.hours),
                                                        '_h_imputation_',
                                                        curr.imp,
                                                        '_keys.csv'))
  }
  # Iterate through unique combinations of repeats and folds
  for (curr.combo.idx in 1:nrow(unique.repeat.fold.combos)){
    # Status update on current repeat/fold combination
    print(paste("Repeat/Fold combination no.",curr.combo.idx,"out of",nrow(unique.repeat.fold.combos), "started."))
    # Extract repeat, fold, and imputation information based on current combination
    curr.repeat <- unique.repeat.fold.combos$Repeat[curr.combo.idx]
    curr.fold <- unique.repeat.fold.combos$Fold[curr.combo.idx]
    curr.imputation <- unique.repeat.fold.combos$Imputation[curr.combo.idx]
    # Create directory base for saving all outputs in current repeat/fold combination
    curr.combo.directory <- paste0('../features/03_formatted_predictor_matrices/predictor_matrices/',
                                   sprintf('%05.2f',curr.obs.window.hours),
                                   '_h_obs_window/repeat',
                                   sprintf('%02.f',curr.repeat),
                                   '/fold',
                                   sprintf('%02.f',curr.fold))
    # Create directory to store current repeat/fold combo information
    dir.create(curr.combo.directory,showWarnings = F,recursive = T)
    # Separate current cross-validation observation information for training and testing 
    curr.train.splits <- curr.validation.splits %>%
      filter(Fold == curr.fold,
             Repeat == curr.repeat,
             Split == 'Train')
    curr.test.splits <- curr.validation.splits %>%
      filter(Fold == curr.fold,
             Repeat == curr.repeat,
             Split == 'Test')
    
    # Extract full matrix based on the drawn imputation of the current repeat/fold combination
    curr.full.matrix <- curr.obs.window.matrices[[curr.imputation]]
    curr.full.matrix.key <- curr.obs.window.keys[[curr.imputation]] %>%
      mutate(MatrixRowIdx = 1:nrow(.))
    
    # Separate training and testing keys of current repeat/fold combination and encode outcome (GCS and GOSE) information
    curr.train.filt.matrix.key <- inner_join(curr.full.matrix.key,curr.train.splits,
                                             by = c("UPI", 
                                                    "HoursFromICUAdmission", 
                                                    "TimeOfDay", 
                                                    "GCST", 
                                                    "GCSm", 
                                                    "GCSe", 
                                                    "GCSv", 
                                                    "CoincidesWithAccelRecording")) %>%
      mutate(GCSm.gt.1 = as.integer(GCSm > 1),
             GCSm.gt.2 = as.integer(GCSm > 2),
             GCSm.gt.3 = as.integer(GCSm > 3),
             GCSm.gt.4 = as.integer(GCSm > 4),
             GCSm.gt.5 = as.integer(GCSm > 5)) %>%
      left_join(patient.outcomes %>% select(UPI,GOSEDischarge),by='UPI') %>%
      mutate(GOSE.gt.1 = as.integer(GOSEDischarge > 1),
             GOSE.gt.2 = as.integer(GOSEDischarge > 2),
             GOSE.gt.3 = as.integer(GOSEDischarge > 3),
             GOSE.gt.4 = as.integer(GOSEDischarge > 4),
             GOSE.gt.5 = as.integer(GOSEDischarge > 5))
    
    curr.test.filt.matrix.key <- inner_join(curr.full.matrix.key,curr.test.splits,
                                            by = c("UPI", 
                                                   "HoursFromICUAdmission", 
                                                   "TimeOfDay", 
                                                   "GCST", 
                                                   "GCSm", 
                                                   "GCSe", 
                                                   "GCSv", 
                                                   "CoincidesWithAccelRecording")) %>%
      mutate(GCSm.gt.1 = as.integer(GCSm > 1),
             GCSm.gt.2 = as.integer(GCSm > 2),
             GCSm.gt.3 = as.integer(GCSm > 3),
             GCSm.gt.4 = as.integer(GCSm > 4),
             GCSm.gt.5 = as.integer(GCSm > 5)) %>%
      left_join(patient.outcomes %>% select(UPI,GOSEDischarge),by='UPI') %>%
      mutate(GOSE.gt.1 = as.integer(GOSEDischarge > 1),
             GOSE.gt.2 = as.integer(GOSEDischarge > 2),
             GOSE.gt.3 = as.integer(GOSEDischarge > 3),
             GOSE.gt.4 = as.integer(GOSEDischarge > 4),
             GOSE.gt.5 = as.integer(GOSEDischarge > 5))
    
    # Save current training and testing keys
    write.csv(curr.train.filt.matrix.key,file.path(curr.combo.directory,'GOSE_train_key.csv'),row.names = F)
    write.csv(curr.test.filt.matrix.key,file.path(curr.combo.directory,'GOSE_test_key.csv'),row.names = F)
    
    # Separate full training and testing matrices and scale columns (each unique combination of feature and sensor type) based on training set
    curr.train.matrix <- abs(curr.full.matrix[curr.train.filt.matrix.key$MatrixRowIdx,])
    curr.test.matrix <- abs(curr.full.matrix[curr.test.filt.matrix.key$MatrixRowIdx,])
    
    for (curr.sf.combo.idx in 1:nrow(sensor.feature.combos)){
      # Get current sensor/feature combo
      curr.sf.combo <- sensor.feature.combos$formatted.feature.name[curr.sf.combo.idx]
      # Find indices of column names with current combination
      curr.sf.col.idx <- which(grepl(curr.sf.combo, colnames(curr.train.matrix), fixed = TRUE))
      # Get mean and standard deviation information from training matrix
      curr.sf.mean <- mean(curr.train.matrix[,curr.sf.col.idx],na.rm = T)
      curr.sf.std <- sd(curr.train.matrix[,curr.sf.col.idx],na.rm = T)
      # Transform training and testing sets accordingly
      curr.train.matrix[,curr.sf.col.idx] <- (curr.train.matrix[,curr.sf.col.idx] - curr.sf.mean)/curr.sf.std
      curr.test.matrix[,curr.sf.col.idx] <- (curr.test.matrix[,curr.sf.col.idx] - curr.sf.mean)/curr.sf.std
    }
    
    # Save full testing and training matrices to appropriate directory
    saveRDS(curr.train.matrix,file.path(curr.combo.directory,'GOSE_full_train_matrix.rds'))
    saveRDS(curr.test.matrix,file.path(curr.combo.directory,'GOSE_full_test_matrix.rds'))
    
    ## Full GOSE (ordinal response)
    # Train and save LOL projection
    curr.full.GOSE.LOL <- lol.project.lol(curr.train.matrix,curr.train.filt.matrix.key$GOSEDischarge,r = 20)
    saveRDS(curr.full.GOSE.LOL,file.path(curr.combo.directory,'LOL_full_GOSE.rds'))
    # Separate LOL-embedded training and testing matrices and save into the directory
    curr.full.GOSE.LOL.train.matrix <- curr.full.GOSE.LOL$Xr
    saveRDS(curr.full.GOSE.LOL.train.matrix,file.path(curr.combo.directory,'LOL_full_GOSE_train_matrix.rds'))
    curr.full.GOSE.LOL.test.matrix <- curr.test.matrix %*% curr.full.GOSE.LOL$A
    saveRDS(curr.full.GOSE.LOL.test.matrix,file.path(curr.combo.directory,'LOL_full_GOSE_test_matrix.rds'))
    
    ## Thresholded GOSE (dichotomous response)
    # Extract GOSE thresholded response variable names
    threshold.GOSE.names <- names(curr.train.filt.matrix.key)[startsWith(names(curr.train.filt.matrix.key),'GOSE.gt.')]
    # Iterate through threshold outcomes and perform LOL on each outcome
    for (curr.threshold.name in threshold.GOSE.names){
      # Skip current threshold if only one outcome case available 
      if (length(unique(curr.train.filt.matrix.key[[curr.threshold.name]])) == 1){
        print(paste0('Only one outcome for ',curr.threshold.name,' in repeat/fold combination no. ',curr.combo.idx))
        next
      }
      # Train and save LOL projection
      curr.thresh.GOSE.LOL <- lol.project.lol(curr.train.matrix,curr.train.filt.matrix.key[[curr.threshold.name]],r = 20)
      saveRDS(curr.thresh.GOSE.LOL,file.path(curr.combo.directory,paste0('LOL_thresh_',curr.threshold.name,'.rds')))
      # Separate LOL-embedded training and testing matrices and save into the directory
      curr.thresh.GOSE.LOL.train.matrix <- curr.thresh.GOSE.LOL$Xr
      saveRDS(curr.thresh.GOSE.LOL.train.matrix,file.path(curr.combo.directory,paste0('LOL_thresh_',curr.threshold.name,'_train_matrix.rds')))
      curr.thresh.GOSE.LOL.test.matrix <- curr.test.matrix %*% curr.thresh.GOSE.LOL$A
      saveRDS(curr.thresh.GOSE.LOL.test.matrix,file.path(curr.combo.directory,paste0('LOL_thresh_',curr.threshold.name,'_test_matrix.rds')))
    }
    
    # Status update on current repeat/fold combination
    print(paste("Repeat/Fold combination no.",curr.combo.idx,"out of",nrow(unique.repeat.fold.combos), "completed."))
  }
  # Status update on current observation window
  print(paste("Observation window",sprintf('%05.2f',curr.obs.window.hours),"completed."))
}

### V. Create predictor matrices for GOSE prediction at 12 months and perform LOL projection to reduce dimensionality
# Create directory to store GOSE12m predictor matrices
dir.create('../features/03_formatted_predictor_matrices/predictor_matrices',showWarnings = F)

# Load patient outcome information
patient.outcomes <- read.csv('../clinical_data/patient_outcomes.csv')

# Acquire list of repeated cross-validation splits for each observation window for GOSE12m
GOSE12m.validation.split.files <- list.files('../validation_sampling','*_h_GOSE12m_folds.csv',full.names = T)

# Acquire list of full bed-corrected imputation files to calculate number of imputations
imp.files <- list.files("../features/02_bed_corrected_imputed_features",glob2rx("bed_corrected_imputation_*.csv"),full.names = T)
m.imputations <- length(imp.files)

# Create dataframe combination of sensors and features
sensor.feature.combos <- as.data.frame(expand.grid(sensor = c('LA','LE','LW','RA','RE','RW'),
                                                   feature = c('BPW','FDE','HLF_h','HLF_l','MFR','SMA','WVL','PhysActivity'))) %>%
  mutate(formatted.feature.name = paste0(sensor,'/',feature,'/'))

# Iterate through through different observation window lengths
for (curr.validation.split.file in GOSE12m.validation.split.files){
  curr.obs.window.hours <- as.numeric(str_match(curr.validation.split.file, "validation_sampling/\\s*(.*?)\\s*_h_GOSE12m_folds.csv")[,2])
  # Status update on current observation window
  print(paste("Observation window",sprintf('%05.2f',curr.obs.window.hours),"started."))
  # Load current repeated cross-validation splits based on observation window length
  curr.validation.splits <- read.csv(curr.validation.split.file)
  # Extract unique repeat/fold combinations and randomly assign (with replacement) an imputation to each one
  unique.repeat.fold.combos <- curr.validation.splits %>%
    select(Repeat,Fold) %>% 
    unique() %>%
    arrange(Repeat,Fold) %>%
    mutate(Imputation = sample(1:m.imputations, nrow(.), replace = T))
  # Save imputation key for repeat/fold combinations in current observation window
  write.csv(unique.repeat.fold.combos,paste0('../validation_sampling/',sprintf('%05.2f',curr.obs.window.hours),'_h_GOSE12m_imputation_key.csv'),row.names = F)
  # Based on random imputation assignments, find imputations that must be loaded
  imps.to.load <- unique(unique.repeat.fold.combos$Imputation)
  # Load and save appropriate imputations (and their keys) into a list
  curr.obs.window.matrices <- vector(mode = "list")
  curr.obs.window.keys <- vector(mode = "list")
  for (curr.imp in imps.to.load){
    curr.obs.window.matrices[[curr.imp]] <- readRDS(paste0('../features/03_formatted_predictor_matrices/full_matrices/',
                                                           sprintf('%05.2f',curr.obs.window.hours),
                                                           '_h_imputation_',
                                                           curr.imp,
                                                           '_full_matrix.rds'))
    curr.obs.window.keys[[curr.imp]] <- read.csv(paste0('../features/03_formatted_predictor_matrices/full_matrices/',
                                                        sprintf('%05.2f',curr.obs.window.hours),
                                                        '_h_imputation_',
                                                        curr.imp,
                                                        '_keys.csv'))
  }
  # Iterate through unique combinations of repeats and folds
  for (curr.combo.idx in 1:nrow(unique.repeat.fold.combos)){
    # Status update on current repeat/fold combination
    print(paste("Repeat/Fold combination no.",curr.combo.idx,"out of",nrow(unique.repeat.fold.combos), "started."))
    # Extract repeat, fold, and imputation information based on current combination
    curr.repeat <- unique.repeat.fold.combos$Repeat[curr.combo.idx]
    curr.fold <- unique.repeat.fold.combos$Fold[curr.combo.idx]
    curr.imputation <- unique.repeat.fold.combos$Imputation[curr.combo.idx]
    # Create directory base for saving all outputs in current repeat/fold combination
    curr.combo.directory <- paste0('../features/03_formatted_predictor_matrices/predictor_matrices/',
                                   sprintf('%05.2f',curr.obs.window.hours),
                                   '_h_obs_window/repeat',
                                   sprintf('%02.f',curr.repeat),
                                   '/fold',
                                   sprintf('%02.f',curr.fold))
    # Create directory to store current repeat/fold combo information
    dir.create(curr.combo.directory,showWarnings = F,recursive = T)
    # Separate current cross-validation observation information for training and testing 
    curr.train.splits <- curr.validation.splits %>%
      filter(Fold == curr.fold,
             Repeat == curr.repeat,
             Split == 'Train')
    curr.test.splits <- curr.validation.splits %>%
      filter(Fold == curr.fold,
             Repeat == curr.repeat,
             Split == 'Test')
    
    # Extract full matrix based on the drawn imputation of the current repeat/fold combination
    curr.full.matrix <- curr.obs.window.matrices[[curr.imputation]]
    curr.full.matrix.key <- curr.obs.window.keys[[curr.imputation]] %>%
      mutate(MatrixRowIdx = 1:nrow(.))
    
    # Separate training and testing keys of current repeat/fold combination and encode outcome (GCS and GOSE12m) information
    curr.train.filt.matrix.key <- inner_join(curr.full.matrix.key,curr.train.splits,
                                             by = c("UPI", 
                                                    "HoursFromICUAdmission", 
                                                    "TimeOfDay", 
                                                    "GCST", 
                                                    "GCSm", 
                                                    "GCSe", 
                                                    "GCSv", 
                                                    "CoincidesWithAccelRecording")) %>%
      mutate(GCSm.gt.1 = as.integer(GCSm > 1),
             GCSm.gt.2 = as.integer(GCSm > 2),
             GCSm.gt.3 = as.integer(GCSm > 3),
             GCSm.gt.4 = as.integer(GCSm > 4),
             GCSm.gt.5 = as.integer(GCSm > 5)) %>%
      left_join(patient.outcomes %>% select(UPI,GOSEDischarge,GOSE12Months),by='UPI') %>%
      mutate(GOSE.gt.1 = as.integer(GOSEDischarge > 1),
             GOSE.gt.2 = as.integer(GOSEDischarge > 2),
             GOSE.gt.3 = as.integer(GOSEDischarge > 3),
             GOSE.gt.4 = as.integer(GOSEDischarge > 4),
             GOSE.gt.5 = as.integer(GOSEDischarge > 5)) %>%
      mutate(GOSE12m.gt.1 = as.integer(GOSE12Months > 1),
             GOSE12m.gt.2 = as.integer(GOSE12Months > 2),
             GOSE12m.gt.3 = as.integer(GOSE12Months > 3),
             GOSE12m.gt.4 = as.integer(GOSE12Months > 4),
             GOSE12m.gt.5 = as.integer(GOSE12Months > 5),
             GOSE12m.gt.6 = as.integer(GOSE12Months > 6),
             GOSE12m.gt.7 = as.integer(GOSE12Months > 7))
    
    curr.test.filt.matrix.key <- inner_join(curr.full.matrix.key,curr.test.splits,
                                            by = c("UPI", 
                                                   "HoursFromICUAdmission", 
                                                   "TimeOfDay", 
                                                   "GCST", 
                                                   "GCSm", 
                                                   "GCSe", 
                                                   "GCSv", 
                                                   "CoincidesWithAccelRecording")) %>%
      mutate(GCSm.gt.1 = as.integer(GCSm > 1),
             GCSm.gt.2 = as.integer(GCSm > 2),
             GCSm.gt.3 = as.integer(GCSm > 3),
             GCSm.gt.4 = as.integer(GCSm > 4),
             GCSm.gt.5 = as.integer(GCSm > 5)) %>%
      left_join(patient.outcomes %>% select(UPI,GOSEDischarge,GOSE12Months),by='UPI') %>%
      mutate(GOSE.gt.1 = as.integer(GOSEDischarge > 1),
             GOSE.gt.2 = as.integer(GOSEDischarge > 2),
             GOSE.gt.3 = as.integer(GOSEDischarge > 3),
             GOSE.gt.4 = as.integer(GOSEDischarge > 4),
             GOSE.gt.5 = as.integer(GOSEDischarge > 5)) %>%
      mutate(GOSE12m.gt.1 = as.integer(GOSE12Months > 1),
             GOSE12m.gt.2 = as.integer(GOSE12Months > 2),
             GOSE12m.gt.3 = as.integer(GOSE12Months > 3),
             GOSE12m.gt.4 = as.integer(GOSE12Months > 4),
             GOSE12m.gt.5 = as.integer(GOSE12Months > 5),
             GOSE12m.gt.6 = as.integer(GOSE12Months > 6),
             GOSE12m.gt.7 = as.integer(GOSE12Months > 7))
    
    # Save current training and testing keys
    write.csv(curr.train.filt.matrix.key,file.path(curr.combo.directory,'GOSE12m_train_key.csv'),row.names = F)
    write.csv(curr.test.filt.matrix.key,file.path(curr.combo.directory,'GOSE12m_test_key.csv'),row.names = F)
    
    # Separate full training and testing matrices and scale columns (each unique combination of feature and sensor type) based on training set
    curr.train.matrix <- abs(curr.full.matrix[curr.train.filt.matrix.key$MatrixRowIdx,])
    curr.test.matrix <- abs(curr.full.matrix[curr.test.filt.matrix.key$MatrixRowIdx,])
    
    for (curr.sf.combo.idx in 1:nrow(sensor.feature.combos)){
      # Get current sensor/feature combo
      curr.sf.combo <- sensor.feature.combos$formatted.feature.name[curr.sf.combo.idx]
      # Find indices of column names with current combination
      curr.sf.col.idx <- which(grepl(curr.sf.combo, colnames(curr.train.matrix), fixed = TRUE))
      # Get mean and standard deviation information from training matrix
      curr.sf.mean <- mean(curr.train.matrix[,curr.sf.col.idx],na.rm = T)
      curr.sf.std <- sd(curr.train.matrix[,curr.sf.col.idx],na.rm = T)
      # Transform training and testing sets accordingly
      curr.train.matrix[,curr.sf.col.idx] <- (curr.train.matrix[,curr.sf.col.idx] - curr.sf.mean)/curr.sf.std
      curr.test.matrix[,curr.sf.col.idx] <- (curr.test.matrix[,curr.sf.col.idx] - curr.sf.mean)/curr.sf.std
    }
    
    # Save full testing and training matrices to appropriate directory
    saveRDS(curr.train.matrix,file.path(curr.combo.directory,'GOSE12m_full_train_matrix.rds'))
    saveRDS(curr.test.matrix,file.path(curr.combo.directory,'GOSE12m_full_test_matrix.rds'))
    
    ## Full GOSE12m (ordinal response)
    # Train and save LOL projection
    curr.full.GOSE12m.LOL <- lol.project.lol(curr.train.matrix,curr.train.filt.matrix.key$GOSE12Months,r = 20)
    saveRDS(curr.full.GOSE12m.LOL,file.path(curr.combo.directory,'LOL_full_GOSE12m.rds'))
    # Separate LOL-embedded training and testing matrices and save into the directory
    curr.full.GOSE12m.LOL.train.matrix <- curr.full.GOSE12m.LOL$Xr
    saveRDS(curr.full.GOSE12m.LOL.train.matrix,file.path(curr.combo.directory,'LOL_full_GOSE12m_train_matrix.rds'))
    curr.full.GOSE12m.LOL.test.matrix <- curr.test.matrix %*% curr.full.GOSE12m.LOL$A
    saveRDS(curr.full.GOSE12m.LOL.test.matrix,file.path(curr.combo.directory,'LOL_full_GOSE12m_test_matrix.rds'))
    
    ## Thresholded GOSE12m (dichotomous response)
    # Extract GOSE12m thresholded response variable names
    threshold.GOSE12m.names <- names(curr.train.filt.matrix.key)[startsWith(names(curr.train.filt.matrix.key),'GOSE12m.gt.')]
    # Iterate through threshold outcomes and perform LOL on each outcome
    for (curr.threshold.name in threshold.GOSE12m.names){
      # Skip current threshold if only one outcome case available 
      if (length(unique(curr.train.filt.matrix.key[[curr.threshold.name]])) == 1){
        print(paste0('Only one outcome for ',curr.threshold.name,' in repeat/fold combination no. ',curr.combo.idx))
        next
      }
      # Train and save LOL projection
      curr.thresh.GOSE12m.LOL <- lol.project.lol(curr.train.matrix,curr.train.filt.matrix.key[[curr.threshold.name]],r = 20)
      saveRDS(curr.thresh.GOSE12m.LOL,file.path(curr.combo.directory,paste0('LOL_thresh_',curr.threshold.name,'.rds')))
      # Separate LOL-embedded training and testing matrices and save into the directory
      curr.thresh.GOSE12m.LOL.train.matrix <- curr.thresh.GOSE12m.LOL$Xr
      saveRDS(curr.thresh.GOSE12m.LOL.train.matrix,file.path(curr.combo.directory,paste0('LOL_thresh_',curr.threshold.name,'_train_matrix.rds')))
      curr.thresh.GOSE12m.LOL.test.matrix <- curr.test.matrix %*% curr.thresh.GOSE12m.LOL$A
      saveRDS(curr.thresh.GOSE12m.LOL.test.matrix,file.path(curr.combo.directory,paste0('LOL_thresh_',curr.threshold.name,'_test_matrix.rds')))
    }
    
    # Status update on current repeat/fold combination
    print(paste("Repeat/Fold combination no.",curr.combo.idx,"out of",nrow(unique.repeat.fold.combos), "completed."))
  }
  # Status update on current observation window
  print(paste("Observation window",sprintf('%05.2f',curr.obs.window.hours),"completed."))
}