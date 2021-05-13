#### Master Script 10: Create manuscript tables and perform statistical analyses ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# University of Cambridge
# Johns Hopkins University
# email address: sb2406@cam.ac.uk
#
### Contents:
# I. Initialization
# II. Table 1: Study population characteristics
# III. Table 2: Classification performance metrics of optimally discriminating models
# IV. Supplementary Table 1: Count distributions of GCSm scores per observation window
# V. Supplementary Table 2: Discrimination of threshold-level GCSm detection models per observation window
# VI. Supplementary Table 3: Count distributions of GOSE scores at hospital discharge per observation window
# VII. Supplementary Table 4: Discrimination of threshold-level GOSE at hospital discharge prediction models per observation window
# VIII. Supplementary Table 5: Count distributions of GOSE scores at 12 months post discharge per observation window
# IX. Supplementary Table 6: Discrimination of threshold-level GOSE at 12 months post discharge prediction models per observation window
# X. Supplementary Table 7: Percentages of missing accelerometry data per sensor and recording duration of each study participant
# XI. Metrics for Figure 5: Feature significance matrices of optimally discriminating motor function detection and functional outcome prediction models
# XII. Metrics for Supplementary Figure 2: Correlation matrices of extracted motion features across different sensor placements
# XIII. Metrics for Supplementary Figure 4: Mean motion feature trajectories in the six hours preceding GCSm evaluation, stratified by GCSm scores and bilateral sensor placement
# XIV. Miscellaneous statistics for manuscript

### I. Initialization
## Load necessary libraries
library(tidyverse)
library(tidymodels)
library(readxl)
library(doParallel)
library(foreach)
library(corrr)
library(R.matlab)

## Compile all motion features into one file for subsequent analysis
# Load all motion features and compile into one dataframe (if not already loaded)
if (!exists("all.motion.features")) {
  feature.files <- list.files('../features',pattern=glob2rx("features_*.csv"),full.names = T)
  all.motion.features <- as.data.frame(matrix(nrow = 0, ncol = 12))
  for (curr.feature.file in feature.files) {
    curr.upi.df <- read.csv(curr.feature.file) 
    names(all.motion.features) <- names(curr.upi.df)
    all.motion.features <- rbind(all.motion.features,curr.upi.df)
    print(paste('Feature file',which(feature.files == curr.feature.file),'out of',length(feature.files),'appended'))
  }
}

# Save compiled dataframe of all features
write.csv(all.motion.features,'../features/all_features.csv',row.names = F)

## Compile bootstrapping metric results across tasks, thresholds, and observation windows
# Acquire observation windows and list of threshold-specific metric files for GCSm detection
obs.window.dirs <- list.files('../results/GCSm_threshold_prediction','*_h_obs_window',include.dirs = T,full.names = T)
metrics.GCSm.files <- list.files(obs.window.dirs,'*_compiled_metrics.csv',include.dirs = F,full.names = T)

# Create a dataframe to store metrics across observation windows and thresholds
compiled.GCSm.metrics.df <- as.data.frame(matrix(ncol=7,nrow=0))

# Iterate through threshold-specific metric files
for (curr.metric.file in metrics.GCSm.files){
  # Load current compiled metric dataframe
  curr.metric.df <- read.csv(curr.metric.file) 
  # Check that compiled metric dataframe has correct number of columns
  if (ncol(curr.metric.df) == 5){
    # Append suitable dataframe summary statistics to compiled metric dataframe
    curr.metric.df <- curr.metric.df %>%
      group_by(Threshold,ObsWindow,Metrics) %>%
      summarise(meanValue = mean(Values,na.rm = T),
                medianValues = quantile(Values,.5,na.rm = T),
                lowerValues = quantile(Values,.025,na.rm = T),
                upperValues = quantile(Values,.975,na.rm = T))
    compiled.GCSm.metrics.df <- rbind(compiled.GCSm.metrics.df,curr.metric.df)
  }
}

# Save compiled GCSm metrics as a dataframe
write.csv(compiled.GCSm.metrics.df,'../results/GCSm_threshold_prediction/compiled_metrics.csv',row.names = F)

# Acquire observation windows and list of threshold-specific metric files for GOSE (discharge) prediction
obs.window.dirs <- list.files('../results/GOSE_threshold_prediction','*_h_obs_window',include.dirs = T,full.names = T)
metrics.GOSE.files <- list.files(obs.window.dirs,'*_compiled_metrics.csv',include.dirs = F,full.names = T)

# Create a dataframe to store metrics across observation windows and thresholds
compiled.GOSE.metrics.df <- as.data.frame(matrix(ncol=7,nrow=0))

# Iterate through threshold-specific metric files
for (curr.metric.file in metrics.GOSE.files){
  # Load current compiled metric dataframe
  curr.metric.df <- read.csv(curr.metric.file) 
  # Check that compiled metric dataframe has correct number of columns
  if (ncol(curr.metric.df) == 5){
    # Append suitable dataframe summary statistics to compiled metric dataframe
    curr.metric.df <- curr.metric.df %>%
      group_by(Threshold,ObsWindow,Metrics) %>%
      summarise(meanValue = mean(Values,na.rm = T),
                medianValues = quantile(Values,.5,na.rm = T),
                lowerValues = quantile(Values,.025,na.rm = T),
                upperValues = quantile(Values,.975,na.rm = T))
    compiled.GOSE.metrics.df <- rbind(compiled.GOSE.metrics.df,curr.metric.df)
  }
}

# Save compiled GOSE metrics as a dataframe
write.csv(compiled.GOSE.metrics.df,'../results/GOSE_threshold_prediction/compiled_metrics.csv',row.names = F)

# Acquire observation windows and list of threshold-specific metric files for GOSE (12 months) prediction
obs.window.dirs <- list.files('../results/GOSE12m_threshold_prediction','*_h_obs_window',include.dirs = T,full.names = T)
metrics.GOSE12m.files <- list.files(obs.window.dirs,'*_compiled_metrics.csv',include.dirs = F,full.names = T)

# Create a dataframe to store metrics across observation windows and thresholds
compiled.GOSE12m.metrics.df <- as.data.frame(matrix(ncol=7,nrow=0))

# Iterate through threshold-specific metric files
for (curr.metric.file in metrics.GOSE12m.files){
  # Load current compiled metric dataframe
  curr.metric.df <- read.csv(curr.metric.file) 
  # Check that compiled metric dataframe has correct number of columns
  if (ncol(curr.metric.df) == 5){
    # Append suitable dataframe summary statistics to compiled metric dataframe
    curr.metric.df <- curr.metric.df %>%
      group_by(Threshold,ObsWindow,Metrics) %>%
      summarise(meanValue = mean(Values,na.rm = T),
                medianValues = quantile(Values,.5,na.rm = T),
                lowerValues = quantile(Values,.025,na.rm = T),
                upperValues = quantile(Values,.975,na.rm = T))
    compiled.GOSE12m.metrics.df <- rbind(compiled.GOSE12m.metrics.df,curr.metric.df)
  }
}

# Save compiled GOSE12m metrics as a dataframe
write.csv(compiled.GOSE12m.metrics.df,'../results/GOSE12m_threshold_prediction/compiled_metrics.csv',row.names = F)

## Compile ROC curves across tasks, thresholds, and observation windows 
# Acquire observation windows and list of threshold-specific ROC files for GCSm detection
obs.window.dirs <- list.files('../results/GCSm_threshold_prediction','*_h_obs_window',include.dirs = T,full.names = T)
ROC.GCSm.files <- list.files(obs.window.dirs,'*_compiled_ROC.csv',include.dirs = F,full.names = T)

# Create a dataframe to store ROC across observation windows and thresholds
compiled.GCSm.ROC.df <- as.data.frame(matrix(ncol=7,nrow=0))

# Iterate through threshold-specific ROC files
for (curr.ROC.file in ROC.GCSm.files){
  # Load current compiled ROC dataframe
  curr.ROC.df <- read.csv(curr.ROC.file) 
  # Check that compiled ROC dataframe has correct number of columns
  if (ncol(curr.ROC.df) == 5){
    # Append suitable dataframe summary statistics to compiled ROC dataframe
    curr.ROC.df <- curr.ROC.df %>%
      group_by(Threshold,ObsWindow,FPR) %>%
      summarise(meanTPR = mean(TPR,na.rm = T),
                medianTPR = quantile(TPR,.5,na.rm = T),
                lowerTPR = quantile(TPR,.025,na.rm = T),
                upperTPR = quantile(TPR,.975,na.rm = T))
    compiled.GCSm.ROC.df <- rbind(compiled.GCSm.ROC.df,curr.ROC.df)
  }
}

# Save compiled GCSm ROC as a dataframe
write.csv(compiled.GCSm.ROC.df,'../results/GCSm_threshold_prediction/compiled_ROC.csv',row.names = F)

# Acquire observation windows and list of threshold-specific ROC files for GOSE (discharge) prediction
obs.window.dirs <- list.files('../results/GOSE_threshold_prediction','*_h_obs_window',include.dirs = T,full.names = T)
ROC.GOSE.files <- list.files(obs.window.dirs,'*_compiled_ROC.csv',include.dirs = F,full.names = T)

# Create a dataframe to store ROC across observation windows and thresholds
compiled.GOSE.ROC.df <- as.data.frame(matrix(ncol=7,nrow=0))

# Iterate through threshold-specific ROC files
for (curr.ROC.file in ROC.GOSE.files){
  # Load current compiled ROC dataframe
  curr.ROC.df <- read.csv(curr.ROC.file) 
  # Check that compiled ROC dataframe has correct number of columns
  if (ncol(curr.ROC.df) == 5){
    # Append suitable dataframe summary statistics to compiled ROC dataframe
    curr.ROC.df <- curr.ROC.df %>%
      group_by(Threshold,ObsWindow,FPR) %>%
      summarise(meanTPR = mean(TPR,na.rm = T),
                medianTPR = quantile(TPR,.5,na.rm = T),
                lowerTPR = quantile(TPR,.025,na.rm = T),
                upperTPR = quantile(TPR,.975,na.rm = T))
    compiled.GOSE.ROC.df <- rbind(compiled.GOSE.ROC.df,curr.ROC.df)
  }
}

# Save compiled GOSE ROC as a dataframe
write.csv(compiled.GOSE.ROC.df,'../results/GOSE_threshold_prediction/compiled_ROC.csv',row.names = F)

# Acquire observation windows and list of threshold-specific ROC files for GOSE (12 months) prediction
obs.window.dirs <- list.files('../results/GOSE12m_threshold_prediction','*_h_obs_window',include.dirs = T,full.names = T)
ROC.GOSE12m.files <- list.files(obs.window.dirs,'*_compiled_ROC.csv',include.dirs = F,full.names = T)

# Create a dataframe to store ROC across observation windows and thresholds
compiled.GOSE12m.ROC.df <- as.data.frame(matrix(ncol=7,nrow=0))

# Iterate through threshold-specific ROC files
for (curr.ROC.file in ROC.GOSE12m.files){
  # Load current compiled ROC dataframe
  curr.ROC.df <- read.csv(curr.ROC.file) 
  # Check that compiled ROC dataframe has correct number of columns
  if (ncol(curr.ROC.df) == 5){
    # Append suitable dataframe summary statistics to compiled ROC dataframe
    curr.ROC.df <- curr.ROC.df %>%
      group_by(Threshold,ObsWindow,FPR) %>%
      summarise(meanTPR = mean(TPR,na.rm = T),
                medianTPR = quantile(TPR,.5,na.rm = T),
                lowerTPR = quantile(TPR,.025,na.rm = T),
                upperTPR = quantile(TPR,.975,na.rm = T))
    compiled.GOSE12m.ROC.df <- rbind(compiled.GOSE12m.ROC.df,curr.ROC.df)
  }
}

# Save compiled GOSE12m ROC as a dataframe
write.csv(compiled.GOSE12m.ROC.df,'../results/GOSE12m_threshold_prediction/compiled_ROC.csv',row.names = F)

### II. Table 1: Study population characteristics

## Load clinical dataframes
# baseline characteristics:
patient.baseline.characteristics <- read.csv('../clinical_data/patient_baseline_characteristics.csv')

# GCS evaluation scores:
gcs.data <- read.csv('../clinical_data/neurological_assessments.csv')

# clinical outcomes:
patient.outcomes <- read.csv('../clinical_data/patient_outcomes.csv')

# temporal information:
patient.temporal.info <- read.csv('../clinical_data/patient_temporal_info.csv')

## Calculate appropriate baseline characteristic information
# Age:
age.info <- patient.baseline.characteristics$Age
age.quantiles <- quantile(age.info,c(.25,.5,.75))
print(paste('Age:',sprintf('%01.f (%01.f–%01.f)',age.quantiles[2],age.quantiles[1],age.quantiles[3])))

# Sex:
sex.info <- patient.baseline.characteristics$Sex
sex.count.distribution <- table(sex.info)
print(paste('M/F:',sprintf('%01.f/%01.f',sex.count.distribution[2],sex.count.distribution[1])))

# Types of SBI:
SBI.types <- data.frame(ICH = sum(patient.baseline.characteristics$ICH),
                        SDHorEDH = sum(patient.baseline.characteristics$SDHOrEDH),
                        SAH = sum(patient.baseline.characteristics$SAH),
                        CVA = sum(patient.baseline.characteristics$CVA),
                        BTorST = sum(patient.baseline.characteristics$BrainTumorOrLesion),
                        TBI = sum(patient.baseline.characteristics$TBI)) %>%
  pivot_longer(cols = everything(),names_to = 'TypeSBI',values_to='count') %>%
  mutate(percentage = 100*count/length(unique(patient.baseline.characteristics$UPI))) %>%
  mutate(Formatted = sprintf('%1.f (%0.2f)',count,percentage))

## Calculate GCSm information
# Join patient temporal data with GCSm evaluations and filter out evalations during ICU stay
icu.gcs.data <- gcs.data %>%
  drop_na(GCSm) %>%
  left_join(patient.temporal.info,by='UPI') %>%
  mutate(DuringICU = HoursFromICUAdmission <= ((DaysInICU+1)*24)) %>%
  filter(DuringICU == T)

# Calculate number of GCS observations during ICU per patient
obs.per.patient <- icu.gcs.data %>%
  group_by(UPI) %>%
  summarise(total.no.evals = n(), evals.per.day = n()/mean(DaysInICU))

# Calculate number of GCS observations coinciding with accelerometry recording per patient
accel.obs.per.patient <- icu.gcs.data %>%
  filter(CoincidesWithAccelRecording == T) %>%
  group_by(UPI) %>%
  summarise(total.no.evals = n(), evals.per.day = n()/mean(HoursDurationAccelRecording/24))

### III. Table 2: Classification performance metrics of optimally discriminating models

## Initialize parameters for parallelized bootstrapping
# Set number of boostrap resamples
NUM.BOOTSTRAPS <- 1000

# Set number of cores to use in parallel
NUM.CORES <- 10

# Initialize `doParallel` cluster
registerDoParallel(cores = NUM.CORES)

## Determine optimally discriminating observation windows for each task-threhsold combination
# GCSm detection:
opt.GCSm.signficant.AUC.df <- read.csv('../results/GCSm_threshold_prediction/compiled_metrics.csv') %>%
  filter(Metrics == 'AUC') %>%
  filter(lowerValues >= 0.50) %>%
  group_by(Threshold) %>%
  summarise(ObsWindow = ObsWindow[which.max(meanValue)],
            meanAUC = max(meanValue),
            medianAUC = medianValues[which.max(meanValue)],
            lowerAUC = lowerValues[which.max(meanValue)],
            upperAUC = upperValues[which.max(meanValue)])

opt.GCSm.nonsignficant.AUC.df <- read.csv('../results/GCSm_threshold_prediction/compiled_metrics.csv') %>%
  filter(Metrics == 'AUC') %>%
  filter(!(Threshold %in% opt.GCSm.signficant.AUC.df$Threshold)) %>%
  group_by(Threshold) %>%
  summarise(ObsWindow = ObsWindow[which.max(meanValue)],
            meanAUC = max(meanValue),
            medianAUC = medianValues[which.max(meanValue)],
            lowerAUC = lowerValues[which.max(meanValue)],
            upperAUC = upperValues[which.max(meanValue)])

opt.GCSm.AUC.df <- rbind(opt.GCSm.signficant.AUC.df,opt.GCSm.nonsignficant.AUC.df)

# GOSE at hospital discharge prediction:
opt.GOSE.signficant.AUC.df <- read.csv('../results/GOSE_threshold_prediction/compiled_metrics.csv') %>%
  filter(Metrics == 'AUC') %>%
  filter(lowerValues >= 0.50) %>%
  group_by(Threshold) %>%
  summarise(ObsWindow = ObsWindow[which.max(meanValue)],
            meanAUC = max(meanValue),
            medianAUC = medianValues[which.max(meanValue)],
            lowerAUC = lowerValues[which.max(meanValue)],
            upperAUC = upperValues[which.max(meanValue)])

opt.GOSE.nonsignficant.AUC.df <- read.csv('../results/GOSE_threshold_prediction/compiled_metrics.csv') %>%
  filter(Metrics == 'AUC') %>%
  filter(!(Threshold %in% opt.GOSE.signficant.AUC.df$Threshold)) %>%
  group_by(Threshold) %>%
  summarise(ObsWindow = ObsWindow[which.max(meanValue)],
            meanAUC = max(meanValue),
            medianAUC = medianValues[which.max(meanValue)],
            lowerAUC = lowerValues[which.max(meanValue)],
            upperAUC = upperValues[which.max(meanValue)])

opt.GOSE.AUC.df <- rbind(opt.GOSE.signficant.AUC.df,opt.GOSE.nonsignficant.AUC.df)

# GOSE at 12 months prediction:
opt.GOSE12m.signficant.AUC.df <- read.csv('../results/GOSE12m_threshold_prediction/compiled_metrics.csv') %>%
  filter(Metrics == 'AUC') %>%
  filter(lowerValues >= 0.50) %>%
  group_by(Threshold) %>%
  summarise(ObsWindow = ObsWindow[which.max(meanValue)],
            meanAUC = max(meanValue),
            medianAUC = medianValues[which.max(meanValue)],
            lowerAUC = lowerValues[which.max(meanValue)],
            upperAUC = upperValues[which.max(meanValue)])

opt.GOSE12m.nonsignficant.AUC.df <- read.csv('../results/GOSE12m_threshold_prediction/compiled_metrics.csv') %>%
  filter(Metrics == 'AUC') %>%
  filter(!(Threshold %in% opt.GOSE12m.signficant.AUC.df$Threshold)) %>%
  group_by(Threshold) %>%
  summarise(ObsWindow = ObsWindow[which.max(meanValue)],
            meanAUC = max(meanValue),
            medianAUC = medianValues[which.max(meanValue)],
            lowerAUC = lowerValues[which.max(meanValue)],
            upperAUC = upperValues[which.max(meanValue)])

opt.GOSE12m.AUC.df <- rbind(opt.GOSE12m.signficant.AUC.df,opt.GOSE12m.nonsignficant.AUC.df)

## Calculate case count distributions at each threshold of each task of the optimally discriminating configurations
# GCSm detection:
compiled.GCSm.threshold.props <- as.data.frame(matrix(ncol = 7, nrow = 0))

for (curr.thresh in opt.GCSm.AUC.df$Threshold){
  # Extract observation window of current optimal configuration
  curr.obs.window <- opt.GCSm.AUC.df$ObsWindow[opt.GCSm.AUC.df$Threshold == curr.thresh]
  
  # Load compiled predictions of current threshold-observation window combination and filter out unique observations
  curr.proportions <-
    read.csv(file.path('../results/GCSm_threshold_prediction',
                       paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
                       paste0(curr.thresh,'_compiled_predictions.csv'))) %>%
    dplyr::select(-c(Prob,ConfigIdx,TargetDim,Repeat,Fold)) %>%
    distinct()
  
  # Paste count distribution of unique observations for Table 2
  print(paste(sprintf('%01.f/%01.f',table(curr.proportions$TrueLabel)[1],table(curr.proportions$TrueLabel)[2]),
              sprintf('(%0.2f)',(table(curr.proportions$TrueLabel)/nrow(curr.proportions))[2])))
}

# GOSE at hospital discharge prediction:
compiled.GOSE.threshold.props <- as.data.frame(matrix(ncol = 7, nrow = 0))

for (curr.thresh in opt.GOSE.AUC.df$Threshold){
  # Extract observation window of current optimal configuration
  curr.obs.window <- opt.GOSE.AUC.df$ObsWindow[opt.GOSE.AUC.df$Threshold == curr.thresh]
  
  # Load compiled predictions of current threshold-observation window combination and filter out unique observations
  curr.proportions <-
    read.csv(file.path('../results/GOSE_threshold_prediction',
                       paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
                       paste0(curr.thresh,'_compiled_predictions.csv'))) %>%
    dplyr::select(-c(Prob,ConfigIdx,TargetDim,Repeat,Fold)) %>%
    distinct()
  
  # Paste count distribution of unique observations for Table 2
  print(paste(sprintf('%01.f/%01.f',table(curr.proportions$TrueLabel)[1],table(curr.proportions$TrueLabel)[2]),
              sprintf('(%0.2f)',(table(curr.proportions$TrueLabel)/nrow(curr.proportions))[2])))
}

# GOSE at 12 months post discharge prediction:
compiled.GOSE12m.threshold.props <- as.data.frame(matrix(ncol = 7, nrow = 0))

for (curr.thresh in opt.GOSE12m.AUC.df$Threshold){
  # Extract observation window of current optimal configuration
  curr.obs.window <- opt.GOSE12m.AUC.df$ObsWindow[opt.GOSE12m.AUC.df$Threshold == curr.thresh]
  
  # Load compiled predictions of current threshold-observation window combination and filter out unique observations
  curr.proportions <-
    read.csv(file.path('../results/GOSE12m_threshold_prediction',
                       paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
                       paste0(curr.thresh,'_compiled_predictions.csv'))) %>%
    dplyr::select(-c(Prob,ConfigIdx,TargetDim,Repeat,Fold)) %>%
    distinct()
  
  # Paste count distribution of unique observations for Table 2
  print(paste(sprintf('%01.f/%01.f',table(curr.proportions$TrueLabel)[1],table(curr.proportions$TrueLabel)[2]),
              sprintf('(%0.2f)',(table(curr.proportions$TrueLabel)/nrow(curr.proportions))[2])))
}

## Calculate classification accuracy for optimally discriminating model configurations and perform bias-corrected bootstrapping to calculate 95% confidence intervals
# GCSm detection:
compiled.GCSm.threshold.metrics <- as.data.frame(matrix(ncol = 7, nrow = 0))

for (curr.thresh in opt.GCSm.AUC.df$Threshold){
  # Extract observation window of current optimal configuration
  curr.obs.window <- opt.GCSm.AUC.df$ObsWindow[opt.GCSm.AUC.df$Threshold == curr.thresh]
  
  # Load compiled predictions of current threshold-observation window combination
  curr.predictions <-
    read.csv(file.path('../results/GCSm_threshold_prediction',
                       paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
                       paste0(curr.thresh,'_compiled_predictions.csv')))
  
  # Identify unique UPIs available
  unique.UPIs <- unique(curr.predictions$UPI)
  
  # Once bootstrap samples have been confirmed, begin parallel bootstrapping
  compiled.GCSm.threshold.accuracy <- foreach(icount(NUM.BOOTSTRAPS), .combine=rbind) %dopar%{
    # Keep drawing sample until both cases are present in both in- and out-sample cases
    fail.condition <- TRUE
    while(fail.condition){
      curr.UPI.resample <- sample(unique.UPIs,length(unique.UPIs),replace = T)
      # Divide in-sample and out-sample predictions
      curr.in.sample.preds <- curr.predictions %>% filter(UPI %in% sort(unique(curr.UPI.resample)))
      curr.out.sample.preds <- curr.predictions %>% filter(UPI %in% sort(unique.UPIs[! unique.UPIs %in% sort(unique(curr.UPI.resample))]))
      # If the necessary condition is met, we may break out of the while loop
      if ((length(unique(curr.in.sample.preds$TrueLabel)) == 2) & 
          (length(unique(curr.out.sample.preds$TrueLabel)) == 2)){
        fail.condition <- FALSE
      }
    }
    
    # Determine optimal configuration in current resample for calibration based on in-sample accuracy
    opt.config <- curr.in.sample.preds %>%
      mutate(CorrectClassif = as.integer(TrueLabel == as.integer(Prob > .5))) %>%
      group_by(ConfigIdx) %>%
      summarise(Accuracy = sum(CorrectClassif)/n()) %>%
      top_n(1,Accuracy)
    
    # Calculate accuracy for current optimal configuration
    curr.config.out.sample.preds <- curr.out.sample.preds %>% 
      filter(ConfigIdx == opt.config$ConfigIdx[1]) %>%
      mutate(CorrectClassif = as.integer(TrueLabel == as.integer(Prob > .5)))
    curr.Accuracy <- sum(curr.config.out.sample.preds$CorrectClassif)/nrow(curr.config.out.sample.preds)
    
    # Return dataframe row of compiled information
    data.frame(
      Threshold = curr.thresh,
      ObsWindow = curr.obs.window,
      Accuracy = curr.Accuracy
    )
  }
  
  # Derive 95% confidence interval of compiled accuracy
  summ.GCSm.threshold.accuracy <- compiled.GCSm.threshold.accuracy %>%
    group_by(Threshold,ObsWindow) %>%
    summarise(Metrics = 'accuracy',
              meanValue = mean(Accuracy,na.rm = T),
              medianValue = median(Accuracy,na.rm = T),
              lowerValue = quantile(Accuracy,.025,na.rm = T,),
              upperValue = quantile(Accuracy,.975,na.rm = T,))
  
  # Append to compiled metric dataframe
  compiled.GCSm.threshold.metrics <- rbind(compiled.GCSm.threshold.metrics,summ.GCSm.threshold.accuracy)
  
  # Status update on completion of threshold
  print(paste(curr.thresh,'complete.'))
}
write.csv(compiled.GCSm.threshold.metrics,'../results/GCSm_threshold_prediction/accuracy.csv',row.names = F)

# GOSE prediction at hospital discharge:
compiled.GOSE.threshold.metrics <- as.data.frame(matrix(ncol = 7, nrow = 0))

for (curr.thresh in opt.GOSE.AUC.df$Threshold){
  # Extract observation window of current optimal configuration
  curr.obs.window <- opt.GOSE.AUC.df$ObsWindow[opt.GOSE.AUC.df$Threshold == curr.thresh]
  
  # Load compiled predictions of current threshold-observation window combination
  curr.predictions <-
    read.csv(file.path('../results/GOSE_threshold_prediction',
                       paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
                       paste0(curr.thresh,'_compiled_predictions.csv')))
  
  # Identify unique UPIs available
  unique.UPIs <- unique(curr.predictions$UPI)
  
  # Once bootstrap samples have been confirmed, begin parallel bootstrapping
  compiled.GOSE.threshold.accuracy <- foreach(icount(NUM.BOOTSTRAPS), .combine=rbind) %dopar%{
    # Keep drawing sample until both cases are present in both in- and out-sample cases
    fail.condition <- TRUE
    while(fail.condition){
      curr.UPI.resample <- sample(unique.UPIs,length(unique.UPIs),replace = T)
      # Divide in-sample and out-sample predictions
      curr.in.sample.preds <- curr.predictions %>% filter(UPI %in% sort(unique(curr.UPI.resample)))
      curr.out.sample.preds <- curr.predictions %>% filter(UPI %in% sort(unique.UPIs[! unique.UPIs %in% sort(unique(curr.UPI.resample))]))
      # If the necessary condition is met, we may break out of the while loop
      if ((length(unique(curr.in.sample.preds$TrueLabel)) == 2) & 
          (length(unique(curr.out.sample.preds$TrueLabel)) == 2)){
        fail.condition <- FALSE
      }
    }
    
    # Determine optimal configuration in current resample for calibration based on in-sample accuracy
    opt.config <- curr.in.sample.preds %>%
      mutate(CorrectClassif = as.integer(TrueLabel == as.integer(Prob > .5))) %>%
      group_by(ConfigIdx) %>%
      summarise(Accuracy = sum(CorrectClassif)/n()) %>%
      top_n(1,Accuracy)
    
    # Calculate accuracy for current optimal configuration
    curr.config.out.sample.preds <- curr.out.sample.preds %>% 
      filter(ConfigIdx == opt.config$ConfigIdx[1]) %>%
      mutate(CorrectClassif = as.integer(TrueLabel == as.integer(Prob > .5)))
    curr.Accuracy <- sum(curr.config.out.sample.preds$CorrectClassif)/nrow(curr.config.out.sample.preds)
    
    # Return dataframe row of compiled information
    data.frame(
      Threshold = curr.thresh,
      ObsWindow = curr.obs.window,
      Accuracy = curr.Accuracy
    )
  }
  
  # Derive 95% of compiled accuracy values as well as optimal configuration index
  summ.GOSE.threshold.accuracy <- compiled.GOSE.threshold.accuracy %>%
    group_by(Threshold,ObsWindow) %>%
    summarise(Metrics = 'accuracy',
              meanValue = mean(Accuracy,na.rm = T),
              medianValue = median(Accuracy,na.rm = T),
              lowerValue = quantile(Accuracy,.025,na.rm = T,),
              upperValue = quantile(Accuracy,.975,na.rm = T,))
  
  # Append to compiled metric dataframe
  compiled.GOSE.threshold.metrics <- rbind(compiled.GOSE.threshold.metrics,summ.GOSE.threshold.accuracy)
  
  # Status update on completion of threshold
  print(paste(curr.thresh,'complete.'))
}
write.csv(compiled.GOSE.threshold.metrics,'../results/GOSE_threshold_prediction/accuracy.csv',row.names = F)

# GOSE prediction at 12 months:
compiled.GOSE12m.threshold.metrics <- as.data.frame(matrix(ncol = 7, nrow = 0))

for (curr.thresh in opt.GOSE12m.AUC.df$Threshold){
  # Extract observation window of current optimal configuration
  curr.obs.window <- opt.GOSE12m.AUC.df$ObsWindow[opt.GOSE12m.AUC.df$Threshold == curr.thresh]
  
  # Load compiled predictions of current threshold-observation window combination
  curr.predictions <-
    read.csv(file.path('../results/GOSE12m_threshold_prediction',
                       paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
                       paste0(curr.thresh,'_compiled_predictions.csv')))
  
  # Identify unique UPIs available
  unique.UPIs <- unique(curr.predictions$UPI)
  
  # Once bootstrap samples have been confirmed, begin parallel bootstrapping
  compiled.GOSE12m.threshold.accuracy <- foreach(icount(NUM.BOOTSTRAPS), .combine=rbind) %dopar%{
    # Keep drawing sample until both cases are present in both in- and out-sample cases
    fail.condition <- TRUE
    while(fail.condition){
      curr.UPI.resample <- sample(unique.UPIs,length(unique.UPIs),replace = T)
      # Divide in-sample and out-sample predictions
      curr.in.sample.preds <- curr.predictions %>% filter(UPI %in% sort(unique(curr.UPI.resample)))
      curr.out.sample.preds <- curr.predictions %>% filter(UPI %in% sort(unique.UPIs[! unique.UPIs %in% sort(unique(curr.UPI.resample))]))
      # If the necessary condition is met, we may break out of the while loop
      if ((length(unique(curr.in.sample.preds$TrueLabel)) == 2) & 
          (length(unique(curr.out.sample.preds$TrueLabel)) == 2)){
        fail.condition <- FALSE
      }
    }
    
    # Determine optimal configuration in current resample for calibration based on in-sample accuracy
    opt.config <- curr.in.sample.preds %>%
      mutate(CorrectClassif = as.integer(TrueLabel == as.integer(Prob > .5))) %>%
      group_by(ConfigIdx) %>%
      summarise(Accuracy = sum(CorrectClassif)/n()) %>%
      top_n(1,Accuracy)
    
    # Calculate accuracy for current optimal configuration
    curr.config.out.sample.preds <- curr.out.sample.preds %>% 
      filter(ConfigIdx == opt.config$ConfigIdx[1]) %>%
      mutate(CorrectClassif = as.integer(TrueLabel == as.integer(Prob > .5)))
    curr.Accuracy <- sum(curr.config.out.sample.preds$CorrectClassif)/nrow(curr.config.out.sample.preds)
    
    # Return dataframe row of compiled information
    data.frame(
      Threshold = curr.thresh,
      ObsWindow = curr.obs.window,
      Accuracy = curr.Accuracy
    )
  }
  
  # Derive 95% of compiled accuracy values as well as optimal configuration index
  summ.GOSE12m.threshold.accuracy <- compiled.GOSE12m.threshold.accuracy %>%
    group_by(Threshold,ObsWindow) %>%
    summarise(Metrics = 'accuracy',
              meanValue = mean(Accuracy,na.rm = T),
              medianValue = median(Accuracy,na.rm = T),
              lowerValue = quantile(Accuracy,.025,na.rm = T,),
              upperValue = quantile(Accuracy,.975,na.rm = T,))
  
  # Append to compiled metric dataframe
  compiled.GOSE12m.threshold.metrics <- rbind(compiled.GOSE12m.threshold.metrics,summ.GOSE12m.threshold.accuracy)
  
  # Status update on completion of threshold
  print(paste(curr.thresh,'complete.'))
}
write.csv(compiled.GOSE12m.threshold.metrics,'../results/GOSE12m_threshold_prediction/accuracy.csv',row.names = F)

# Compile accuracies across tasks and thresholds and format for manuscript table
compiled.accuracies <- rbind(compiled.GCSm.threshold.metrics %>% mutate(Task = 'GCSm'),
                             compiled.GOSE.threshold.metrics %>% mutate(Task = 'GOSE'),
                             compiled.GOSE12m.threshold.metrics %>% mutate(Task = 'GOSE12m')) %>%
  arrange(Task, Threshold) %>%
  relocate(Task, Threshold) %>%
  mutate(FormattedAccuracy = sprintf('%0.2f (%0.2f–%0.2f)',meanValue,lowerValue,upperValue))

# Stop implicit cluster from parallel processing
stopImplicitCluster()

## Extract other classification metrics of optimally discriminating model configurations
# GCSm detection:
table.GCSm.metrics <- read.csv('../results/GCSm_threshold_prediction/compiled_metrics.csv') %>%
  mutate(formattedValues = sprintf('%0.2f (%0.2f–%0.2f)',meanValue,lowerValues,upperValues)) %>%
  inner_join(opt.GCSm.AUC.df,by = c('Threshold','ObsWindow')) %>%
  dplyr::select(-ends_with('.y'),-c(meanValue,medianValues,lowerValues,upperValues)) %>%
  pivot_wider(id_cols = c(Threshold,ObsWindow),names_from = Metrics.x,values_from = formattedValues.x) %>%
  arrange(Threshold) %>%
  relocate(Threshold,ObsWindow,precision,recall,specificity,f1_score)

# GOSE at hospital discharge prediction:
table.GOSE.metrics <- read.csv('../results/GOSE_threshold_prediction/compiled_metrics.csv') %>%
  mutate(formattedValues = sprintf('%0.2f (%0.2f–%0.2f)',meanValue,lowerValues,upperValues)) %>%
  inner_join(opt.GOSE.AUC.df,by = c('Threshold','ObsWindow')) %>%
  dplyr::select(-ends_with('.y'),-c(meanValue,medianValues,lowerValues,upperValues)) %>%
  pivot_wider(id_cols = c(Threshold,ObsWindow),names_from = Metrics.x,values_from = formattedValues.x) %>%
  arrange(Threshold) %>%
  relocate(Threshold,ObsWindow,precision,recall,specificity,f1_score)

# GOSE at 12 months prediction:
table.GOSE12m.metrics <- read.csv('../results/GOSE12m_threshold_prediction/compiled_metrics.csv') %>%
  mutate(formattedValues = sprintf('%0.2f (%0.2f–%0.2f)',meanValue,lowerValues,upperValues)) %>%
  inner_join(opt.GOSE12m.AUC.df,by = c('Threshold','ObsWindow')) %>%
  dplyr::select(-ends_with('.y'),-c(meanValue,medianValues,lowerValues,upperValues)) %>%
  pivot_wider(id_cols = c(Threshold,ObsWindow),names_from = Metrics.x,values_from = formattedValues.x) %>%
  arrange(Threshold) %>%
  relocate(Threshold,ObsWindow,precision,recall,specificity,f1_score)

### IV. Supplementary Table 1: Count distributions of GCSm scores per observation window

# Load full matrix keys (one imputation is sufficient) to get GCSm scores per observation
full.matrix.keys <- list.files('../features/03_formatted_predictor_matrices/full_matrices',
                               pattern = glob2rx('*_h_imputation_1_keys.csv'),
                               full.names = T)

# Initialize empty dataframe to store count distributions for table
compiled.count.dist.table <- data.frame(matrix(ncol = 3,nrow = 0))

# Iterate through full matrix keys
for (curr.full.matrix.key.file in full.matrix.keys){
  
  # Load current full matrix key
  curr.key <- read.csv(curr.full.matrix.key.file)
  
  # Extract current observation window from key file name
  curr.obs.window <- as.numeric(str_match(curr.full.matrix.key.file, "full_matrices/\\s*(.*?)\\s*_h_imputation_")[,2])
  
  # Calculate current GCSm frequencies from the full matrix key
  GCSm.freqs <- curr.key %>%
    group_by(GCSm) %>%
    tally() %>%
    drop_na() %>%
    mutate(GCSm = paste0('GCSm=',GCSm)) %>%
    rename(Label = GCSm) %>%
    mutate(n = sprintf('%01.f (%0.2f&)',n,(100*n/sum(.$n))))
  
  # Add information about unique patient count and total observation count and append to compiled distribution dataframe
  compiled.count.dist.table <- rbind(compiled.count.dist.table,
                                     rbind(data.frame(Label = 'n_unique_patients',n = length(unique(curr.key$UPI))),
                                           data.frame(Label = 'n_total_observations',n = nrow(curr.key)),
                                           GCSm.freqs) %>%
                                       mutate(ObsWindow = curr.obs.window))
}

# Pivot compiled count distribution wider for manuscript supplementary table
compiled.count.dist.table <- compiled.count.dist.table %>%
  pivot_wider(id_cols = ObsWindow, names_from = Label,values_from = n)

### V. Supplementary Table 2: Discrimination of threshold-level GCSm detection models per observation window

GCSm.AUC.table <- read.csv('../results/GCSm_threshold_prediction/compiled_metrics.csv') %>%
  filter(Metrics == 'AUC') %>%
  mutate(FormattedAUC = sprintf('%0.2f (%0.2f–%0.2f)',meanValue,lowerValues,upperValues)) %>%
  pivot_wider(id_cols = ObsWindow, names_from = Threshold, values_from = FormattedAUC) %>%
  arrange(ObsWindow)

### VI. Supplementary Table 3: Count distributions of GOSE scores at hospital discharge per observation window

# Load full matrix keys (one imputation is sufficient) to get GOSE scores per observation
full.matrix.keys <- list.files('../features/03_formatted_predictor_matrices/full_matrices',
                               pattern = glob2rx('*_h_imputation_1_keys.csv'),
                               full.names = T)

# Initialize empty dataframe to store count distributions for table
patient.outcomes <- read.csv('../clinical_data/patient_outcomes.csv')

# Initialize empty dataframe to store count distributions for table
compiled.count.dist.table <- data.frame(matrix(ncol = 3,nrow = 0))

# Iterate through full matrix keys
for (curr.full.matrix.key.file in full.matrix.keys){
  
  # Load current full matrix key
  curr.key <- read.csv(curr.full.matrix.key.file)
  
  # Extract current observation window from key file name
  curr.obs.window <- as.numeric(str_match(curr.full.matrix.key.file, "full_matrices/\\s*(.*?)\\s*_h_imputation_")[,2])
  
  # Calculate current GOSE (discharge) frequencies from the full matrix key
  GOSE.freqs <- curr.key %>%
    left_join(patient.outcomes %>% select(UPI,GOSEDischarge), by = 'UPI') %>%
    group_by(GOSEDischarge) %>%
    tally() %>%
    drop_na() %>%
    mutate(GOSEDischarge = paste0('GOSE=',GOSEDischarge)) %>%
    rename(Label = GOSEDischarge) %>%
    mutate(n = sprintf('%01.f (%0.2f&)',n,(100*n/sum(.$n))))
  
  # Add information about unique patient count and total observation count and append to compiled distribution dataframe
  compiled.count.dist.table <- rbind(compiled.count.dist.table,
                                     rbind(data.frame(Label = 'n_unique_patients',n = length(unique(curr.key$UPI))),
                                           data.frame(Label = 'n_total_observations',n = nrow(curr.key)),
                                           GOSE.freqs) %>%
                                       mutate(ObsWindow = curr.obs.window))
}

# Pivot compiled count distribution wider for manuscript supplementary table
compiled.count.dist.table <- compiled.count.dist.table %>%
  pivot_wider(id_cols = ObsWindow, names_from = Label,values_from = n)

### VII. Supplementary Table 4: Discrimination of threshold-level GOSE at hospital discharge prediction models per observation window

GOSE.AUC.table <- read.csv('../results/GOSE_threshold_prediction/compiled_metrics.csv') %>%
  filter(Metrics == 'AUC') %>%
  mutate(FormattedAUC = sprintf('%0.2f (%0.2f–%0.2f)',meanValue,lowerValues,upperValues)) %>%
  pivot_wider(id_cols = ObsWindow, names_from = Threshold, values_from = FormattedAUC) %>%
  arrange(ObsWindow)

### VIII. Supplementary Table 5: Count distributions of GOSE scores at 12 months post discharge per observation window

# Load full matrix keys (one imputation is sufficient) to get GOSE (12m) scores per observation
full.matrix.keys <- list.files('../features/03_formatted_predictor_matrices/full_matrices',
                               pattern = glob2rx('*_h_imputation_1_keys.csv'),
                               full.names = T)

# Initialize empty dataframe to store count distributions for table
patient.outcomes <- read.csv('../clinical_data/patient_outcomes.csv')

# Initialize empty dataframe to store count distributions for table
compiled.count.dist.table <- data.frame(matrix(ncol = 3,nrow = 0))

# Iterate through full matrix keys
for (curr.full.matrix.key.file in full.matrix.keys){
  
  # Load current full matrix key
  curr.key <- read.csv(curr.full.matrix.key.file)
  
  # Extract current observation window from key file name
  curr.obs.window <- as.numeric(str_match(curr.full.matrix.key.file, "full_matrices/\\s*(.*?)\\s*_h_imputation_")[,2])
  
  # Calculate current GOSE (12 months) frequencies from the full matrix key
  GOSE12m.freqs <- curr.key %>%
    left_join(patient.outcomes %>% select(UPI,GOSE12Months), by = 'UPI') %>%
    group_by(GOSE12Months) %>%
    tally() %>%
    drop_na() %>%
    mutate(GOSE12Months = paste0('GOSE12m=',GOSE12Months)) %>%
    rename(Label = GOSE12Months) %>%
    mutate(n = sprintf('%01.f (%0.2f&)',n,(100*n/sum(.$n))))
  
  # Add information about unique patient count and total observation count and append to compiled distribution dataframe
  unique.key.UPIs <- nrow(curr.key %>% left_join(patient.outcomes %>% 
                                                   select(UPI,GOSE12Months), by = 'UPI') %>% 
                            drop_na(GOSE12Months) %>%
                            select(UPI,GOSE12Months) %>%
                            distinct())
  
  compiled.count.dist.table <- rbind(compiled.count.dist.table,
                                     rbind(data.frame(Label = 'n_unique_patients',n = (unique.key.UPIs)),
                                           data.frame(Label = 'n_total_observations',n = nrow(curr.key)),
                                           GOSE12m.freqs) %>%
                                       mutate(ObsWindow = curr.obs.window))
}

# Pivot compiled count distribution wider for manuscript supplementary table
compiled.count.dist.table <- compiled.count.dist.table %>%
  pivot_wider(id_cols = ObsWindow, names_from = Label,values_from = n)

### IX. Supplementary Table 6: Discrimination of threshold-level GOSE at 12 months post discharge prediction models per observation window

GOSE12m.AUC.table <- read.csv('../results/GOSE12m_threshold_prediction/compiled_metrics.csv') %>%
  filter(Metrics == 'AUC') %>%
  mutate(FormattedAUC = sprintf('%0.2f (%0.2f–%0.2f)',meanValue,lowerValues,upperValues)) %>%
  pivot_wider(id_cols = ObsWindow, names_from = Threshold, values_from = FormattedAUC) %>%
  arrange(ObsWindow)

### X. Supplementary Table 7: Percentages of missing accelerometry data per sensor and recording duration of each study participant
## Calculate patient-specific missing and static activity information
patient.specific.recording.missing.static.info <- read_csv('../features/all_features.csv') %>%
  filter(Feature == 'SMA') %>%
  dplyr::select(-Feature) %>%
  pivot_longer(cols = -c(UPI,RecordingIdx,HoursFromICUAdmission,TimeOfDay),names_to = 'Sensor') %>%
  group_by(UPI,Sensor) %>%
  summarise(total.duration.hours = n()/720,
            total.missing.perc = 100*sum(is.na(value))/n(),
            total.static.perc = 100*sum(value < .135,na.rm = T)/n())

missingness.table <- patient.specific.recording.missing.static.info %>%
  dplyr::select(UPI,Sensor,total.duration.hours,total.missing.perc) %>%
  pivot_wider(id_cols = c(UPI,total.duration.hours),names_from = 'Sensor',values_from = total.missing.perc)

### XI. Metrics for Figure 5: Feature significance matrices of optimally discriminating motor function detection and functional outcome prediction models

## Initialize parallel bootstrapping parameters
# Number of boostrap resamples
NUM.BOOTSTRAPS <- 1000

# Number of cores to use in parallel
NUM.CORES <- 10

# Initialize `doParallel` cluster
registerDoParallel(cores = NUM.CORES)

## Feature significance values of threshold-level GCSm detection (GCSm > 4, Obs. Window = 6 hr)
# Load feature significance values of current optimal-AUC configuration
GCSm.feature.sig.values <-  read.csv('../features/03_formatted_predictor_matrices/feature_analysis/06.00_h_obs_window/GCSm.gt.4_feature_analysis_values.csv')

# Identify unique sensor-feature combinations
unique.sensor.feature.combos <- GCSm.feature.sig.values %>% 
  dplyr::select(Sensor, Feature) %>%
  distinct()

# Initialize dataframe to store bootstrapped results across combinations
bs.GCSm.feature.sig.values <- data.frame(matrix(ncol = 5,nrow = 0))

# Iterate through unique sensor-feature combinations
for (curr.combo.idx in 1:nrow(unique.sensor.feature.combos)){
  # Filter out feature significance values of current combination
  curr.sensor <- unique.sensor.feature.combos$Sensor[curr.combo.idx]
  curr.feature <- unique.sensor.feature.combos$Feature[curr.combo.idx]
  filt.GCSm.feature.sig.values <- GCSm.feature.sig.values %>%
    filter(Sensor == curr.sensor, Feature == curr.feature)
  
  # In parallel, bootstrap N feature significance means, medians, and max values
  curr.bs.values <- foreach(icount(NUM.BOOTSTRAPS), .combine=rbind) %dopar% {
    ind <- sample(nrow(filt.GCSm.feature.sig.values), nrow(filt.GCSm.feature.sig.values), replace=TRUE)
    meanSignificance <- mean(filt.GCSm.feature.sig.values[ind,'Significance'],na.rm = T)
    maxSignificance <- max(filt.GCSm.feature.sig.values[ind,'Significance'],na.rm = T)
    medianSignificance <- median(filt.GCSm.feature.sig.values[ind,'Significance'],na.rm = T)
    data.frame(meanSignificance,maxSignificance,medianSignificance)
  }
  
  # Add information about current sensor-feature combination
  curr.bs.values <- curr.bs.values %>%
    mutate(Sensor = curr.sensor,
           Feature = curr.feature) %>%
    relocate(Sensor, Feature)
  
  # Append to running dataframe
  bs.GCSm.feature.sig.values <- rbind(bs.GCSm.feature.sig.values,curr.bs.values)
  
  # Status update
  print(paste0('Combination no. ',curr.combo.idx,' out of ',nrow(unique.sensor.feature.combos),' completed.'))
}

# Save bootstrapped feature significance results
write.csv(bs.GCSm.feature.sig.values,'../results/GCSm_threshold_prediction/feature_significance.csv',row.names = F)

## Feature significance values of threshold-level GOSE (discharge) prediction (GOSE > 5, Obs. Window = 6 hr)
# Load feature significance values of current optimal-AUC configuration
GOSE.feature.sig.values <-  read.csv('../features/03_formatted_predictor_matrices/feature_analysis/06.00_h_obs_window/GOSE.gt.5_feature_analysis_values.csv')

# Identify unique sensor-feature combinations
unique.sensor.feature.combos <- GOSE.feature.sig.values %>% 
  dplyr::select(Sensor, Feature) %>%
  distinct()

# Initialize dataframe to store bootstrapped results across combinations
bs.GOSE.feature.sig.values <- data.frame(matrix(ncol = 5,nrow = 0))

# Iterate through unique sensor-feature combinations
for (curr.combo.idx in 1:nrow(unique.sensor.feature.combos)){
  # Filter out feature significance values of current combination
  curr.sensor <- unique.sensor.feature.combos$Sensor[curr.combo.idx]
  curr.feature <- unique.sensor.feature.combos$Feature[curr.combo.idx]
  filt.GOSE.feature.sig.values <- GOSE.feature.sig.values %>%
    filter(Sensor == curr.sensor, Feature == curr.feature)
  
  # In parallel, bootstrap N feature significance means, medians, and max values
  curr.bs.values <- foreach(icount(NUM.BOOTSTRAPS), .combine=rbind) %dopar% {
    ind <- sample(nrow(filt.GOSE.feature.sig.values), nrow(filt.GOSE.feature.sig.values), replace=TRUE)
    meanSignificance <- mean(filt.GOSE.feature.sig.values[ind,'Significance'],na.rm = T)
    maxSignificance <- max(filt.GOSE.feature.sig.values[ind,'Significance'],na.rm = T)
    medianSignificance <- median(filt.GOSE.feature.sig.values[ind,'Significance'],na.rm = T)
    data.frame(meanSignificance,maxSignificance,medianSignificance)
  }
  
  # Add information about current sensor-feature combination
  curr.bs.values <- curr.bs.values %>%
    mutate(Sensor = curr.sensor,
           Feature = curr.feature) %>%
    relocate(Sensor, Feature)
  
  # Append to running dataframe
  bs.GOSE.feature.sig.values <- rbind(bs.GOSE.feature.sig.values,curr.bs.values)
  
  # Status update
  print(paste0('Combination no. ',curr.combo.idx,' out of ',nrow(unique.sensor.feature.combos),' completed.'))
}

# Save bootstrapped feature significance results
write.csv(bs.GOSE.feature.sig.values,'../results/GOSE_threshold_prediction/feature_significance.csv',row.names = F)

## Stop implicit cluster used in parallel processing
stopImplicitCluster()

### XII. Metrics for Supplementary Figure 2: Correlation matrices of extracted motion features across different sensor placements

## Initialize parallel bootstrapping parameters
# Number of boostrap resamples
NUM.BOOTSTRAPS <- 1000

# Number of cores to use in parallel
NUM.CORES <- 10

# Initialize `doParallel` cluster
registerDoParallel(cores = NUM.CORES)

## Calculate correlation of features across pairwise sensor combinations per patient
# Load all motion features
if (!exists("all.motion.features")) {
  all.motion.features <- read_csv('../features/all_features.csv')
}

# Calcualte pairwise (sensor) Spearman correlation coefficients for each feature for each patient
sensor.correlation.df <- all.motion.features %>%
  dplyr::select(-c(RecordingIdx,HoursFromICUAdmission,TimeOfDay)) %>%
  group_by(UPI,Feature) %>%
  nest() %>%
  mutate(corr.df = map(data,correlate,method = "spearman",diagonal = 1,quiet = T)) %>% 
  unnest(corr.df) %>%
  dplyr::select(-data) %>%
  rename(term1 = term) %>%
  pivot_longer(cols = -c(UPI,Feature,term1),names_to = 'term2',values_to = 'rho')

## Bootstrap to calculate 95% confidence intervals for correlation coefficients
# In parallel, bootstrap N correlation coefficient means
curr.bs.coeffs <- foreach(icount(NUM.BOOTSTRAPS), .combine=rbind) %dopar% {
  curr.UPIs <- unique(sample(unique(sensor.correlation.df$UPI),length(unique(sensor.correlation.df$UPI)),replace = T))
  
  sensor.correlation.df %>%
    filter(UPI %in% curr.UPIs) %>%
    group_by(Feature,term1,term2) %>%
    summarise(rho = mean(rho,na.rm=T))
}

# Calculate confidence intervals for each coefficient
sensor.correlations.CI <- curr.bs.coeffs %>%
  group_by(Feature,term1,term2) %>%
  summarise(meanRho = mean(rho,na.rm = T),
            lowerRho = quantile(rho,.025,na.rm=T),
            upperRho = quantile(rho,.975,na.rm=T)) %>%
  mutate(FormattedRho = sprintf('%0.2f \n (%0.2f–%0.2f)',meanRho,lowerRho,upperRho))

# Fix diagonal correlation coefficients to 1
sensor.correlations.CI$FormattedRho[sensor.correlations.CI$term1 == sensor.correlations.CI$term2] <- '1.00'

# Reformat feature names for subsequent plotting
sensor.correlations.CI$Feature <- plyr::mapvalues(
  sensor.correlations.CI$Feature,
  from = c("SMA","HLF_h","HLF_l","MFR","FDE","BPW","WVL"),
  to = c("SMA","HLF (h)","HLF (l)","MFR","FDE","BPW","WVL")
)

# Save sensor correlation matrix and confidence interval for subsequent plotting
write.csv(sensor.correlations.CI,'../summary_statistics/sensor_correlations.csv',row.names = F)

### XIII. Metrics for Supplementary Figure 4: Mean motion feature trajectories in the six hours preceding GCSm evaluation, stratified by GCSm scores and bilateral sensor placement

## Initialize parallel bootstrapping parameters
# Number of boostrap resamples
NUM.BOOTSTRAPS <- 1000

# Number of cores to use in parallel
NUM.CORES <- 10

# Initialize `doParallel` cluster
registerDoParallel(cores = NUM.CORES)

## Prepare features and feature metadata for mean value calculation
# Load GCSm and UPI keys corresponding to 6-hour observation window motion features
keys.6hr <- read.csv('../features/03_formatted_predictor_matrices/full_matrices/06.00_h_imputation_1_keys.csv')

# Load 6-hour observation window motion features (only first imputation)
motion.features.6hr <- readRDS('../features/03_formatted_predictor_matrices/full_matrices/06.00_h_imputation_1_full_matrix.rds') %>%
  as.data.frame() %>%
  mutate(UPI = keys.6hr$UPI,
         GCSm = factor(keys.6hr$GCSm),
         ObsIdx = 1:nrow(.)) %>%
  pivot_longer(cols = -c(UPI,GCSm,ObsIdx)) %>%
  mutate(Sensor = sub("\\/.*", "", name),
         WindowsBeforeEvaluation = as.numeric(sub(".*/", "", name)),
         Feature = str_match(name, "/\\s*(.*?)\\s*/")[,2]) %>%
  dplyr::select(-name) %>%
  relocate(ObsIdx,UPI,GCSm,Sensor,Feature,WindowsBeforeEvaluation,value) %>%
  filter(Feature != 'PhysActivity') %>%
  mutate(Feature = factor(Feature,
                          levels = c("SMA","HLF_h","HLF_l","MFR","FDE","BPW","WVL")))

# Create new variable for bilateral sensor placement
motion.features.6hr$Placement <-
  factor(plyr::mapvalues(
    motion.features.6hr$Sensor,
    from = c("RE", "LE", "RW", "LW", "RA", "LA"),
    to = c("Elbows", "Elbows", "Wrists", "Wrists", "Ankles", "Ankles")
  ),levels = c("Elbows","Wrists","Ankles"))

# Discretize `WindowsBeforeEvaluation` into uniform 10 minute bins before GCSm evaluation
motion.features.6hr$WindowsBeforeEvaluation <- discretize(motion.features.6hr$WindowsBeforeEvaluation,36) %>% predict(motion.features.6hr$WindowsBeforeEvaluation)
motion.features.6hr$WindowsBeforeEvaluation <- as.numeric(sub(".*bin", "", motion.features.6hr$WindowsBeforeEvaluation))

# Mean feature values in each 10-minute bin for each observation, GCSm, feature, and placement combination
motion.features.6hr <- motion.features.6hr %>%
  group_by(ObsIdx,GCSm,Feature,Placement,WindowsBeforeEvaluation) %>%
  summarise(value = mean(value,na.rm=T))

# Define new function to remove outliers beyond 3*IQR of 3rd quartile
remove.outliers <- function(df) {
  x <- df$value
  y <- x
  H <- 3 * IQR(x, na.rm = T)
  qnt <- quantile(x, probs=c(.25, .75), na.rm = T)
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  df$value <- y
  return(df)
}

# Clean out outliers of features for each combination of feature and bilateral placement
cleaned.motion.features.6hr <- motion.features.6hr %>%
  group_by(Feature,Placement) %>%
  nest() %>%
  mutate(data = map(data, remove.outliers)) %>%
  unnest(cols = c(data)) %>% 
  drop_na(value)

# Identify unique observation indices
uniq.ObsIdx <- unique(cleaned.motion.features.6hr$ObsIdx)

## Calculate mean trajectories and accompanying 95% confidence intervals
# Parallelized foreach loop:
curr.bs.means <- foreach(icount(NUM.BOOTSTRAPS), .combine=rbind, .inorder=FALSE) %dopar% {
  
  curr.ObsIdx.resample <- unique(sample(uniq.ObsIdx,length(uniq.ObsIdx),replace = T))
  
  cleaned.motion.features.6hr %>%
    filter(ObsIdx %in% curr.ObsIdx.resample) %>%
    group_by(GCSm,Feature,Placement,WindowsBeforeEvaluation) %>%
    summarise(meanValue = mean(value,na.rm=T))
}

# Calculate 95% confidence intervals from bootstrapped mean values
trajectory.means.CI <- curr.bs.means %>%
  group_by(GCSm,Feature,Placement,WindowsBeforeEvaluation) %>%
  summarise(meanValues = mean(meanValue,na.rm = T),
            lowerValues = quantile(meanValue,.025,na.rm=T),
            upperValues = quantile(meanValue,.975,na.rm=T)) %>%
  mutate(HoursBeforeEvaluation = (1/6)*(WindowsBeforeEvaluation - .5)) %>%
  arrange(Feature,Placement,GCSm,desc(HoursBeforeEvaluation)) %>%
  relocate(Feature,Placement,GCSm,WindowsBeforeEvaluation,HoursBeforeEvaluation)

# Save feature mean trajectories and associated confidence intervals for subsequent plotting
write.csv(trajectory.means.CI,'../summary_statistics/feature_mean_trajectories.csv',row.names = F)

### XIV. Miscellaneous statistics for manuscript