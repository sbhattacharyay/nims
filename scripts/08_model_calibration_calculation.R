#### Master Script 8: Calculate model calibration on validation set predictions ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# University of Cambridge
# Johns Hopkins University
# email address: sb2406@cam.ac.uk
#
### Contents:
# I. Initialization
# II. Calculate calibration curves of optimal (maximum AUC) GCSm-threshold detection model types
# III. Calculate calibration curves of optimal (maximum AUC) GOSE-threshold detection model types 
# IV. Calculate calibration curves of optimal (maximum AUC) GOSE12m-threshold detection model types 

## NOTE: Requires compilation of AUC metrics (performed in Master script 10)

### I. Initialization
# Load necessary packages
library(tidyverse)
library(rms)
library(foreach)
library(doParallel)

# Set the number of parallel cores
no.parallel.cores <- 10
registerDoParallel(cores = no.parallel.cores)

# Define number of bootstrap resamples
NUM.BOOTSTRAPS <- 1000

# Define common x-axis for calibration curve interpolation
X.INTERP <- seq(0,1,length = 200)

# Define function to calculate Integrated Calibration Index (ICI)
ICI <- function(probs,labels){
  calibrated.model <- lowess(probs,labels,iter=0)
  return(mean(abs(calibrated.model$y - calibrated.model$x)))
}

# Define function to calculate E50
E50 <- function(probs,labels){
  calibrated.model <- lowess(probs,labels,iter=0)
  return(median(abs(calibrated.model$y - calibrated.model$x)))
}

# Define function to calculate E90
E90 <- function(probs,labels){
  calibrated.model <- lowess(probs,labels,iter=0)
  return(quantile(abs(calibrated.model$y - calibrated.model$x),0.9))
}

# Define function to calculate Emax
Emax <- function(probs,labels){
  calibrated.model <- lowess(probs,labels,iter=0)
  return(max(abs(calibrated.model$y - calibrated.model$x)))
}

### II. Calculate calibration curves of optimal (maximum AUC) GCSm-threshold detection model types 
## Determine optimal observation window for each threshold
# Load compiled AUC metrics
compiled.GCSm.AUC.df <- read.csv('../results/GCSm_threshold_prediction/compiled_metrics.csv') %>%
  filter(Metrics == 'AUC')

# First, filter out significantly discriminating AUC, and determine maximum per observation window
opt.GCSm.signficant.AUC.df <- compiled.GCSm.AUC.df %>%
  filter(lowerValues >= 0.50) %>%
  group_by(Threshold) %>%
  summarise(ObsWindow = ObsWindow[which.max(meanValue)],
            meanAUC = max(meanValue),
            medianAUC = medianValues[which.max(meanValue)],
            lowerAUC = lowerValues[which.max(meanValue)],
            upperAUC = upperValues[which.max(meanValue)])

# Then, for observation windows for which we cannot achieve significant discrimination, determine maximum AUC
opt.GCSm.nonsignficant.AUC.df <- compiled.GCSm.AUC.df %>%
  filter(!(Threshold %in% opt.GCSm.signficant.AUC.df$Threshold)) %>%
  group_by(Threshold) %>%
  summarise(ObsWindow = ObsWindow[which.max(meanValue)],
            meanAUC = max(meanValue),
            medianAUC = medianValues[which.max(meanValue)],
            lowerAUC = lowerValues[which.max(meanValue)],
            upperAUC = upperValues[which.max(meanValue)])

# Compile the significant and non-significant AUC values to determine ROCs for plotting
opt.GCSm.AUC.df <- rbind(opt.GCSm.signficant.AUC.df,opt.GCSm.nonsignficant.AUC.df)

## Calculate calibration curve (with 95% CI - BBC-CV) of the optimal configurations
# Initialize empty dataframe to store ICIs across thresholds
compiled.GCSm.threshold.metrics <- as.data.frame(matrix(ncol = 8, nrow = 0))

# Iterate through each GCSm threshold
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
  compiled.GCSm.threshold.calibration <- foreach(icount(NUM.BOOTSTRAPS), .combine=rbind) %dopar%{
    # Keep drawing sample until both cases are present in both in- and out-sample cases
    fail.condition <- TRUE
    while(fail.condition){
      curr.UPI.resample <- sample(unique.UPIs,length(unique.UPIs),replace = T)
      # Extract in-sample and out-sample UPIs
      curr.in.sample <- sort(unique(curr.UPI.resample))
      curr.out.sample <- sort(unique.UPIs[! unique.UPIs %in% curr.in.sample])
      # Divide in-sample and out-sample predictions
      curr.in.sample.preds <- curr.predictions %>% filter(UPI %in% curr.in.sample)
      curr.out.sample.preds <- curr.predictions %>% filter(UPI %in% curr.out.sample)
      # If the necessary condition is met, we may break out of the while loop
      if ((length(unique(curr.in.sample.preds$TrueLabel)) == 2) & 
          (length(unique(curr.out.sample.preds$TrueLabel)) == 2)){
        fail.condition <- FALSE
      }
    }
    
    # Determine optimal configuration in current resample for calibration based on in-sample ICI
    opt.config <- curr.in.sample.preds %>%
      group_by(ConfigIdx) %>%
      summarise(ICIvalues = ICI(Prob,TrueLabel)) %>%
      top_n(-1,ICIvalues)
    
    # Calculate ICI and calibration curve for current optimal configuration
    curr.config.out.sample.preds <- curr.out.sample.preds %>% filter(ConfigIdx == opt.config$ConfigIdx)
    curr.ICI <- ICI(curr.config.out.sample.preds$Prob,curr.config.out.sample.preds$TrueLabel)
    curr.E50 <- E50(curr.config.out.sample.preds$Prob,curr.config.out.sample.preds$TrueLabel)
    curr.E90 <- E90(curr.config.out.sample.preds$Prob,curr.config.out.sample.preds$TrueLabel)
    curr.Emax <- Emax(curr.config.out.sample.preds$Prob,curr.config.out.sample.preds$TrueLabel)
    curr.config.out.calibrated.model <- lowess(curr.config.out.sample.preds$Prob,curr.config.out.sample.preds$TrueLabel,iter=0)
    curr.calib.curve <- approx(x = curr.config.out.calibrated.model$x,
                               y = curr.config.out.calibrated.model$y,
                               xout = X.INTERP)
    
    # Return dataframe row of compiled information
    data.frame(
      Threshold = curr.thresh,
      ObsWindow = curr.obs.window,
      ConfigIdx = opt.config$ConfigIdx,
      PredProb = curr.calib.curve$x,
      TrueProb = curr.calib.curve$y,
      ICI = curr.ICI,
      E50 = curr.E50,
      E90 = curr.E90,
      Emax = curr.Emax
    )
  }
  
  # Derive 95% confidence intervals from compiled calibration axes
  summ.GCSm.threshold.calibration <- compiled.GCSm.threshold.calibration %>%
    group_by(Threshold,ObsWindow,PredProb) %>%
    summarise(meanTrueProb = mean(TrueProb,na.rm = T),
              medianTrueProb = median(TrueProb,na.rm = T),
              lowerTrueProb = quantile(TrueProb,.025,na.rm = T,),
              upperTrueProb = quantile(TrueProb,.975,na.rm = T,))
  
  # Save summarized calibration curve
  write.csv(summ.GCSm.threshold.calibration,
            file.path('../results/GCSm_threshold_prediction',
                        paste0(curr.thresh, '_calibration_curve.csv')),
            row.names = F)
  
  # Derive 95% of compiled ICI values as well as optimal configuration index
  summ.GCSm.threshold.metrics <- compiled.GCSm.threshold.calibration %>%
    dplyr::select(-c(PredProb,TrueProb)) %>%
    distinct() %>%
    pivot_longer(cols = -c(Threshold,ObsWindow,ConfigIdx),names_to = 'Metric',values_to = 'Value') %>%
    group_by(Threshold,ObsWindow,Metric) %>%
    summarise(meanValue = mean(Value,na.rm = T),
              medianValue = median(Value,na.rm = T),
              lowerValue = quantile(Value,.025,na.rm = T,),
              upperValue = quantile(Value,.975,na.rm = T,),
              optConfigIdx = unique(ConfigIdx)[which.max(table(ConfigIdx))])
  
  # Append to compiled metric dataframe
  compiled.GCSm.threshold.metrics <- rbind(compiled.GCSm.threshold.metrics,summ.GCSm.threshold.metrics)
  
  # Status update on completion of threshold
  print(paste(curr.thresh,'complete.'))
}

# Save compiled metrics
write.csv(compiled.GCSm.threshold.metrics,'../results/GCSm_threshold_prediction/calibration_metrics.csv',row.names = F)

### III. Calculate calibration curves of optimal (maximum AUC) GOSE-threshold detection model types 
## Determine optimal observation window for each threshold
# Load compiled AUC metrics
compiled.GOSE.AUC.df <- read.csv('../results/GOSE_threshold_prediction/compiled_metrics.csv') %>%
  filter(Metrics == 'AUC')

# First, filter out significantly discriminating AUC, and determine maximum per observation window
opt.GOSE.signficant.AUC.df <- compiled.GOSE.AUC.df %>%
  filter(lowerValues >= 0.50) %>%
  group_by(Threshold) %>%
  summarise(ObsWindow = ObsWindow[which.max(meanValue)],
            meanAUC = max(meanValue),
            medianAUC = medianValues[which.max(meanValue)],
            lowerAUC = lowerValues[which.max(meanValue)],
            upperAUC = upperValues[which.max(meanValue)])

# Then, for observation windows for which we cannot achieve significant discrimination, determine maximum AUC
opt.GOSE.nonsignficant.AUC.df <- compiled.GOSE.AUC.df %>%
  filter(!(Threshold %in% opt.GOSE.signficant.AUC.df$Threshold)) %>%
  group_by(Threshold) %>%
  summarise(ObsWindow = ObsWindow[which.max(meanValue)],
            meanAUC = max(meanValue),
            medianAUC = medianValues[which.max(meanValue)],
            lowerAUC = lowerValues[which.max(meanValue)],
            upperAUC = upperValues[which.max(meanValue)])

# Compile the significant and non-significant AUC values to determine ROCs for plotting
opt.GOSE.AUC.df <- rbind(opt.GOSE.signficant.AUC.df,opt.GOSE.nonsignficant.AUC.df)

## Calculate calibration curve (with 95% CI - BBC-CV) of the optimal configurations
# Initialize empty dataframe to store ICIs across thresholds
compiled.GOSE.threshold.metrics <- as.data.frame(matrix(ncol = 8, nrow = 0))

# Iterate through each GOSE threshold
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
  compiled.GOSE.threshold.calibration <- foreach(icount(NUM.BOOTSTRAPS), .combine=rbind) %dopar%{
    # Keep drawing sample until both cases are present in both in- and out-sample cases
    fail.condition <- TRUE
    while(fail.condition){
      curr.UPI.resample <- sample(unique.UPIs,length(unique.UPIs),replace = T)
      # Extract in-sample and out-sample UPIs
      curr.in.sample <- sort(unique(curr.UPI.resample))
      curr.out.sample <- sort(unique.UPIs[! unique.UPIs %in% curr.in.sample])
      # Divide in-sample and out-sample predictions
      curr.in.sample.preds <- curr.predictions %>% filter(UPI %in% curr.in.sample)
      curr.out.sample.preds <- curr.predictions %>% filter(UPI %in% curr.out.sample)
      # If the necessary condition is met, we may break out of the while loop
      if ((length(unique(curr.in.sample.preds$TrueLabel)) == 2) & 
          (length(unique(curr.out.sample.preds$TrueLabel)) == 2)){
        fail.condition <- FALSE
      }
    }
    
    # Determine optimal configuration in current resample for calibration based on in-sample ICI
    opt.config <- curr.in.sample.preds %>%
      group_by(ConfigIdx) %>%
      summarise(ICIvalues = ICI(Prob,TrueLabel)) %>%
      top_n(-1,ICIvalues)
    
    # Calculate ICI and calibration curve for current optimal configuration
    curr.config.out.sample.preds <- curr.out.sample.preds %>% filter(ConfigIdx == opt.config$ConfigIdx)
    curr.ICI <- ICI(curr.config.out.sample.preds$Prob,curr.config.out.sample.preds$TrueLabel)
    curr.E50 <- E50(curr.config.out.sample.preds$Prob,curr.config.out.sample.preds$TrueLabel)
    curr.E90 <- E90(curr.config.out.sample.preds$Prob,curr.config.out.sample.preds$TrueLabel)
    curr.Emax <- Emax(curr.config.out.sample.preds$Prob,curr.config.out.sample.preds$TrueLabel)
    curr.config.out.calibrated.model <- lowess(curr.config.out.sample.preds$Prob,curr.config.out.sample.preds$TrueLabel,iter=0)
    curr.calib.curve <- approx(x = curr.config.out.calibrated.model$x,
                               y = curr.config.out.calibrated.model$y,
                               xout = X.INTERP)
    
    # Return dataframe row of compiled information
    data.frame(
      Threshold = curr.thresh,
      ObsWindow = curr.obs.window,
      ConfigIdx = opt.config$ConfigIdx,
      PredProb = curr.calib.curve$x,
      TrueProb = curr.calib.curve$y,
      ICI = curr.ICI,
      E50 = curr.E50,
      E90 = curr.E90,
      Emax = curr.Emax
    )
  }
  
  # Derive 95% confidence intervals from compiled calibration axes
  summ.GOSE.threshold.calibration <- compiled.GOSE.threshold.calibration %>%
    group_by(Threshold,ObsWindow,PredProb) %>%
    summarise(meanTrueProb = mean(TrueProb,na.rm = T),
              medianTrueProb = median(TrueProb,na.rm = T),
              lowerTrueProb = quantile(TrueProb,.025,na.rm = T,),
              upperTrueProb = quantile(TrueProb,.975,na.rm = T,))
  
  # Save summarized calibration curve
  write.csv(summ.GOSE.threshold.calibration,
            file.path('../results/GOSE_threshold_prediction',
                      paste0(curr.thresh, '_calibration_curve.csv')),
            row.names = F)
  
  # Derive 95% of compiled ICI values as well as optimal configuration index
  summ.GOSE.threshold.metrics <- compiled.GOSE.threshold.calibration %>%
    dplyr::select(-c(PredProb,TrueProb)) %>%
    distinct() %>%
    pivot_longer(cols = -c(Threshold,ObsWindow,ConfigIdx),names_to = 'Metric',values_to = 'Value') %>%
    group_by(Threshold,ObsWindow,Metric) %>%
    summarise(meanValue = mean(Value,na.rm = T),
              medianValue = median(Value,na.rm = T),
              lowerValue = quantile(Value,.025,na.rm = T,),
              upperValue = quantile(Value,.975,na.rm = T,),
              optConfigIdx = unique(ConfigIdx)[which.max(table(ConfigIdx))])
  
  # Append to compiled metric dataframe
  compiled.GOSE.threshold.metrics <- rbind(compiled.GOSE.threshold.metrics,summ.GOSE.threshold.metrics)
  
  # Status update on completion of threshold
  print(paste(curr.thresh,'complete.'))
}

# Save compiled metrics
write.csv(compiled.GOSE.threshold.metrics,'../results/GOSE_threshold_prediction/calibration_metrics.csv',row.names = F)

### IV. Calculate calibration curves of optimal (maximum AUC) GOSE12m-threshold detection model types 
## Determine optimal observation window for each threshold
# Load compiled AUC metrics
compiled.GOSE12m.AUC.df <- read.csv('../results/GOSE12m_threshold_prediction/compiled_metrics.csv') %>%
  filter(Metrics == 'AUC')

# First, filter out significantly discriminating AUC, and determine maximum per observation window
opt.GOSE12m.signficant.AUC.df <- compiled.GOSE12m.AUC.df %>%
  filter(lowerValues >= 0.50) %>%
  group_by(Threshold) %>%
  summarise(ObsWindow = ObsWindow[which.max(meanValue)],
            meanAUC = max(meanValue),
            medianAUC = medianValues[which.max(meanValue)],
            lowerAUC = lowerValues[which.max(meanValue)],
            upperAUC = upperValues[which.max(meanValue)])

# Then, for observation windows for which we cannot achieve significant discrimination, determine maximum AUC
opt.GOSE12m.nonsignficant.AUC.df <- compiled.GOSE12m.AUC.df %>%
  filter(!(Threshold %in% opt.GOSE12m.signficant.AUC.df$Threshold)) %>%
  group_by(Threshold) %>%
  summarise(ObsWindow = ObsWindow[which.max(meanValue)],
            meanAUC = max(meanValue),
            medianAUC = medianValues[which.max(meanValue)],
            lowerAUC = lowerValues[which.max(meanValue)],
            upperAUC = upperValues[which.max(meanValue)])

# Compile the significant and non-significant AUC values to determine ROCs for plotting
opt.GOSE12m.AUC.df <- rbind(opt.GOSE12m.signficant.AUC.df,opt.GOSE12m.nonsignficant.AUC.df)

## Calculate calibration curve (with 95% CI - BBC-CV) of the optimal configurations
# Initialize empty dataframe to store ICIs across thresholds
compiled.GOSE12m.threshold.metrics <- as.data.frame(matrix(ncol = 8, nrow = 0))

# Iterate through each GOSE12m threshold
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
  compiled.GOSE12m.threshold.calibration <- foreach(icount(NUM.BOOTSTRAPS), .combine=rbind) %dopar%{
    # Keep drawing sample until both cases are present in both in- and out-sample cases
    fail.condition <- TRUE
    while(fail.condition){
      curr.UPI.resample <- sample(unique.UPIs,length(unique.UPIs),replace = T)
      # Extract in-sample and out-sample UPIs
      curr.in.sample <- sort(unique(curr.UPI.resample))
      curr.out.sample <- sort(unique.UPIs[! unique.UPIs %in% curr.in.sample])
      # Divide in-sample and out-sample predictions
      curr.in.sample.preds <- curr.predictions %>% filter(UPI %in% curr.in.sample)
      curr.out.sample.preds <- curr.predictions %>% filter(UPI %in% curr.out.sample)
      # If the necessary condition is met, we may break out of the while loop
      if ((length(unique(curr.in.sample.preds$TrueLabel)) == 2) & 
          (length(unique(curr.out.sample.preds$TrueLabel)) == 2)){
        fail.condition <- FALSE
      }
    }
    
    # Determine optimal configuration in current resample for calibration based on in-sample ICI
    opt.config <- curr.in.sample.preds %>%
      group_by(ConfigIdx) %>%
      summarise(ICIvalues = ICI(Prob,TrueLabel)) %>%
      top_n(-1,ICIvalues)
    
    # Calculate ICI and calibration curve for current optimal configuration
    curr.config.out.sample.preds <- curr.out.sample.preds %>% filter(ConfigIdx == opt.config$ConfigIdx)
    curr.ICI <- ICI(curr.config.out.sample.preds$Prob,curr.config.out.sample.preds$TrueLabel)
    curr.E50 <- E50(curr.config.out.sample.preds$Prob,curr.config.out.sample.preds$TrueLabel)
    curr.E90 <- E90(curr.config.out.sample.preds$Prob,curr.config.out.sample.preds$TrueLabel)
    curr.Emax <- Emax(curr.config.out.sample.preds$Prob,curr.config.out.sample.preds$TrueLabel)
    curr.config.out.calibrated.model <- lowess(curr.config.out.sample.preds$Prob,curr.config.out.sample.preds$TrueLabel,iter=0)
    curr.calib.curve <- approx(x = curr.config.out.calibrated.model$x,
                               y = curr.config.out.calibrated.model$y,
                               xout = X.INTERP)
    
    # Return dataframe row of compiled information
    data.frame(
      Threshold = curr.thresh,
      ObsWindow = curr.obs.window,
      ConfigIdx = opt.config$ConfigIdx,
      PredProb = curr.calib.curve$x,
      TrueProb = curr.calib.curve$y,
      ICI = curr.ICI,
      E50 = curr.E50,
      E90 = curr.E90,
      Emax = curr.Emax
    )
  }
  
  # Derive 95% confidence intervals from compiled calibration axes
  summ.GOSE12m.threshold.calibration <- compiled.GOSE12m.threshold.calibration %>%
    group_by(Threshold,ObsWindow,PredProb) %>%
    summarise(meanTrueProb = mean(TrueProb,na.rm = T),
              medianTrueProb = median(TrueProb,na.rm = T),
              lowerTrueProb = quantile(TrueProb,.025,na.rm = T,),
              upperTrueProb = quantile(TrueProb,.975,na.rm = T,))
  
  # Save summarized calibration curve
  write.csv(summ.GOSE12m.threshold.calibration,
            file.path('../results/GOSE12m_threshold_prediction',
                      paste0(curr.thresh, '_calibration_curve.csv')),
            row.names = F)
  
  # Derive 95% of compiled ICI values as well as optimal configuration index
  summ.GOSE12m.threshold.metrics <- compiled.GOSE12m.threshold.calibration %>%
    dplyr::select(-c(PredProb,TrueProb)) %>%
    distinct() %>%
    pivot_longer(cols = -c(Threshold,ObsWindow,ConfigIdx),names_to = 'Metric',values_to = 'Value') %>%
    group_by(Threshold,ObsWindow,Metric) %>%
    summarise(meanValue = mean(Value,na.rm = T),
              medianValue = median(Value,na.rm = T),
              lowerValue = quantile(Value,.025,na.rm = T,),
              upperValue = quantile(Value,.975,na.rm = T,),
              optConfigIdx = unique(ConfigIdx)[which.max(table(ConfigIdx))])
  
  # Append to compiled metric dataframe
  compiled.GOSE12m.threshold.metrics <- rbind(compiled.GOSE12m.threshold.metrics,summ.GOSE12m.threshold.metrics)
  
  # Status update on completion of threshold
  print(paste(curr.thresh,'complete.'))
}

# Save compiled metrics
write.csv(compiled.GOSE12m.threshold.metrics,'../results/GOSE12m_threshold_prediction/calibration_metrics.csv',row.names = F)

## Stop Implicit Cluster
stopImplicitCluster()