#### Master Script 11: Plot figures for the manuscript ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# University of Cambridge
# Johns Hopkins University
# email address: sb2406@cam.ac.uk
#
### Contents:
# I. Initialization
# II. Figure 1: Accelerometry processing and feature extraction pipeline and experimental paradigm
# III. Figure 2: Discrimination performance of motor function detection models on validation sets
# IV. Figure 3: Discrimination performance of functional outcome at hospital discharge prediction models on validation sets
# V. Figure 4: Probability calibration of optimally discriminating motor function detection and functional outcome prediction models on validation sets
# VI. Figure 5: Feature significance matrices of optimally discriminating motor function detection and functional outcome prediction models
# VII. Figure 6: Retrospective case study analysis of accelerometry-based detection of motor function in six patients who experienced relevant transition
# VIII. Supplementary Figure 1: Discrimination performance of functional outcome at 12 months post discharge prediction models on validation sets
# IX. Supplementary Figure 2: Correlation matrices of extracted motion features across different sensor placements
# X. Supplementary Figure 3: Violin plots of extracted motion feature values (30 min observation window), stratified by bilateral sensor placement and GCSm scores
# XI. Supplementary Figure 4: Mean motion feature trajectories in the six hours preceding GCSm evaluation, stratified by GCSm scores and bilateral sensor placement
# XII. Supplementary Figure 5: Trajectories of motor component scores of the Glasgow Coma Scale (GCSm) of each study participant during ICU stay
# XIII. Supplementary Figure 7: Percentages of missing, static, and dynamic accelerometry data by time of day of recording and sensor placement
# XIV. Supplementary Figure 8: Count histograms of accelerometry recording information

### I. Initialization
## Import necessary packages
library(tidyverse)
library(ggpubr)
library(viridis)
library(shadowtext)
library(officer)
library(rvg)
library(egg)
library(svglite)
library(ggfittext)
library(gridExtra)
library(grid)
library(lemon)
library(lubridate)
library(latex2exp)

### II. Figure 1: Accelerometry processing and feature extraction pipeline and experimental paradigm

## First, determine candidate profiles to plot as example for figure 1
# Load patient temporal information to get recording duration for patient
patient.temporal.info <- read.csv('../clinical_data/patient_temporal_info.csv')

# Load compiled motion features of all patients (if not already loaded)
if (!exists("all.motion.features")) {
  all.motion.features <- read.csv('../features/all_features.csv')
}

# Filter out cases of high bed activity
active.bed.motion.features <- all.motion.features %>%
  filter(Feature == 'SMA') %>%
  filter(Bed >= 0.135) %>%
  group_by(UPI) %>%
  count() %>%
  left_join(patient.temporal.info %>% select(UPI,HoursDurationAccelRecording))

# Based on number of active bed motion cases and recording duration, we select UPI: '4oCkC1hc' as the example

## Display motion feature values around point of interest for UPI: '4oCkC1hc'
# Load motion features specific to patient '4oCkC1hc' and filter out motion features around point of interest
chosen.case.motion.features <- read.csv('../features/features_4oCkC1hc.csv') %>%
  filter(TimeOfDay >= '14:26:15', TimeOfDay <= '14:26:40') %>%
  mutate(Feature = factor(Feature,levels = c('SMA','HLF_h','HLF_l','FDE','MFR','BPW','WVL'))) %>%
  dplyr::select(-c(UPI,RecordingIdx,HoursFromICUAdmission)) %>%
  pivot_longer(cols = -c(TimeOfDay,Feature),names_to = 'Sensor',values_to = 'value') %>%
  pivot_wider(id_cols = c(Sensor,Feature), names_from = 'TimeOfDay', values_from = 'value') %>%
  arrange(Sensor, Feature)

# Load bed-corrected, imputed motion features specific to patient '4oCkC1hc' and filter out motion features around point of interest
chosen.case.bed.corrected.motion.features <- read.csv('../features/02_bed_corrected_imputed_features/bed_corrected_imputation_1.csv') %>%
  filter(UPI == '4oCkC1hc') %>%
  filter(TimeOfDay >= '14:26:15', TimeOfDay <= '14:26:40') %>%
  mutate(Feature = factor(Feature,levels = c('SMA','HLF_h','HLF_l','FDE','MFR','BPW','WVL'))) %>%
  dplyr::select(-c(UPI,RecordingIdx,HoursFromICUAdmission)) %>%
  pivot_longer(cols = -c(TimeOfDay,Feature),names_to = 'Sensor',values_to = 'value') %>%
  pivot_wider(id_cols = c(Sensor,Feature), names_from = 'TimeOfDay', values_from = 'value') %>%
  arrange(Sensor, Feature)

### III. Figure 2: Discrimination performance of motor function detection models on validation sets

## (a) ROC curves of threshold-level GCSm detection
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

# Compile the significant and non-significant AUC values to determine optimally discriminating ROCs for plotting
opt.GCSm.AUC.df <- rbind(opt.GCSm.signficant.AUC.df,opt.GCSm.nonsignficant.AUC.df)

# Load compiled ROC axes dataframe
compiled.GCSm.roc.axes.df <- read.csv('../results/GCSm_threshold_prediction/compiled_ROC.csv')

# Filter out observation window/threshold combinations in `opt.GCSm.AUC.df`
plot.GCSm.roc.axes.df <- compiled.GCSm.roc.axes.df %>%
  inner_join(opt.GCSm.AUC.df,by=c('Threshold','ObsWindow')) %>%
  arrange(Threshold,ObsWindow)

# Change GCSm labels to proper format for plot
plot.GCSm.roc.axes.df$Threshold <- plyr::mapvalues(plot.GCSm.roc.axes.df$Threshold,
                                                   from = c(
                                                     "GCSm.gt.1",
                                                     "GCSm.gt.2",
                                                     "GCSm.gt.3",
                                                     "GCSm.gt.4",
                                                     "GCSm.gt.5"
                                                   ),
                                                   to = c("GCSm > 1", 
                                                          "GCSm > 2", 
                                                          "GCSm > 3", 
                                                          "GCSm > 4", 
                                                          "GCSm > 5"))

# Fix endpoints to corners of AUC plot
plot.GCSm.roc.axes.df[plot.GCSm.roc.axes.df$FPR == 0,c('meanTPR','medianTPR','lowerTPR','upperTPR')] <- 0
plot.GCSm.roc.axes.df[plot.GCSm.roc.axes.df$FPR == 1,c('meanTPR','medianTPR','lowerTPR','upperTPR')] <- 1

# Use ggplot to visualize optimal ROC curves per threshold
source('functions/plot_ROC.R')
GCSm.roc.curves.plot <- plot.ROC(plot.GCSm.roc.axes.df,axis.text.font.size = 5)

# Create directory for current date and save GCSm ROC plot
dir.create(file.path('../plots',Sys.Date()),showWarnings = F,recursive = T)
ggsave(file.path('../plots',Sys.Date(),'GCSm_ROC.svg'),GCSm.roc.curves.plot,device= svg,units='in',dpi=300,width=6.5,height = 4.61)

# Determine text labels for optimal observation window, AUC, and confidence interval
GCSm.roc.text.df <- plot.GCSm.roc.axes.df %>%
  dplyr::select(Threshold,ObsWindow,meanAUC,lowerAUC,upperAUC) %>%
  distinct() %>%
  mutate(formatted.label = sprintf('Optimal Obs. Window: %s hr \n AUC: %0.2f (%0.2f – %0.2f)',ObsWindow,meanAUC,lowerAUC,upperAUC))

## (b) AUC vs. observation window curves of threshold-level GCSm detection
# Load compiled AUC metrics
compiled.GCSm.AUC.df <- read.csv('../results/GCSm_threshold_prediction/compiled_metrics.csv') %>%
  filter(Metrics == 'AUC')

# Reformat threshold names for figure
compiled.GCSm.AUC.df$Threshold <- plyr::mapvalues(compiled.GCSm.AUC.df$Threshold, from = c("GCSm.gt.1","GCSm.gt.2","GCSm.gt.3","GCSm.gt.4","GCSm.gt.5"), to = c("GCSm > 1","GCSm > 2","GCSm > 3","GCSm > 4","GCSm > 5"))

# Use ggplot to visualize AUC vs. observation window per threshold
source('functions/plot_AUC_v_ObsWindow.R')
GCSm.AUC.curves.plot <- plot.AUC.ObsWindow(compiled.GCSm.AUC.df,ow.cutoff = 30,axis.text.font.size = 5)

# Create directory for current date and save GCSm AUC plot
dir.create(file.path('../plots',Sys.Date()),showWarnings = F,recursive = T)
ggsave(file.path('../plots',Sys.Date(),'GCSm_AUC.svg'),GCSm.AUC.curves.plot,device= svg,units='in',dpi=300,width=6.5,height = 2.75)

### IV. Figure 3: Discrimination performance of functional outcome at hospital discharge prediction models on validation sets

## (a) ROC curves of threshold-level GOSE (discharge) prediction
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

# Load compiled ROC axes dataframe
compiled.GOSE.roc.axes.df <- read.csv('../results/GOSE_threshold_prediction/compiled_ROC.csv')

# Filter out observation window/threshold combinations in `opt.GOSE.AUC.df`
plot.GOSE.roc.axes.df <- compiled.GOSE.roc.axes.df %>%
  inner_join(opt.GOSE.AUC.df,by=c('Threshold','ObsWindow')) %>%
  arrange(Threshold,ObsWindow)

# Change GOSE labels to proper format for plot
plot.GOSE.roc.axes.df$Threshold <- plyr::mapvalues(plot.GOSE.roc.axes.df$Threshold,
                                                   from = c(
                                                     "GOSE.gt.1",
                                                     "GOSE.gt.2",
                                                     "GOSE.gt.3",
                                                     "GOSE.gt.4",
                                                     "GOSE.gt.5"
                                                   ),
                                                   to = c("GOSE > 1", 
                                                          "GOSE > 2", 
                                                          "GOSE > 3", 
                                                          "GOSE > 4", 
                                                          "GOSE > 5"))

# Fix endpoints to corners of AUC plot
plot.GOSE.roc.axes.df[plot.GOSE.roc.axes.df$FPR == 0,c('meanTPR','medianTPR','lowerTPR','upperTPR')] <- 0
plot.GOSE.roc.axes.df[plot.GOSE.roc.axes.df$FPR == 1,c('meanTPR','medianTPR','lowerTPR','upperTPR')] <- 1

# Use ggplot to visualize optimal ROC curves per threshold
source('functions/plot_ROC.R')
GOSE.roc.curves.plot <- plot.ROC(plot.GOSE.roc.axes.df,axis.text.font.size = 5)

# Create directory for current date and save GOSE ROC plot
dir.create(file.path('../plots',Sys.Date()),showWarnings = F,recursive = T)
ggsave(file.path('../plots',Sys.Date(),'GOSE_ROC.svg'),GOSE.roc.curves.plot,device= svg,units='in',dpi=300,width=6.5,height = 4.61)

# Determine text labels for optimal observation window, AUC, and confidence interval
GOSE.roc.text.df <- plot.GOSE.roc.axes.df %>%
  dplyr::select(Threshold,ObsWindow,meanAUC,lowerAUC,upperAUC) %>%
  distinct() %>%
  mutate(formatted.label = sprintf('Optimal Obs. Window: %s hr \n AUC: %0.2f (%0.2f – %0.2f)',ObsWindow,meanAUC,lowerAUC,upperAUC))

## (b) AUC vs. observation window curves of threshold-level GOSE (discharge) prediction
# Load compiled AUC metrics
compiled.GOSE.AUC.df <- read.csv('../results/GOSE_threshold_prediction/compiled_metrics.csv') %>%
  filter(Metrics == 'AUC')

# Reformat threshold names for figure
compiled.GOSE.AUC.df$Threshold <- plyr::mapvalues(compiled.GOSE.AUC.df$Threshold, from = c("GOSE.gt.1","GOSE.gt.2","GOSE.gt.3","GOSE.gt.4","GOSE.gt.5"), to = c("GOSE > 1","GOSE > 2","GOSE > 3","GOSE > 4","GOSE > 5"))

# Use ggplot to visualize AUC vs. observation window per threshold
source('functions/plot_AUC_v_ObsWindow.R')
GOSE.AUC.curves.plot <- plot.AUC.ObsWindow(compiled.GOSE.AUC.df,
                                           ow.cutoff = 6,
                                           ow.units = 'hr',
                                           axis.text.font.size = 5,
                                           auc.min = 0,
                                           auc.max = 1,
                                           step.size = 1)

# Create directory for current date and save GOSE AUC plot
dir.create(file.path('../plots',Sys.Date()),showWarnings = F,recursive = T)
ggsave(file.path('../plots',Sys.Date(),'GOSE_AUC.svg'),GOSE.AUC.curves.plot,device= svg,units='in',dpi=300,width=6.5,height = 2.75)

### V. Figure 4: Probability calibration of optimally discriminating motor function detection and functional outcome prediction models on validation sets

## (a) Calibration curves for AUC-optimal GCSm detection models
# Load GCSm threshold calibration metrics
GCSm.calibration.metrics <- read.csv('../results/GCSm_threshold_prediction/calibration_metrics.csv')

# Isolate optimal configuration indices based on integrated calibration indices (ICI)
optimal.ConfigIdx <- GCSm.calibration.metrics %>%
  dplyr::select(Threshold,ObsWindow,optConfigIdx) %>%
  distinct()

# Reformat threshold names for plot
optimal.ConfigIdx$FormattedThresh <- plyr::mapvalues(optimal.ConfigIdx$Threshold, 
                                                     from = c("GCSm.gt.1","GCSm.gt.2","GCSm.gt.3","GCSm.gt.4","GCSm.gt.5"), 
                                                     to = c("GCSm > 1","GCSm > 2","GCSm > 3","GCSm > 4","GCSm > 5"))

# Create directory for current date to save GCSm feature significance plots
dir.create(file.path('../plots',Sys.Date()),showWarnings = F,recursive = T)

# Iterate through thresholds
for (curr.thresh.idx in 1:nrow(optimal.ConfigIdx)){
  # Identify current threshold, observation window, and configuration index information
  curr.thresh <- optimal.ConfigIdx$Threshold[curr.thresh.idx]
  curr.formatted.thresh <- optimal.ConfigIdx$FormattedThresh[curr.thresh.idx]
  curr.obs.window <- optimal.ConfigIdx$ObsWindow[curr.thresh.idx]
  curr.opt.config.idx <- optimal.ConfigIdx$optConfigIdx[curr.thresh.idx]
  
  # Load compiled predictions of current threshold-observation window combination
  # and filter out predictions of current optimal configuration index
  curr.predictions <-
    read.csv(file.path('../results/GCSm_threshold_prediction',
                       paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
                       paste0(curr.thresh,'_compiled_predictions.csv'))) %>%
    filter(ConfigIdx == curr.opt.config.idx)
  
  # Separate true-positive and true-negative predictions
  true.pos.preds <- curr.predictions %>%
    filter(TrueLabel == 1)
  true.neg.preds <- curr.predictions %>%
    filter(TrueLabel == 0)
  
  # Print out number of distinct positive and negative samples in current configuration
  print(paste(curr.thresh,'Positive n =',nrow(true.pos.preds %>% dplyr::select(-c(Prob,Repeat,Fold)) %>% distinct(UPI,HoursFromICUAdmission))))
  print(paste(curr.thresh,'Negative n =',nrow(true.neg.preds %>% dplyr::select(-c(Prob,Repeat,Fold)) %>% distinct(UPI,HoursFromICUAdmission))))
  
  # Load current calibration curve information
  curr.GCSm.calibration.curve <- read.csv(file.path('../results/GCSm_threshold_prediction',
                                                    paste0(curr.thresh,'_calibration_curve.csv'))) %>%
    drop_na(meanTrueProb)
  
  # Create ggplot object of calibration plot
  calib.plot <- ggplot(data = curr.GCSm.calibration.curve) +
    coord_cartesian(xlim = c(0,1),ylim = c(0,1),expand = T)+
    geom_segment(x = 0, y = 0, xend = 1, yend = 1,alpha=0.5,linetype = "dashed",size=.75/.pt, color = 'gray') +
    geom_ribbon(aes(x = PredProb, ymin = lowerTrueProb, ymax = upperTrueProb), alpha = 0.1,fill='red',linetype = "dotdash",size=.75/.pt,color='black') +
    geom_line(aes(x = PredProb,y = meanTrueProb),size=1.3/.pt,color='red') +
    xlab("") +
    ylab("") +
    ggtitle(curr.formatted.thresh)+
    theme_classic()+
    theme(
      strip.text = element_text(size=7, color = "black",face = 'bold'), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 5, color = "black"),
      axis.text.y = element_text(size = 5, color = "black"),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size = 2/.pt),
      aspect.ratio = 1,
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 7, color = "black",face = 'bold',hjust = 0.5),
      plot.margin=grid::unit(c(0,0,0,0), "mm")
    )
  
  # Save calibration plot of current GCSm threshold
  ggsave(file.path('../plots',Sys.Date(),paste0(curr.thresh,'_calibration.svg')),calib.plot,device= svg,units='in',dpi=300,width=2.12,height = 2.19)
  
  # Create ggplot object of prediction distribution plot
  dist.plot <- ggplot(data = NULL,mapping = aes(x)) +
    geom_histogram(data = true.pos.preds,mapping = aes(x = Prob, y = ..density..),fill = "black", bins = 200) +
    geom_histogram(data = true.neg.preds,mapping = aes(x = Prob, y = -..density..),fill = "black", bins = 200) +
    geom_segment(inherit.aes = F,aes(x = 0, y = 0, xend = 1, yend = 0),size=.75/.pt, color = 'gray') +
    coord_cartesian(xlim = c(0,1), expand = T)+
    scale_y_symmetric(mid = 0) +
    xlab("") +
    ylab("") +
    theme_minimal()+
    theme(
      rect = element_rect(fill = "transparent"),
      strip.text = element_text(size=7, color = "black",face = 'bold'), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 5, color = NA),
      axis.text.y = element_text(size = 5, color = NA),
      axis.ticks.y = element_blank(),
      strip.background = element_blank(),
      panel.border = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin=grid::unit(c(0,0,0,0), "mm")
    )
  
  # Save distribution histogram of current GCSm threshold
  ggsave(file.path('../plots',Sys.Date(),paste0(curr.thresh,'_distribution.svg')),dist.plot,device= svg,units='in',dpi=300,width=2.1,height = .5,bg ="transparent")
}

## (b) Calibration curves for AUC-optimal GOSE (discharge) prediction models
# Load GOSE threshold calibration metrics
GOSE.calibration.metrics <- read.csv('../results/GOSE_threshold_prediction/calibration_metrics.csv')

# Isolate optimal configuration indices based on integrated calibration indices (ICI)
optimal.ConfigIdx <- GOSE.calibration.metrics %>%
  dplyr::select(Threshold,ObsWindow,optConfigIdx) %>%
  distinct()

# Reformat threshold names for plot
optimal.ConfigIdx$FormattedThresh <- plyr::mapvalues(optimal.ConfigIdx$Threshold, 
                                                     from = c("GOSE.gt.1","GOSE.gt.2","GOSE.gt.3","GOSE.gt.4","GOSE.gt.5"), 
                                                     to = c("GOSE > 1","GOSE > 2","GOSE > 3","GOSE > 4","GOSE > 5"))

# Create directory for current date to save GOSE feature significance plots
dir.create(file.path('../plots',Sys.Date()),showWarnings = F,recursive = T)

# Iterate through thresholds
for (curr.thresh.idx in 1:nrow(optimal.ConfigIdx)){
  # Identify current threshold, observation window, and configuration index information
  curr.thresh <- optimal.ConfigIdx$Threshold[curr.thresh.idx]
  curr.formatted.thresh <- optimal.ConfigIdx$FormattedThresh[curr.thresh.idx]
  curr.obs.window <- optimal.ConfigIdx$ObsWindow[curr.thresh.idx]
  curr.opt.config.idx <- optimal.ConfigIdx$optConfigIdx[curr.thresh.idx]
  
  # Load compiled predictions of current threshold-observation window combination
  # and filter out predictions of current optimal configuration index
  curr.predictions <-
    read.csv(file.path('../results/GOSE_threshold_prediction',
                       paste0(sprintf('%05.2f', curr.obs.window), '_h_obs_window'),
                       paste0(curr.thresh,'_compiled_predictions.csv'))) %>%
    filter(ConfigIdx == curr.opt.config.idx)
  
  # Separate true-positive and true-negative predictions
  true.pos.preds <- curr.predictions %>%
    filter(TrueLabel == 1)
  true.neg.preds <- curr.predictions %>%
    filter(TrueLabel == 0)
  
  # Print out number of distinct positive and negative samples in current configuration
  print(paste(curr.thresh,'Positive n =',nrow(true.pos.preds %>% dplyr::select(-c(Prob,Repeat,Fold)) %>% distinct(UPI,HoursFromICUAdmission))))
  print(paste(curr.thresh,'Negative n =',nrow(true.neg.preds %>% dplyr::select(-c(Prob,Repeat,Fold)) %>% distinct(UPI,HoursFromICUAdmission))))
  
  # Load current calibration curve information
  curr.GOSE.calibration.curve <- read.csv(file.path('../results/GOSE_threshold_prediction',
                                                    paste0(curr.thresh,'_calibration_curve.csv'))) %>%
    drop_na(meanTrueProb)
  
  # Create ggplot object of calibration plot
  calib.plot <- ggplot(data = curr.GOSE.calibration.curve) +
    coord_cartesian(xlim = c(0,1),ylim = c(0,1),expand = T)+
    geom_segment(x = 0, y = 0, xend = 1, yend = 1,alpha=0.5,linetype = "dashed",size=.75/.pt, color = 'gray') +
    geom_ribbon(aes(x = PredProb, ymin = lowerTrueProb, ymax = upperTrueProb), alpha = 0.1,fill='red',linetype = "dotdash",size=.75/.pt,color='black') +
    geom_line(aes(x = PredProb,y = meanTrueProb),size=1.3/.pt,color='red') +
    xlab("") +
    ylab("") +
    ggtitle(curr.formatted.thresh)+
    theme_classic()+
    theme(
      strip.text = element_text(size=7, color = "black",face = 'bold'), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 5, color = "black"),
      axis.text.y = element_text(size = 5, color = "black"),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size = 2/.pt),
      aspect.ratio = 1,
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 7, color = "black",face = 'bold',hjust = 0.5),
      plot.margin=grid::unit(c(0,0,0,0), "mm")
    )
  
  # Save calibration plot of current GOSE threshold
  ggsave(file.path('../plots',Sys.Date(),paste0(curr.thresh,'_calibration.svg')),calib.plot,device= svg,units='in',dpi=300,width=2.12,height = 2.19)
  
  # Create ggplot object of prediction distribution plot
  dist.plot <- ggplot(data = NULL,mapping = aes(x)) +
    geom_histogram(data = true.pos.preds,mapping = aes(x = Prob, y = ..density..),fill = "black", bins = 200) +
    geom_histogram(data = true.neg.preds,mapping = aes(x = Prob, y = -..density..),fill = "black", bins = 200) +
    geom_segment(inherit.aes = F,aes(x = 0, y = 0, xend = 1, yend = 0),size=.75/.pt, color = 'gray') +
    coord_cartesian(xlim = c(0,1), expand = T)+
    scale_y_symmetric(mid = 0) +
    xlab("") +
    ylab("") +
    theme_minimal()+
    theme(
      rect = element_rect(fill = "transparent"),
      strip.text = element_text(size=7, color = "black",face = 'bold'), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 5, color = NA),
      axis.text.y = element_text(size = 5, color = NA),
      axis.ticks.y = element_blank(),
      strip.background = element_blank(),
      panel.border = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin=grid::unit(c(0,0,0,0), "mm")
    )
  
  # Save distribution histogram of current GOSE threshold
  ggsave(file.path('../plots',Sys.Date(),paste0(curr.thresh,'_distribution.svg')),dist.plot,device= svg,units='in',dpi=300,width=2.1,height = .5,bg ="transparent")
}

### VI. Figure 5: Feature significance matrices of optimally discriminating motor function detection and functional outcome prediction models

## (a) Feature significance plot for AUC-optimal GCSm detection model (GCSm > 4, Obs. Window = 6 hr)
# Load bootstrapped feature significance value information for AUC-optimal model
GCSm.feature.sig.values <- read.csv('../results/GCSm_threshold_prediction/feature_significance.csv')

# Calculate summary statistics of the mean significance score for each sensor-feature combination
summ.GCSm.feature.sig.values <- GCSm.feature.sig.values %>%
  group_by(Sensor,Feature) %>%
  summarise(sampleMean = mean(meanSignificance,na.rm = T),
            lowerMean = quantile(meanSignificance,0.025,na.rm = T),
            upperMean = quantile(meanSignificance,0.975,na.rm = T)) %>%
  mutate(formattedMean = sprintf('%0.2f \n (%0.2f – %0.2f)',sampleMean,lowerMean,upperMean))

# Remap feature values to formatted feature labels
summ.GCSm.feature.sig.values$Feature <- plyr::mapvalues(
  summ.GCSm.feature.sig.values$Feature,
  from = c("SMA","HLF_h","HLF_l","MFR","FDE","BPW","WVL","PhysActivity"),
  to = c("SMA","HLF (h)","HLF (l)","MFR","FDE","BPW","WVL","PDA")
)

# Create `ggplot` heatmap of feature significance scores
GCSm.feat.sig.matrix <- summ.GCSm.feature.sig.values %>%
  ggplot(aes(x = Sensor,y = Feature,fill = sampleMean))+
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  geom_text(aes(label = formattedMean,color = as.factor(as.integer(sampleMean>0.60))),show.legend = F,size = 5/.pt)+
  guides(fill = guide_colourbar(title = sprintf('Mean significance score'),
                                barwidth = grid::unit(2.75,'inches'),
                                barheight = grid::unit(.15,'inches'),
                                direction="horizontal",
                                title.position = 'top',
                                frame.colour=c("black"),
                                frame.linewidth = 1.5/.pt,
                                title.hjust = .5))+
  xlab(label = "Sensor") +
  ylab(label = "Feature") +
  scale_color_manual(values = c('white','black')) +
  scale_y_discrete(limits = rev(c("PDA","SMA","HLF (h)","HLF (l)","FDE","MFR","BPW","WVL")))+
  scale_x_discrete(limits = c("RE","LE","RW","LW","RA","LA"))+
  theme_classic()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black",angle=90,vjust = 0,hjust = .5),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 8, color = "black",face = 'bold'),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size = 2/.pt),
    legend.position = 'bottom',
    legend.title = element_text(size = 7, color = "black"),
    legend.text=element_text(size=5),
    aspect.ratio = 7/6
  )

# Create directory for current date and save GCSm feature significance plot
dir.create(file.path('../plots',Sys.Date()),showWarnings = F,recursive = T)
ggsave(file.path('../plots',Sys.Date(),'GCSm_feat_sig_matrix.svg'),GCSm.feat.sig.matrix,device= svg,units='in',dpi=300,width=3.15,height = 4.44)

## (b) Feature significance plot for AUC-optimal GOSE detection model (GOSE > 5, Obs. Window = 6 hr)
# Load bootstrapped feature significance value information for AUC-optimal model
GOSE.feature.sig.values <- read.csv('../results/GOSE_threshold_prediction/feature_significance.csv')

# Calculate summary statistics of the mean significance score for each sensor-feature combination
summ.GOSE.feature.sig.values <- GOSE.feature.sig.values %>%
  group_by(Sensor,Feature) %>%
  summarise(sampleMean = mean(meanSignificance,na.rm = T),
            lowerMean = quantile(meanSignificance,0.025,na.rm = T),
            upperMean = quantile(meanSignificance,0.975,na.rm = T)) %>%
  mutate(formattedMean = sprintf('%0.2f \n (%0.2f – %0.2f)',sampleMean,lowerMean,upperMean))

# Remap feature values to formatted feature labels
summ.GOSE.feature.sig.values$Feature <- plyr::mapvalues(
  summ.GOSE.feature.sig.values$Feature,
  from = c("SMA","HLF_h","HLF_l","MFR","FDE","BPW","WVL","PhysActivity"),
  to = c("SMA","HLF (h)","HLF (l)","MFR","FDE","BPW","WVL","PDA")
)

# Create `ggplot` heatmap of feature significance scores
GOSE.feat.sig.matrix <- summ.GOSE.feature.sig.values %>%
  ggplot(aes(x = Sensor,y = Feature,fill = sampleMean))+
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  geom_text(aes(label = formattedMean,color = as.factor(as.integer(sampleMean>0.80))),show.legend = F,size = 5/.pt)+
  guides(fill = guide_colourbar(title = sprintf('Mean significance score'),
                                barwidth = grid::unit(2.75,'inches'),
                                barheight = grid::unit(.15,'inches'),
                                direction="horizontal",
                                title.position = 'top',
                                frame.colour=c("black"),
                                frame.linewidth = 1.5/.pt,
                                title.hjust = .5))+
  xlab(label = "Sensor") +
  ylab(label = "Feature") +
  scale_color_manual(values = c('white','black')) +
  scale_y_discrete(limits = rev(c("PDA","SMA","HLF (h)","HLF (l)","FDE","MFR","BPW","WVL")))+
  scale_x_discrete(limits = c("RE","LE","RW","LW","RA","LA"))+
  theme_classic()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black",angle=90,vjust = 0,hjust = .5),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 8, color = "black",face = 'bold'),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size = 2/.pt),
    legend.position = 'bottom',
    legend.title = element_text(size = 7, color = "black"),
    legend.text=element_text(size=5),
    aspect.ratio = 7/6
  )

# Create directory for current date and save GOSE AUC plot
dir.create(file.path('../plots',Sys.Date()),showWarnings = F,recursive = T)
ggsave(file.path('../plots',Sys.Date(),'GOSE_feat_sig_matrix.svg'),GOSE.feat.sig.matrix,device= svg,units='in',dpi=300,width=3.15,height = 4.44)

### VII. Figure 6: Retrospective case study analysis of accelerometry-based detection of motor function in six patients who experienced relevant transition

## Load case study predictions and associated GCSm data
# Load time limit information
case.study.time.limits <- read.csv('../results/case_study_analysis_v2/case_study_time_limits.csv')

# Compiled case study predictions (6 hr):
compiled.case.study.predictions.06.00hr <- read.csv('../results/case_study_analysis/compiled_case_study_predictions.csv') %>%
  mutate(ObsWindow = 6)

# Compiled case study predictions (27 min):
compiled.case.study.predictions.00.45hr <- read.csv('../results/case_study_analysis_v2/compiled_case_study_predictions.csv') %>%
  mutate(ObsWindow = 0.45)

# Compile both case study prediction sets
compiled.case.study.predictions <- rbind(compiled.case.study.predictions.00.45hr,
                                         compiled.case.study.predictions.06.00hr) %>%
  relocate(ObsWindow,Case,UPI)

# Load neurological assessment data and isolate case study patient cases
case.study.gcs.data <- read.csv('../clinical_data/neurological_assessments.csv') %>%
  filter(UPI %in% case.study.time.limits$UPI) %>%
  drop_na(GCSm) %>%
  left_join(case.study.time.limits %>% dplyr::select(UPI,Case),
            by = 'UPI') %>%
  relocate(Case)

# Transform GCS scores to common axis with probability scores
case.study.gcs.data$TransformedGCSm = plyr::mapvalues(case.study.gcs.data$GCSm,from = c(1,2,3,4,5,6),to = c(.1,.2,.3,.4,2/3,5/6))

# Load patient temporal information of case study patients and create dummy variable for ICU admission timestamp
case.study.patient.temporal.info <- read.csv('../clinical_data/patient_temporal_info.csv') %>%
  filter(UPI %in% case.study.time.limits$UPI) %>%
  mutate(TimeStampICUAdmission = as.POSIXct(paste('1970-01-01',TimeOfDayICUAdmission),format = '%Y-%m-%d %H:%M:%S',tz = 'UTC'))

# Add a new dummy variable (`TimeStamp`) to compiled predictions and GCS information for plotting purposes
compiled.case.study.predictions <- left_join(compiled.case.study.predictions,
                                             case.study.patient.temporal.info %>% 
                                               dplyr::select(UPI,TimeStampICUAdmission),
                                             by = 'UPI') %>%
  mutate(TimeStamp = TimeStampICUAdmission + seconds(round(HoursFromICUAdmission*3600))) %>%
  dplyr::select(-c(TimeOfDay,TimeStampICUAdmission)) %>%
  relocate(TimeStamp,.after = UPI)

case.study.gcs.data <- left_join(case.study.gcs.data,
                                 case.study.patient.temporal.info %>% 
                                   dplyr::select(UPI,TimeStampICUAdmission),
                                 by = 'UPI') %>%
  mutate(TimeStamp = TimeStampICUAdmission + seconds(round(HoursFromICUAdmission*3600))) %>%
  dplyr::select(-c(TimeOfDay,TimeStampICUAdmission)) %>%
  relocate(TimeStamp,.after = UPI)

case.study.time.limits <- left_join(case.study.time.limits,
                                    case.study.patient.temporal.info %>% 
                                      dplyr::select(UPI,TimeStampICUAdmission),
                                    by = 'UPI') %>%
  mutate(startTimeStamp = TimeStampICUAdmission + seconds(round(startHoursFromICUAdmission*3600)),
         endTimeStamp = TimeStampICUAdmission + seconds(round(endHoursFromICUAdmission*3600))) %>%
  dplyr::select(-TimeStampICUAdmission)

# Filter out GCS information within the time limits (plus one more and less) for plotting purposes
filtered.case.study.gcs.data <- data.frame(matrix(ncol=9,nrow=0))
for (curr.UPI in case.study.time.limits$UPI){
  
  #Extract current start and end limits
  curr.startHoursFromICUAdmission <- case.study.time.limits$startHoursFromICUAdmission[case.study.time.limits$UPI == curr.UPI]
  curr.endHoursFromICUAdmission <- case.study.time.limits$endHoursFromICUAdmission[case.study.time.limits$UPI == curr.UPI]
  
  # Identify row indices designating current GCSm selection
  curr.row.indices <- which((case.study.gcs.data$UPI == curr.UPI) &
                              (case.study.gcs.data$HoursFromICUAdmission >= curr.startHoursFromICUAdmission) &
                              (case.study.gcs.data$HoursFromICUAdmission <= curr.endHoursFromICUAdmission))
  
  # Add one to the end of each row index limit
  #curr.row.indices <- c(min(curr.row.indices) - 1,curr.row.indices,max(curr.row.indices) + 1)
  
  # Append GCS data of the desired row indices to the filtered dataframe
  filtered.case.study.gcs.data <- rbind(filtered.case.study.gcs.data,
                                        case.study.gcs.data[curr.row.indices,])
}

## Create prediction trajectory plots with corresponding GCSm information
# Use `ggplot`:
dir.create(file.path('../plots',Sys.Date()),showWarnings = F,recursive = T)

for (curr.UPI in case.study.time.limits$UPI){
  
  curr.UPI.start.limit <- case.study.time.limits$startTimeStamp[case.study.time.limits$UPI == curr.UPI]
  curr.UPI.end.limit <- case.study.time.limits$endTimeStamp[case.study.time.limits$UPI == curr.UPI]
  
  curr.UPI.case.study.predictions <- compiled.case.study.predictions %>%
    filter(UPI == curr.UPI)
  curr.case.label <- case.study.time.limits$Case[case.study.time.limits$UPI == curr.UPI]
  
  curr.UPI.filtered.gcs.data <- filtered.case.study.gcs.data %>% filter(UPI == curr.UPI)
  
  curr.trajectories.plot <- ggplot(data = NULL,mapping = aes(x = TimeStamp)) +
    geom_line(data = curr.UPI.case.study.predictions, mapping = aes(y=meanProb,color=as.factor(ObsWindow)),size = 1.15/.pt) +
    geom_ribbon(data = curr.UPI.case.study.predictions, mapping = aes(ymin = lowerProb,ymax = upperProb,fill=as.factor(ObsWindow)),alpha = 0.2) + 
    geom_hline(yintercept = .5,color='gray',size = .75/.pt,alpha = 0.6) +
    scale_x_datetime("Time of Day \n (Day of ICU Stay)", date_labels = "%H:%M \n (Day %e)",
                     date_breaks = "3 hour")+
    xlab('Days from ICU Admission') +
    ylab('Pr(GCSm > 4)') +
    guides(fill = guide_legend(nrow=1, title.position = 'top'),
           color = guide_legend(nrow=1, title.position = 'top')) +
    scale_fill_discrete(name = 'Observation Window', labels = c('27 min','6 hr')) +
    scale_color_discrete(name = 'Observation Window', labels = c('27 min','6 hr')) +
    geom_vline(data = curr.UPI.filtered.gcs.data, mapping = aes(xintercept = TimeStamp),color='gray',linetype = 'dashed',size = 1/.pt) +
    coord_cartesian(ylim=c(0,1),
                    xlim=c(curr.UPI.start.limit,curr.UPI.end.limit),
                    expand = F)+
    ggtitle(curr.case.label)+
    theme_bw() + 
    theme(panel.grid = element_blank(),
          legend.position = 'bottom',
          legend.title = element_text(size = 7, color = "black"),
          legend.text=element_text(size=5),
          strip.text = element_text(size=7, color = "black",face = 'bold'), 
          axis.text.x = element_text(size = 5, color = "black"),
          axis.text.y = element_text(size = 5, color = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size = 2/.pt),
          axis.title.y = element_text(size = 6, color = "black",face = 'bold'),
          axis.title.x = element_text(size = 6, color = "black",face = 'bold'),
          plot.title = element_text(size = 7, color = "black",face = 'bold',hjust = 0.5),
          plot.margin=grid::unit(c(0,0,0,0), "mm"))
  
  ggsave(file.path('../plots',Sys.Date(),paste0(curr.UPI,'legend_trajectory_plot.svg')),curr.trajectories.plot,device= svg,units='in',dpi=300,width=6,height = 1.5,bg ="transparent")
}

### VIII. Supplementary Figure 1: Discrimination performance of functional outcome at 12 months post discharge prediction models on validation sets

## (a) ROC curves of threshold-level GOSE (12 months) prediction
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

# Load compiled ROC axes dataframe
compiled.GOSE12m.roc.axes.df <- read.csv('../results/GOSE12m_threshold_prediction/compiled_ROC.csv')

# Filter out observation window/threshold combinations in `opt.GOSE12m.AUC.df`
plot.GOSE12m.roc.axes.df <- compiled.GOSE12m.roc.axes.df %>%
  inner_join(opt.GOSE12m.AUC.df,by=c('Threshold','ObsWindow')) %>%
  arrange(Threshold,ObsWindow)

# Change GOSE12m labels to proper format for plot
plot.GOSE12m.roc.axes.df$Threshold <- plyr::mapvalues(plot.GOSE12m.roc.axes.df$Threshold,
                                                      from = c(
                                                        "GOSE12m.gt.1",
                                                        "GOSE12m.gt.2",
                                                        "GOSE12m.gt.3",
                                                        "GOSE12m.gt.4",
                                                        "GOSE12m.gt.5",
                                                        "GOSE12m.gt.6",
                                                        "GOSE12m.gt.7"
                                                      ),
                                                      to = c("GOSE (12m) > 1", 
                                                             "GOSE (12m) > 2", 
                                                             "GOSE (12m) > 3", 
                                                             "GOSE (12m) > 4", 
                                                             "GOSE (12m) > 5", 
                                                             "GOSE (12m) > 6", 
                                                             "GOSE (12m) > 7"))

# Fix endpoints to corners of AUC plot
plot.GOSE12m.roc.axes.df[plot.GOSE12m.roc.axes.df$FPR == 0,c('meanTPR','medianTPR','lowerTPR','upperTPR')] <- 0
plot.GOSE12m.roc.axes.df[plot.GOSE12m.roc.axes.df$FPR == 1,c('meanTPR','medianTPR','lowerTPR','upperTPR')] <- 1

# Use ggplot to visualize optimal ROC curves per threshold
source('functions/plot_ROC.R')
GOSE12m.roc.curves.plot <- plot.ROC(plot.GOSE12m.roc.axes.df,axis.text.font.size = 5,num.col = 4)

# Create directory for current date and save GOSE12m ROC plot
dir.create(file.path('../plots',Sys.Date()),showWarnings = F,recursive = T)
ggsave(file.path('../plots',Sys.Date(),'GOSE12m_ROC.svg'),GOSE12m.roc.curves.plot,device= svg,units='in',dpi=300,width=6.5,height = 3.48)

# Determine text labels for optimal observation window, AUC, and confidence interval
GOSE12m.roc.text.df <- plot.GOSE12m.roc.axes.df %>%
  dplyr::select(Threshold,ObsWindow,meanAUC,lowerAUC,upperAUC) %>%
  distinct() %>%
  mutate(formatted.label = sprintf('Optimal Obs. Window: %s hr \n AUC: %0.2f (%0.2f – %0.2f)',ObsWindow,meanAUC,lowerAUC,upperAUC))

## (b) AUC vs. observation window curves of threshold-level GOSE (12 months) prediction
# Load compiled AUC metrics
compiled.GOSE12m.AUC.df <- read.csv('../results/GOSE12m_threshold_prediction/compiled_metrics.csv') %>%
  filter(Metrics == 'AUC')

# Reformat threshold names for figure
compiled.GOSE12m.AUC.df$Threshold <- plyr::mapvalues(compiled.GOSE12m.AUC.df$Threshold,
                                                     from = c(
                                                       "GOSE12m.gt.1",
                                                       "GOSE12m.gt.2",
                                                       "GOSE12m.gt.3",
                                                       "GOSE12m.gt.4",
                                                       "GOSE12m.gt.5",
                                                       "GOSE12m.gt.6",
                                                       "GOSE12m.gt.7"
                                                     ),
                                                     to = c("GOSE (12m) > 1", 
                                                            "GOSE (12m) > 2", 
                                                            "GOSE (12m) > 3", 
                                                            "GOSE (12m) > 4", 
                                                            "GOSE (12m) > 5", 
                                                            "GOSE (12m) > 6", 
                                                            "GOSE (12m) > 7"))

# Use ggplot to visualize AUC vs. observation window per threshold
source('functions/plot_AUC_v_ObsWindow.R')
GOSE12m.AUC.curves.plot <- plot.AUC.ObsWindow(compiled.GOSE12m.AUC.df,
                                              ow.cutoff = 18,
                                              ow.units = 'hr',
                                              axis.text.font.size = 5,
                                              num.col = 4,
                                              auc.min = 0.1,
                                              auc.max = 0.9,
                                              step.size = 2)

# Create directory for current date and save GOSE12m AUC plot
dir.create(file.path('../plots',Sys.Date()),showWarnings = F,recursive = T)
ggsave(file.path('../plots',Sys.Date(),'GOSE12m_AUC.svg'),GOSE12m.AUC.curves.plot,device= svg,units='in',dpi=300,width=6.5,height = 3)

### IX. Supplementary Figure 2: Correlation matrices of extracted motion features across different sensor placements
# Load correlations of motion features across sensors with associated confidence intervals
sensor.correlations.CI <- read.csv('../summary_statistics/sensor_correlations.csv')

# Create a correlation heatmap with `ggplot`
sensor.correlations.plot <- sensor.correlations.CI %>%
  mutate(Feature = factor(Feature,levels = c("SMA","HLF (h)","HLF (l)","MFR","FDE","BPW","WVL"))) %>%
  ggplot(aes(x = term1,y = term2,fill = meanRho))+
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  geom_text(aes(label = FormattedRho,color = as.factor(as.integer(meanRho>0.60))),show.legend = F,size = 3.5/.pt)+
  guides(fill = guide_colourbar(title = 'Spearmans rank correlation coefficient (p)',
                                barwidth = grid::unit(5.5,'inches'),
                                barheight = grid::unit(.15,'inches'),
                                direction="horizontal",
                                title.position = 'top',
                                frame.colour=c("black"),
                                frame.linewidth = 1.5/.pt,
                                title.hjust = .5) )+
  facet_rep_wrap(~Feature,ncol = 3, nrow = 3, scales='free') +
  scale_color_manual(values = c('white','black')) +
  scale_y_discrete(limits = rev(c("Bed","RE","RW","RA","LE","LW","LA")))+
  scale_x_discrete(limits = c("Bed","RE","RW","RA","LE","LW","LA"))+
  theme_classic()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 4, color = "black"),
    axis.text.y = element_text(size = 4, color = "black",angle=90,vjust = 0,hjust = .5),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size=5, color = "black",face = 'bold'),
    panel.border = element_rect(colour = "black", fill=NA, size = 2/.pt),
    legend.position = 'bottom',
    legend.title = element_text(size = 5, color = "black",face = 'bold'),
    legend.text=element_text(size=4),
    aspect.ratio = 1
  )

# Create directory for current date and save GCSm ROC plot
dir.create(file.path('../plots',Sys.Date()),showWarnings = F,recursive = T)
ggsave(file.path('../plots',Sys.Date(),'sensor_correlation.svg'),sensor.correlations.plot,device= svg,units='in',dpi=300,width=6.5,height = 7.43)

### X. Supplementary Figure 3: Violin plots of extracted motion feature values (30 min observation window), stratified by bilateral sensor placement and GCSm scores
## Prepare motion features corresponding to 30 min before GCS evaluations
# Load GCSm and UPI keys corresponding to 30-min observation window motion features
keys.30min <- read.csv('../features/03_formatted_predictor_matrices/full_matrices/00.50_h_imputation_1_keys.csv')

# Load 30-min observation window motion features (only first imputation)
motion.features.30min <- readRDS('../features/03_formatted_predictor_matrices/full_matrices/00.50_h_imputation_1_full_matrix.rds') %>%
  as.data.frame() %>%
  mutate(UPI = keys.30min$UPI,
         GCSm = factor(keys.30min$GCSm)) %>%
  pivot_longer(cols = -c(UPI,GCSm)) %>%
  mutate(Sensor = sub("\\/.*", "", name),
         WindowsBeforeEvaluation = sub(".*/", "", name),
         Feature = str_match(name, "/\\s*(.*?)\\s*/")[,2]) %>%
  dplyr::select(-name) %>%
  relocate(UPI,GCSm,Sensor,Feature,WindowsBeforeEvaluation,value) %>%
  filter(Feature != 'PhysActivity') %>%
  mutate(Feature = factor(Feature,
                          levels = c("SMA","HLF_h","HLF_l","MFR","FDE","BPW","WVL")))

# Create new variable for bilateral sensor placement
motion.features.30min$Placement <-
  factor(plyr::mapvalues(
    motion.features.30min$Sensor,
    from = c("RE", "LE", "RW", "LW", "RA", "LA"),
    to = c("Elbows", "Elbows", "Wrists", "Wrists", "Ankles", "Ankles")
  ),levels = c("Elbows","Wrists","Ankles"))

## Remove outliers (outlier coefficient factor = 2)
# Identify unique features, placements, and GCSm scores
unique.Features <- unique(motion.features.30min$Feature)
unique.Placements <- unique(motion.features.30min$Placement)
unique.GCSm <- unique(motion.features.30min$GCSm)

# Iterate through unique feature, placement, and GCSm combinations and replacee outliers with NA
for (i in unique.Features){
  for (j in unique.Placements){
    for(k in unique.GCSm){
      curr.idx <- motion.features.30min$Feature == i & 
        motion.features.30min$Placement == j & 
        motion.features.30min$GCSm == k
      outliers <- 
        boxplot.stats(motion.features.30min$value[curr.idx],coef = 2)$out
      motion.features.30min[motion.features.30min$value %in% outliers, "value"] = NA
    }
  }
}

# Delete outlier values
motion.features.30min <- motion.features.30min %>% drop_na(value)

## Produce violin plots of motion features stratified by bilateral placement, feature type, and GCSm
# Create `ggplot` object for violin plots
violin.plots <- motion.features.30min %>%
  ggplot(aes(x=GCSm,y=value,fill=GCSm,group=GCSm)) +
  geom_violin() +
  geom_boxplot(width=0.15, 
               fill="white",
               outlier.shape = NA) +
  facet_grid(featureType ~ placement, 
             switch = "y",
             scales = 'free') +
  theme_classic() + 
  theme(
    strip.text = element_text(size=22, color = "black"), 
    panel.grid.major.y = element_line(colour = "grey70"), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 22, color = "black"),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", fill=NA, size = 2),
    legend.position = "none"
  )

# Calculate p-values for violin plots
unique.Features <- unique(motion.features.30min$Feature)
unique.Placements <- unique(motion.features.30min$Placement)
wilcox.p.value.df <- as.data.frame(matrix(ncol = 6,nrow=0))
for (i in unique.Features){
  print(paste(i,"started"))
  for (j in unique.Placements){
    print(paste(j,"started"))
    curr.idx <- motion.features.30min$Feature == i & 
      motion.features.30min$Placement == j
    curr.test <- compare_means(value~GCSm,motion.features.30min[curr.idx,],ref.group = '.all.')
    wilcox.p.value.df <- rbind(wilcox.p.value.df,data.frame(Feature = i, 
                                                            Placement = j, 
                                                            GCSm = curr.test$group2, 
                                                            p = curr.test$p, 
                                                            p.adj = curr.test$p.adj,
                                                            p.signif = curr.test$p.signif))
  }
}

### XI. Supplementary Figure 4: Mean motion feature trajectories in the six hours preceding GCSm evaluation, stratified by GCSm scores and bilateral sensor placement

## Prepare calculated trajectory mean information with 95% confidence intervals
# Load calculated trajectory mean information
trajectory.means.CI <- read.csv('../summary_statistics/feature_mean_trajectories.csv')

# Rename HLF feature names for plot and order feature names for the plot
trajectory.means.CI$Feature <- plyr::mapvalues(trajectory.means.CI$Feature,
                                               from = c('HLF_h','HLF_l'), 
                                               to = c('HLF (h)','HLF (l)'))
trajectory.means.CI$Feature <- factor(trajectory.means.CI$Feature,levels = c("SMA","HLF (h)","HLF (l)","MFR","FDE","BPW","WVL"))

## Produce mean line plots of motion features stratified by extremity, feature type, and GCSm
# Create `ggplot` object for the trajectory plots
GCSm.trajectory.means.plots <- trajectory.means.CI %>%
  mutate(Placement = factor(Placement,levels = c('Elbows','Wrists','Ankles'))) %>%
  ggplot(aes(x = HoursBeforeEvaluation, y= meanValues)) +
  scale_x_reverse(limits=c(6,0),
                  breaks = 0:6) +
  facet_grid(rows = vars(Feature), 
             cols = vars(Placement), 
             scales = "free",
             switch = "y") +
  xlab('Hours before GCSm Evaluation')+
  geom_ribbon(aes(ymin = lowerValues, ymax = upperValues, fill = factor(GCSm)),alpha = 0.2) +
  coord_cartesian(expand = F) +
  guides(fill = guide_legend(title = 'GCSm'),
         color = guide_legend(title = 'GCSm')) +
  geom_vline(xintercept = 0, size = 1, color = 'firebrick') +
  geom_line(aes(color = factor(GCSm)),size=.75/.pt)+
  theme_bw() +
  theme(
    strip.text = element_text(size=7, color = "black",face = 'bold'), 
    strip.placement = "outside",
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 5, color = "black"),
    axis.text.y = element_text(size = 5, color = "black"),
    axis.title.x = element_text(size = 7, color = "black",face = 'bold'),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    panel.spacing = unit(1/.pt, "lines"),
    legend.position = 'bottom',
    legend.title = element_text(size = 7, color = "black",face = 'bold'),
    legend.text=element_text(size=6),
    legend.key.size = unit(1/.pt,"line")
  )

# Create directory for current date and save feature means trajectory plots
dir.create(file.path('../plots',Sys.Date()),showWarnings = F,recursive = T)
ggsave(file.path('../plots',Sys.Date(),'trajectory_means.svg'),GCSm.trajectory.means.plots,device= svg,units='in',dpi=300,width=7.5,height = 6)

### XII. Supplementary Figure 5: Trajectories of motor component scores of the Glasgow Coma Scale (GCSm) of each study participant during ICU stay
## Load relevant clinical information
# Load automatically extracted GCS labels with indicator labels:
gcs.data <- read.csv('../clinical_data/neurological_assessments.csv')

# Load patient temporal information
patient.temporal.info <- read.csv('../clinical_data/patient_temporal_info.csv')

# Extract list of unique UPIs and partition into 3 pages
unique.UPIs <- sort(unique(patient.temporal.info$UPI))
partition.UPIs <- split(unique.UPIs, ceiling(seq_along(unique.UPIs)/24))
partition.patIdx <- split(seq_along(unique.UPIs),ceiling(seq_along(unique.UPIs)/24))

## Produce GCSm trajectory and concurrent accelerometry recording indicator plots for partitions of patients
# Create directory for current date to save GCSm trajectory plots
dir.create(file.path('../plots',Sys.Date()),showWarnings = F,recursive = T)

# Iterate through UPI partitions
for (partition.idx in 1:length(partition.UPIs)){
  # Extract current partition UPIs and patient indices
  curr.partition.UPIs <- partition.UPIs[[partition.idx]]
  curr.partition.patIdxs <- partition.patIdx[[partition.idx]]
  
  # Filter out UPIs in current partition and non-missing GCSm values
  curr.partition.GCS.data <- gcs.data %>%
    filter(UPI %in% curr.partition.UPIs) %>%
    drop_na(GCSm) %>% 
    left_join(data.frame(UPI = curr.partition.UPIs, 
                         patIdx = curr.partition.patIdxs),
              by = c('UPI'))
  
  # With `ggplot`, produce trajectory of GCSm scores for current partition
  curr.partition.GCS.plot <- curr.partition.GCS.data %>%
    ggplot(aes(x = HoursFromICUAdmission/24, y = GCSm)) +
    geom_rect(inherit.aes = FALSE, 
              data = patient.temporal.info %>%
                filter(UPI %in% curr.partition.UPIs) %>%
                left_join(data.frame(UPI = curr.partition.UPIs, 
                                     patIdx = curr.partition.patIdxs),
                          by = c('UPI')),
              aes(xmin=HoursFromICUAdmissionAccelRecordingStart/24,
                  xmax=HoursFromICUAdmissionAccelRecordingEnd/24,
                  ymin=0, ymax=Inf),
              fill='red',
              alpha = .2) +
    geom_line(size = .75/.pt,) +
    geom_point(size = .85/.pt) +
    ylab('GCSm')+
    xlab('Time from Admission (days)')+
    coord_cartesian(ylim=c(1, 6)) +
    scale_y_continuous(breaks = 1:6)+
    scale_x_continuous(expand = c(0,0))+
    facet_wrap(~patIdx,
               scales = 'free_x',
               ncol = 2) +
    theme_bw() + 
    theme(
      axis.text.x = element_text(size = 5, color = "black"),
      axis.text.y = element_text(size = 5, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(size=7, color = "black",face = 'bold'),
      axis.title.x = element_text(size = 8, color = "black",face = 'bold'),
      axis.title.y = element_text(size = 8, color = "black",face = 'bold'),
      panel.border = element_rect(colour = "black", fill=NA, size = 2/.pt),
      plot.margin=grid::unit(c(0,0,0,0), "mm")
    )
  
  # Save current partition trajectory
  ggsave(file.path('../plots',Sys.Date(),paste0('GCSm_trajectory_partition',partition.idx,'.svg')),curr.partition.GCS.plot,device= svg,units='in',dpi=300,width=6.5,height = 10)
}

### XIII. Supplementary Figure 7: Percentages of missing, static, and dynamic accelerometry data by time of day of recording and sensor placement
# Load motion features prior to imputation and bed movement correction
all.motion.features <- read_csv('../features/all_features.csv')
n <- length(unique(all.motion.features$UPI))

# Filter dataset for just SMA
all.SMA.features <-all.motion.features %>% 
  filter(Feature=="SMA") %>% 
  pivot_longer(cols = -c(Feature,UPI,RecordingIdx,HoursFromICUAdmission,TimeOfDay), names_to="Sensor") %>% 
  mutate(na.indicator = as.numeric(is.na(value)), TimeOfDay = format(TimeOfDay,tz="UTC",format ="%H:%M:%S"))

# Calculate missingness by time of day
miss.NM.TOD <- all.SMA.features %>% 
  mutate(no.motion.indicator = as.numeric(value < 0.135)) %>%
  group_by(Sensor,TimeOfDay) %>%
  summarise(n = n(),
            num.Missing = sum(na.indicator), 
            perc.Missing = sum(na.indicator)/n(), 
            num.NM = sum(no.motion.indicator,na.rm = TRUE), 
            perc.NM = sum(no.motion.indicator,na.rm = TRUE)/n())

# `ggplot`: produce missing and static activity percentage line plots by time of day
time.breaks <- unique(miss.NM.TOD$TimeOfDay)[1] + seq(0,86400,by=21600)
time.labels <- c("00:00","06:00","12:00","18:00","24:00")
miss.NM.curves <- miss.NM.TOD %>% 
  filter(sensor != 'Bed') %>%
  mutate(sensor = factor(sensor,levels=c("RE","LE","RW","LW","RA","LA"))) %>%
  ggplot(aes(x = TimeOfDay)) +
  geom_area(aes(y = perc.Missing),fill="grey93") +
  geom_ribbon(aes(ymin = perc.Missing, ymax = perc.Missing + perc.NM), fill = 'lightcyan', alpha = 0.5) +
  geom_ribbon(aes(ymin = perc.Missing + perc.NM, ymax = 1), fill = 'darkolivegreen1', alpha = 0.5) +
  geom_line(aes(y = perc.Missing),color = "red") + 
  geom_line(aes(y = perc.Missing + perc.NM),color="blue") + 
  facet_wrap(~sensor,ncol = 2) + 
  scale_x_datetime(breaks = time.breaks, labels = time.labels) + 
  scale_y_continuous(labels = scales::percent) +
  xlab('Time of Day')+
  theme_classic() +
  theme(
    strip.text = element_text(size=22, color = "black"), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 18, color = "black",angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 22, color = "black"),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size = 1),
    aspect.ratio = 1
  )

### XIV. Supplementary Figure 8: Count histograms of accelerometry recording information
## Load patient temporal information
patient.temporal.info <- read.csv('../clinical_data/patient_temporal_info.csv')

## Format variables into neat dataframe for `ggplot`
histogram.df <- rbind(data.frame(value = as.numeric(patient.temporal.info$HoursDurationAccelRecording), 
                                 type = 'Recording Duration (h)'),
                      data.frame(value = patient.temporal.info$HoursFromICUAdmissionAccelRecordingStart/24, 
                                 type = 'NCCU Admission to Recording (d)'),
                      data.frame(value = patient.temporal.info$HoursDurationAccelRecording/(patient.temporal.info$DaysInICU*24), 
                                 type = 'Proportion of NCCU Stay Recorded'),
                      data.frame(value = patient.temporal.info$HoursFromICUAdmissionAccelRecordingStart/(patient.temporal.info$DaysInICU*24), 
                                 type = 'Proportion of NCCU Stay before Recording'))

## `ggplot`: produce characteristic histograms
supplementary.histograms <- histogram.df %>%
  mutate(type = factor(type, levels = c('Recording Duration (h)',
                                        'NCCU Admission to Recording (d)',
                                        'Proportion of NCCU Stay Recorded',
                                        'Proportion of NCCU Stay before Recording'))) %>%
  ggplot(aes(x = value)) +
  geom_histogram(colour="black", 
                 bins = 25,
                 fill="grey60") +
  facet_wrap(~type,
             scales = 'free',
             strip.position = "bottom") +
  coord_cartesian(ylim = c(0, 23)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size=14, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.placement = 'outside',
    strip.background = element_blank(),
    aspect.ratio = 1,
    panel.spacing = unit(2, "lines")
  )