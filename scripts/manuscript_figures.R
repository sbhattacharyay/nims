#### Plot figures for the manuscript ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# University of Cambridge
# Johns Hopkins University
# email address: sb2406@cam.ac.uk
#
### Contents:
# I. Initialization
# II. Figure 2: Missingness and activity curves by time of day
# III. Figure 3: Violin plots of motion features
# IV. Figure 4: Smoothed feature means vs. time before recording, stratified by GCSm
# V. Figure 5: Heatmaps of standardized mean absolute LOL coefficients
# VI. Figure 6: ROC and PR curves of the best-performing models
# VII. Figure 7: Mean predicted probability by observation window
# VIII. Supplementary Figure 2: Plot complete GCS trajectory of patients in the study
# IX. Supplementary Figure 3: Plot histograms characterizing recording duration and start time of recording
# X. Supplementary Figure 5: Plot hyperparameter optimization tuning results
# XI. Supplementary Figure 6: Plot probability calibration curves

### I. Initialization
## Import necessary packages
library(tidyverse)
library(ggpubr)
library(viridis)
library(shadowtext)

## Load clinical patient data
source('./functions/load_patient_clinical_data.R')
patient.clinical.data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv') %>% arrange(AccelPatientNo_) %>% mutate(ptIdx = 1:nrow(.))

### II. Figure 2: Missingness and activity curves by time of day
## Load motion features prior to imputation and bed movement correction
source('./functions/get_motion_features.R')
if (!exists("all_motion_features")) {
  all_sensors <-
    readMat('../all_motion_feature_data/01_features/complete_sensor_data.mat')$sensors
  all_motion_features <- do.call(rbind, all_sensors)
  featureLabels <- read.csv('../all_motion_feature_data/01_features/feature_names.csv',header = FALSE)
  all_mf_times <- read.csv('../all_motion_feature_data/01_features/indexed_times.csv') %>% mutate(times = as.POSIXct(times,tz="America/New_York",format ="%d-%b-%Y %H:%M:%S"))
}
n <- length(all_motion_features)/length(featureLabels)

## Recode all missing values to NA, find "totally missing" data-streams, and store all feature values in single dataframe
source('./functions/mf_to_dataframe.R')
out.MF2DF <- mf_to_dataframe(all_motion_features,n,verbose = TRUE) 
completeFeatureSet <- out.MF2DF[[1]]
totallyMissingSet <- out.MF2DF[[2]]

## Add corresponding timestamps to compiled motion feature dataset
featureLabels <- unlist(featureLabels[1,])
complete_timestamps <- as.data.frame(matrix(ncol = 3, nrow = 0))
for (i in 1:n){
  curr_timestamps <- all_mf_times %>% filter(ptIdx == i)
  curr_timestamps$timeCount <- seq(nrow(curr_timestamps))
  complete_timestamps <- rbind(complete_timestamps,curr_timestamps)
}
complete_timestamps <- rename(complete_timestamps,timeStamps = times)
completeFeatureSet <- inner_join(completeFeatureSet,complete_timestamps,by=c("ptIdx","timeCount"))
rm(all_motion_features,all_sensors,out.MF2DF,complete_timestamps,all_mf_times)
gc()

## Filter dataset for just SMA
completeSMASet <-completeFeatureSet %>% 
  filter(featureType=="sma") %>% 
  pivot_longer(cols = -c(featureType,ptIdx,timeCount,timeStamps), names_to="sensor") %>% 
  mutate(na.indicator = as.numeric(is.na(value)), timeValue = format(timeStamps,tz="America/New_York",format ="%H:%M:%S")) %>%
  mutate(time_of_day = as.POSIXct(timeValue,tz="America/New_York",format ="%H:%M:%S"))

## Calculate missingness by time of day
miss.NM.TOD <- completeSMASet %>% 
  mutate(no.motion.indicator = as.numeric(value < 0.135)) %>%
  group_by(sensor,timeValue) %>%
  summarise(n = n(),
            num.Missing = sum(na.indicator), 
            perc.Missing = sum(na.indicator)/n(), 
            num.NM = sum(no.motion.indicator,na.rm = TRUE), 
            perc.NM = sum(no.motion.indicator,na.rm = TRUE)/n()) %>%
  mutate(time_of_day = as.POSIXct(timeValue,tz="America/New_York",format ="%H:%M:%S"))

## ggplot: produce missing and static activity pewrcentage line plots by time of day
time.breaks <- unique(miss.NM.TOD$time_of_day)[1] + seq(0,86400,by=21600)
time.labels <- c("00:00","06:00","12:00","18:00","24:00")
miss.NM.curves <- miss.NM.TOD %>% 
  filter(sensor != 'Bed') %>%
  mutate(sensor = factor(sensor,levels=c("RE","LE","RW","LW","RA","LA"))) %>%
  ggplot(aes(x = time_of_day)) +
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

### III. Figure 3: Violin plots of motion features
## Load motion features corresponding to 30 min before GCS evaluations (only first imputation)
long.form.motor.train.df <- readRDS('../all_motion_feature_data/04_formatted_predictor_matrices/imp1/detection_window_0.5/long_form_motor_train_values.rds') %>% mutate(partition = 'train')
long.form.motor.test.df <- readRDS('../all_motion_feature_data/04_formatted_predictor_matrices/imp1/detection_window_0.5/long_form_motor_test_values.rds') %>% mutate(partition = 'test')
long.form.30min <- rbind(long.form.motor.train.df,long.form.motor.test.df)
rm(long.form.motor.train.df,long.form.motor.test.df)
gc()
long.form.30min <- long.form.30min %>%
  mutate(GCSm = factor(GCSm), placement = plyr::mapvalues(sensor, from = c("RE","LE","RW","LW","RA","LA"), to = c("Elbows","Elbows","Wrists","Wrists","Ankles","Ankles"))) %>%
  mutate(placement = factor(placement,levels = c("Elbows","Wrists","Ankles"))) %>% 
  mutate(GCSm = factor(GCSm),
         featureType = plyr::mapvalues(featureType, 
                                       from = c("band_power","freq_entropy","freq_pairs1","freq_pairs2","med_freq","sma","wavelets"),
                                       to = c("BPW","FDE","HLF (h)","HLF (l)","MFR","SMA","WVL"))) %>%
  mutate(featureType = relevel(featureType,"SMA","HLF (h)","HLF (l)","MFR","FDE","BPW","WVL"))

## Remove outliers from long form dataframe of motion features
unique.feature.types <- unique(long.form.30min$featureType)
unique.placements <- unique(long.form.30min$placement)
unique.GCSm <- unique(long.form.30min$GCSm)
for (i in unique.feature.types){
  for (j in unique.placements){
    for(k in unique.GCSm){
      curr.idx <- long.form.30min$featureType == i & 
        long.form.30min$placement == j & 
        long.form.30min$GCSm == k
      outliers <- 
        boxplot.stats(long.form.30min$value[curr.idx],coef = 2)$out
      long.form.30min[long.form.30min$value %in% outliers, "value"] = NA
    }
  }
}
long.form.30min <- long.form.30min %>% drop_na(value)

## ggplot: produce violin plots of motion features stratified by extremity, feature type, and GCSm
violin.plots <- long.form.30min %>%
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

## Calculate p-values for violin plots
unique.feature.types <- unique(long.form.30min$featureType)
unique.placements <- unique(long.form.30min$placement)
wilcox.p.value.df <- as.data.frame(matrix(ncol = 6,nrow=0))
for (i in unique.feature.types){
  print(paste(i,"started"))
  for (j in unique.placements){
    print(paste(j,"started"))
    curr.idx <- long.form.30min$featureType == i & 
      long.form.30min$placement == j
    curr.test <- compare_means(value~GCSm,long.form.30min[curr.idx,],ref.group = '.all.')
    wilcox.p.value.df <- rbind(wilcox.p.value.df,data.frame(feature.type = i, 
                                                            placement = j, 
                                                            GCSm = curr.test$group2, 
                                                            p = curr.test$p, 
                                                            p.adj = curr.test$p.adj,
                                                            p.signif = curr.test$p.signif))
  }
}

### IV. Figure 4: Smoothed feature means vs. time before recording, stratified by GCSm
## Load automatically extracted GCS data and create markers for inclusion in the figure
gcs.data <- read.csv('../clinical_data/clean_auto_GCS_table.csv') %>% select(-X, -Coincide.with.Accel.Recording) %>% mutate(TakenInstant = as.POSIXct(TakenInstant, tz = "America/New_York")) %>% mutate(TakenDay = as.Date(TakenInstant, tz = "America/New_York")) %>% arrange(AccelPatientNo_,TakenInstant)
gcs.data$keep.detection <- FALSE

## Load accelerometry missingness information:
missing.time.info <- read_xlsx('../all_motion_feature_data/MissingPercentTable.xlsx',.name_repair = "universal") %>% 
  mutate(Start.Timestamp = as.POSIXct(Start.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"), End.Timestamp = as.POSIXct(End.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York")) %>% 
  mutate(Recording.Duration = End.Timestamp - Start.Timestamp) %>%
  rename(AccelPatientNo_=Accel.Patient.No.) %>%
  mutate(Bed.miss.hours = Bed*Recording.Duration,LA.miss.hours = LA*Recording.Duration,LE.miss.hours = LE*Recording.Duration,LW.miss.hours = LW*Recording.Duration,RA.miss.hours = RA*Recording.Duration,RE.miss.hours = RE*Recording.Duration,RW.miss.hours = RW*Recording.Duration)

## Load long-form NA-free complete feature set
na.free.complete.feature.set.long <- readRDS('../all_motion_feature_data/01_features/long_feature_set.rds') %>% drop_na(value)
na.free.complete.feature.set.long$keep.detection <- FALSE

## Find GCS data that corresponds to figure 4 visualization
for (i in 1:nrow(missing.time.info)){
  curr.start.time <- missing.time.info$Start.Timestamp[i]
  curr.end.time <- missing.time.info$End.Timestamp[i]
  keep.detection.idx <- which(gcs.data$AccelPatientNo_ == missing.time.info$AccelPatientNo_[i] & (gcs.data$TakenInstant>=curr.start.time & gcs.data$TakenInstant<=curr.end.time+(6*60*60)))
  gcs.data$keep.detection[keep.detection.idx] <- TRUE
}

## Extract and store relevant motion feature values for Figure 4
detection.plot.df <- as.data.frame(matrix(nrow = 0,ncol = 10))
detection.gcs.data <- gcs.data %>% filter(keep.detection == TRUE)
for (i in 1:nrow(detection.gcs.data)){
  print(paste("Detection row",i,"out of",nrow(detection.gcs.data),"started"))
  filt.feat.set <- na.free.complete.feature.set.long %>% filter(AccelPatientNo_ == detection.gcs.data$AccelPatientNo_[i] & timeStamps >= detection.gcs.data$TakenInstant[i] - (6*60*60) & timeStamps <= detection.gcs.data$TakenInstant[i])
  if (nrow(filt.feat.set) == 0){
    print(paste("Detection row",i,"out of",nrow(detection.gcs.data),"skipped"))
    next
  }
  time.seq <- seq(from = detection.gcs.data$TakenInstant[i] - (6*60*60),to = detection.gcs.data$TakenInstant[i],by= 5)
  time.seq.idx <- 1:length(time.sequence)
  merged.filt.feat.set <- left_join(filt.feat.set,data.frame(time.seq,time.seq.idx), by = c("timeStamps" = "time.seq")) %>% 
    mutate(GCST = detection.gcs.data$Glasgow.Coma.Scale.Score[i], GCSm = detection.gcs.data$Best.Motor.Response[i], GCSv = detection.gcs.data$Best.Verbal.Response[i], GCSe = detection.gcs.data$Eye.Opening[i])
  detection.plot.df <- rbind(detection.plot.df, merged.filt.feat.set)
  names(detection.plot.df) <- names(merged.filt.feat.set)
  print(paste("Detection row",i,"out of",nrow(detection.gcs.data),"completed"))
}
saveRDS(object = detection.plot.df, file = '../all_motion_feature_data/detection_plot_df.rds')

## Load filtered, Figure 4-specific motion features
detection.plot.df <- readRDS('../all_motion_feature_data/detection_plot_df.rds') %>% 
  dplyr::select(-c(keep.detection,GCST,GCSv,GCSe))

## Group filtered motion features by GCSm and calculate and save mean feature values
GCSm.detection.grouping <- detection.plot.df %>% drop_na(GCSm) %>% 
  group_by(GCSm,featureType,sensor,time.seq.idx) %>% 
  summarise(mean.feat.value = mean(value), sd.feat.value = sd(value), lower.bound = max((mean(value) - sd(value)),0), upper.bound = mean(value) + sd(value))
saveRDS(object = GCSm.detection.grouping, file = '../all_motion_feature_data/gcsm_detection_group_df.rds')

## Load stratified mean feature values and clean dataframe for plotting
GCSm.detection.grouping <- readRDS('../all_motion_feature_data/gcsm_detection_group_df.rds') %>% 
  filter(sensor != 'Bed') %>%
  mutate(sensor = as.factor(sensor), GCSm = as.factor(GCSm)) %>%
  mutate(placement = plyr::mapvalues(sensor, from = c("RE","LE","RW","LW","RA","LA"), to = c("Elbows","Elbows","Wrists","Wrists","Ankles","Ankles"))) %>%
  mutate(placement = factor(placement,levels = c("Elbows","Wrists","Ankles"))) %>% 
  mutate(featureType = plyr::mapvalues(featureType, 
                                       from = c("band_power","freq_entropy","freq_pairs1","freq_pairs2","med_freq","sma","wavelets"),
                                       to = c("BPW","FDE","HLF (h)","HLF (l)","MFR","SMA","WVL"))) %>%
  mutate(featureType = factor(featureType,levels = c("SMA","HLF (h)","HLF (l)","MFR","FDE","BPW","WVL"))) %>%
  group_by(GCSm,featureType,placement,time.seq.idx) %>%
  summarise(mean.feat.value = mean(mean.feat.value))

## Load moving average function and apply to filtered mean motion feature curves
source('./functions/movingAverage.R')
GCSm.detection.grouping$filt.mean.feat.values <- NA
for (curr.feat.type in unique(GCSm.detection.grouping$featureType)){
  for (curr.placement in unique(GCSm.detection.grouping$placement)){
    for (curr.GCSm in unique(GCSm.detection.grouping$GCSm)){
      row.idx <- GCSm.detection.grouping$featureType == curr.feat.type & GCSm.detection.grouping$placement ==  curr.placement & GCSm.detection.grouping$GCSm ==  curr.GCSm
      curr.time.seq.idx <- GCSm.detection.grouping$time.seq.idx[row.idx]
      curr.mean.feat.value <- GCSm.detection.grouping$mean.feat.value[row.idx]
      curr.df <- data.frame(curr.time.seq.idx, curr.mean.feat.value) %>% arrange(curr.time.seq.idx)
      GCSm.detection.grouping$filt.mean.feat.values[row.idx] <- movingAverage(curr.df$curr.mean.feat.value, n=(10*12+1),fill.type = NA)
    }
  }
}
GCSm.detection.grouping <-
  left_join(GCSm.detection.grouping,
            data.frame(
              time.seq.idx = 1:4321,
              time.before.gcs = rev(seq(
                0,
                6,
                by = 1/60/12,
              ))
            ),by = "time.seq.idx")

## ggplot: produce mean line plots of motion features stratified by extremity, feature type, and GCSm
GCSm.mean.feature.plots <- GCSm.detection.grouping %>%
  ggplot(aes(x = time.before.gcs, y= filt.mean.feat.values)) +
  scale_x_reverse(limits=c(6,0)) +
  facet_grid(rows = vars(featureType), 
             cols = vars(placement), 
             scales = "free",
             switch = "y") +
  xlab('Hours before GCSm Evaluation')+
  geom_vline(xintercept = 6, size = 1,linetype='dotdash',color='slategray3') +
  geom_vline(xintercept = 3, size = 1,linetype='dotdash',color='slategray3') +
  geom_vline(xintercept = 1, size = 1,linetype='dotdash',color='slategray3') +
  geom_vline(xintercept = .5, size = 1,linetype='dotdash',color='slategray3') +
  geom_vline(xintercept = 0, size = 1, color = 'firebrick') +
  geom_line(aes(color = GCSm),size=.75)+
  theme_bw() +
  theme(
    strip.text = element_text(size=22, color = "black"), 
    strip.placement = "outside",
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 22, color = "black"),
    axis.title.y = element_blank(),
    strip.background = element_blank(),
    legend.position = 'bottom',
    legend.title = element_text(size = 18, color = "black"),
    legend.text=element_text(size=18)
  )

### V. Figure 5: Heatmaps of standardized mean absolute LOL coefficients
## Load dataframe of compiled LOL coefficients and average the first three coefficients
lol.coeff.df <- readRDS('../all_motion_feature_data/04_formatted_predictor_matrices/compiled_lol_coefficients.rds')
lol.coeff.df$coeff <- (abs(lol.coeff.df$first.lol.coeff) + abs(lol.coeff.df$second.lol.coeff) + abs(lol.coeff.df$third.lol.coeff))/3
lol.coeff.df$extrem <- str_sub(lol.coeff.df$sensor,2,2)
lol.coeff.df$side <- str_sub(lol.coeff.df$sensor,1,1)
lol.coeff.df$upper.lower <- plyr::mapvalues(lol.coeff.df$sensor,
                                            from = c("RE","LE","RW","LW","RA","LA"),
                                            to = c('U','U','U','U','L','L'))
lol.coeff.df <- lol.coeff.df[,-c(3,4,5)]

## Group LOL coefficients by commmon variables and average across imputations
pairwise.group.lol.coeff <- lol.coeff.df %>%
  group_by(observation.window, imp, sensor, feature.type) %>%
  summarise(abs.mean.lol.coeff = mean(coeff),
            abs.sd.lol.coeff = sd(coeff)) %>%
  group_by(observation.window, sensor, feature.type) %>%
  summarise(pooled.abs.mean.lol.coeff = mean(abs.mean.lol.coeff),
            pooled.abs.sd.lol.coeff = sqrt( (1/n())*sum(abs.sd.lol.coeff^2)  + ((n()+1)/(n()))*(sd(abs.mean.lol.coeff)^2)))

total.group.lol.coeff <- lol.coeff.df %>%
  group_by(observation.window, imp) %>%
  summarise(abs.mean.lol.coeff = mean(coeff),
            abs.sd.lol.coeff = sd(coeff)) %>%
  group_by(observation.window) %>%
  summarise(pooled.abs.mean.lol.coeff = mean(abs.mean.lol.coeff),
            pooled.abs.sd.lol.coeff = sqrt( (1/n())*sum(abs.sd.lol.coeff^2)  + ((n()+1)/(n()))*(sd(abs.mean.lol.coeff)^2))) %>%
  mutate(sensor = "All",feature.type = "All")

sensors.group.lol.coeff <- lol.coeff.df %>%
  group_by(observation.window, imp, feature.type) %>%
  summarise(abs.mean.lol.coeff = mean(coeff),
            abs.sd.lol.coeff = sd(coeff)) %>%
  group_by(observation.window, feature.type) %>%
  summarise(pooled.abs.mean.lol.coeff = mean(abs.mean.lol.coeff),
            pooled.abs.sd.lol.coeff = sqrt( (1/n())*sum(abs.sd.lol.coeff^2)  + ((n()+1)/(n()))*(sd(abs.mean.lol.coeff)^2))) %>%
  mutate(sensor = "All")

features.group.lol.coeff <- lol.coeff.df %>%
  group_by(observation.window, imp, sensor) %>%
  summarise(abs.mean.lol.coeff = mean(coeff),
            abs.sd.lol.coeff = sd(coeff)) %>%
  group_by(observation.window, sensor) %>%
  summarise(pooled.abs.mean.lol.coeff = mean(abs.mean.lol.coeff),
            pooled.abs.sd.lol.coeff = sqrt( (1/n())*sum(abs.sd.lol.coeff^2)  + ((n()+1)/(n()))*(sd(abs.mean.lol.coeff)^2))) %>%
  mutate(feature.type = "All")

# Pool all resultant averaged dataframes together for plotting
pooled.lol.coeff <- rbind(pairwise.group.lol.coeff, total.group.lol.coeff, sensors.group.lol.coeff, features.group.lol.coeff)

## ggplot: feature significance heatmaps
feat.sig.matrix <- pooled.lol.coeff %>%
  mutate(observation.window = plyr::mapvalues(as.factor(observation.window), 
                                              from = c('0.5','1','3','6'), 
                                              to = c("0.5 hr window", "1 hr window","3 hr window","6 hr window"))) %>%
  ggplot(aes(x = sensor,y = feature.type,fill = pooled.abs.mean.lol.coeff))+
  geom_tile() +
  geom_shadowtext(aes(label= sprintf("%0.2f\n(%0.2f)",pooled.abs.mean.lol.coeff,pooled.abs.sd.lol.coeff)),color="white", size = 4.5)+
  scale_fill_viridis(discrete=FALSE) +
  guides(fill = guide_colourbar(title = sprintf('Mean absolute LOL coefficients'),
                                barwidth = 40,
                                direction="horizontal",
                                title.position = 'top',
                                frame.colour=c("black"),
                                frame.linewidth = 1.5,
                                title.hjust = .5)) +
  xlab(label = "Sensor") +
  ylab(label = "Feature") +
  scale_y_discrete(limits = rev(c("SMA","HLF (h)","HLF (l)","MFR","FDE","BPW","WVL","All")))+
  scale_x_discrete(limits = c("RE","LE","RW","LW","RA","LA","All"))+
  facet_wrap(~observation.window, nrow = 2, ncol = 2) +
  theme_classic()+
  theme(
    strip.text = element_text(size=22), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size = 2),
    legend.position = 'bottom',
    legend.title = element_text(size = 18, color = "black"),
    legend.text=element_text(size=16),
    aspect.ratio = 1
  )

### VI. Figure 6: ROC and PR curves of the best-performing models
## Load filtered ROC axes
filtered.group.curve.axes <- read.csv('../results/detection_results/regression_results/tuned_model_metrics/plot_averaged_curve_axes.csv')

## Load filtered PRC axes
filtered.group.prc.axes <- read.csv('../results/detection_results/regression_results/tuned_model_metrics/plot_averaged_prc_axes.csv')

## ggplot: ROC curves
roc.plots <- filtered.group.curve.axes %>%
  mutate(model.name = str_sub(model,1,4),
         GCSm = paste('GCSm:',GCSm)) %>%
  ggplot(aes(x = fpr)) +
  facet_wrap( ~ GCSm,
              ncol = 2,
              scales = 'free') +
  xlab("False Positive Rate") +
  ylab("True Positive Rate") +
  coord_cartesian(ylim = c(0,1),xlim = c(0,1))+
  geom_ribbon(aes(ymin = lowerTPR, ymax = upperTPR, fill = factor(obs.window)), alpha = 0.3) +
  geom_line(aes(y = meanTPR, color = factor(obs.window),linetype = factor(model.name)), alpha = 0.75, size=1.20) +
  guides(linetype = FALSE, color = guide_legend(nrow = 2)) +
  geom_abline(
    intercept = 0,
    slope = 1,
    alpha = 0.5,
    linetype = "dashed",
    size=1
  ) +
  theme_classic()+
  theme(
    strip.text = element_text(size=22), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size = 2),
    aspect.ratio = 1
  )

## ggplot: PR curves
load('../validation_resampling/detection_labels.RData')
get.prop <- function(x){
  table(x$Best.Motor.Response)/sum(table(x$Best.Motor.Response))
}
props.of.classes <- as.data.frame(matrix(ncol = 2,nrow = 0))
for(i in 1:length(det_gcs_labels)){
  props.of.classes <- rbind(props.of.classes,as.data.frame(get.prop(det_gcs_labels[[i]])))
}
grouped.freqs <- props.of.classes %>%
  group_by(Var1) %>%
  summarise(freq = min(Freq)) %>%
  rename(GCSm = Var1)
filtered.group.prc.axes <- left_join(filtered.group.prc.axes, grouped.freqs, by = 'GCSm')
prc.plots <- filtered.group.prc.axes %>%
  mutate(model.name = as.factor(toupper(str_sub(model,1,4))),
         GCSm = paste('GCSm:',GCSm)) %>%
  ggplot(aes(x = rec)) +
  facet_wrap( ~ GCSm,
              ncol = 2,
              scales = 'free') +
  xlab("Recall") +
  ylab("Precision") +
  coord_cartesian(ylim = c(0,1),xlim = c(0,1))+
  geom_ribbon(aes(ymin = lowerPREC, ymax = upperPREC, fill = factor(obs.window)), alpha = 0.3) +
  geom_line(aes(y = meanPREC, color = factor(obs.window),linetype=model.name), alpha = 0.75, size=1.20) +
  guides(linetype = guide_legend(title = 'Model Type',nrow = 1, order = 2), 
         color = guide_legend(title = 'Observation Window (h)',nrow = 1,order = 1),
         fill = FALSE) +
  geom_hline(
    aes(yintercept = freq),
    alpha = 0.5,
    linetype = "dashed",
    size=1
  ) +
  theme_classic()+
  theme(
    strip.text = element_text(size=22), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size = 2),
    legend.position = "bottom",
    legend.title = element_text(size=22),
    legend.text=element_text(size=16),
    aspect.ratio = 1
  )

### VII. Figure 7: Mean predicted probability by observation window

## Load partition and resampling information
load('../validation_resampling/detection_labels.RData')
load('../validation_resampling/detection_partitions.RData')

## Load optimal (by AUROC) model information and group by observation window, model type, and GCSm
group.aucs <- read.csv('../results/detection_results/regression_results/tuned_model_metrics/averaged_aucs.csv') %>% 
  mutate(model.id = paste0(model,'_',obs.window),
         model.type = str_sub(model,1,4))
opt.r.df <- read.csv('../results/detection_results/regression_results/optimal_tuning_results.csv') %>% 
  mutate(model.id = paste0(model,'_',obs.window)) %>%
  select(model.id,maxAUROCdim)
opt.group.aucs <- group.aucs %>%
  group_by(obs.window, model.type, GCSm) %>%
  summarise(
    model = model[which.max(meanAUROC)],
    model.type = model.type[which.max(meanAUROC)],
    maxAUROC = max(meanAUROC),
    model.id = model.id[which.max(meanAUROC)],
    sdAUROC = sdAUROC[which.max(meanAUROC)]
  ) %>% 
  left_join(opt.r.df,by = "model.id")

## Load imputation directories
imp.dirs <- list.files('../results/detection_results/regression_results',pattern = '^imp*',include.dirs = TRUE,full.names = TRUE)

## Cycle through optimal models for plot (and through imputations) to source and store mean predicted probabilities
stored.prob.values <- as.data.frame(matrix(ncol = 6, nrow = 0))
for (i in 1:nrow(opt.group.aucs)){
  print(paste("Model no.",i,"out of",nrow(opt.group.aucs),"started."))

  if (opt.group.aucs$model[i] == 'mnlr'){
    model.file.name <- paste0('mnlr_mdl_',opt.group.aucs$maxAUROCdim[i],'_preds.csv')
  } else if (opt.group.aucs$model[i] == 'mnlr.smote') {
    model.file.name <- paste0('SMOTE_mnlr_mdl_',opt.group.aucs$maxAUROCdim[i],'_preds.csv')
  } else if (opt.group.aucs$model[i] == 'polr') {
    model.file.name <- paste0('polr_mdl_',opt.group.aucs$maxAUROCdim[i],'_preds.csv')
  } else if (opt.group.aucs$model[i] == 'polr.smote') {
    model.file.name <- paste0('SMOTE_polr_mdl_',opt.group.aucs$maxAUROCdim[i],'_preds.csv')
  }
  curr.mdl.prob.values <- as.data.frame(matrix(ncol = 6, nrow = 0))
  for (j in 1:length(imp.dirs)){
    print(paste("Imputation no.",j,"out of",length(imp.dirs),"started."))
    curr.mdl.results <-
      read.csv(file.path(
        '../results/detection_results/regression_results',
        paste0('imp', j),
        paste0('detection_window_',opt.group.aucs$obs.window[i]),
        model.file.name
      )) %>%
      mutate(true.labels = factor(true.labels, levels = 1:6),
             pred.labels = factor(pred.labels, levels = 1:6))
    curr.label <- opt.group.aucs$GCSm[i]
    temp.results.df <- curr.mdl.results %>% 
      mutate(temp.label = as.numeric(true.labels == curr.label))
    curr.label.prob.name <- names(temp.results.df)[grep(paste('GCSm',curr.label,sep = "."),names(temp.results.df))]
    
    curr.mdl.prob.values <- rbind(curr.mdl.prob.values,
                                  data.frame(model.type = opt.group.aucs$model.type[i],
                                             GCSm = opt.group.aucs$GCSm[i],
                                             obs.window = opt.group.aucs$obs.window[i],
                                             prob = temp.results.df[[curr.label.prob.name]],
                                             pos.indicator = temp.results.df$temp.label,
                                             imp.no = j,
                                             test.idx = 1:nrow(temp.results.df)
                                  ))
  }
  grouped.curr.mdl.prob.values <- curr.mdl.prob.values %>%
    group_by(model.type, GCSm, obs.window, pos.indicator,test.idx) %>%
    summarise(mean.prob = mean(prob,na.rm = TRUE), sd.prob = sd(prob,na.rm = TRUE)) %>%
    dplyr::select(-test.idx)
  stored.prob.values <- rbind(stored.prob.values, grouped.curr.mdl.prob.values)
}

## Average the stored probability values across imputations
grouped.stored.prob.values <- stored.prob.values %>%
  group_by(model.type,GCSm,obs.window,pos.indicator) %>%
  summarise(mean.macro.prob = mean(mean.prob,na.rm = TRUE),
            sd.macro.prob = sqrt((1/n())*sum(sd.prob^2)+((n()+1)/(n()))*(sd(mean.prob)^2) )) %>%
  rowwise() %>%
  mutate(upper.bound = min(mean.macro.prob + sd.macro.prob,1),
         lower.bound = max(mean.macro.prob - sd.macro.prob,0))

## ggplot: MNLR probability means by observation window
mnlr.ow.plot <- grouped.stored.prob.values %>%
  filter(model.type == 'mnlr') %>%
  mutate(GCSm = paste('GCSm:',GCSm)) %>%
  ggplot(aes(x = as.factor(obs.window), color = as.factor(pos.indicator), y = mean.macro.prob, group = as.factor(pos.indicator))) +
  geom_line(position=position_dodge(0.2)) + 
  geom_pointrange(aes(ymin = lower.bound, ymax = upper.bound),position=position_dodge(0.2)) +
  facet_wrap(~GCSm,
             ncol = 2,
             scales = 'free') +
  ylab("Mean Predicted Probability") +
  xlab("Observation Window (h)") +
  coord_cartesian(ylim = c(0,1)) +
  scale_color_manual(values = c('red','blue'))+
  theme_classic()+
  theme(
    strip.text = element_text(size=22), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size = 2),
    legend.position = "none",
    legend.title = element_text(size=22),
    legend.text=element_text(size=16),
    plot.title = element_text(hjust = 0.5),
    aspect.ratio = 1
  )

## ggplot: POLR probability means by observation window
polr.ow.plot <- grouped.stored.prob.values %>%
  filter(model.type == 'polr') %>%
  mutate(GCSm = paste('GCSm:',GCSm)) %>%
  ggplot(aes(x = as.factor(obs.window), color = as.factor(pos.indicator), y = mean.macro.prob, group = as.factor(pos.indicator))) +
  geom_line(position=position_dodge(0.2)) + 
  geom_pointrange(aes(ymin = lower.bound, ymax = upper.bound),position=position_dodge(0.2)) +
  facet_wrap(~GCSm,
             ncol = 2,
             scales = 'free') +
  ylab("Mean Predicted Probability") +
  xlab("Observation Window (h)") +
  coord_cartesian(ylim = c(0,1)) +
  scale_color_manual(values = c('red','blue'))+
  theme_classic()+
  theme(
    strip.text = element_text(size=22), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size = 2),
    legend.position = "none",
    legend.title = element_text(size=22),
    legend.text=element_text(size=16),
    plot.title = element_text(hjust = 0.5),
    aspect.ratio = 1
  )

### VIII. Supplementary Figure 2: Plot complete GCS trajectory of patients in the study
## Load automatically extracted GCS labels with indicator labels:
gcs.data <- read.csv('../clinical_data/02_clean_auto_GCS_table.csv') %>% 
  select(-X) %>% 
  mutate(TakenInstant = as.POSIXct(TakenInstant, tz = "America/New_York"))

## Load missing time info (includes information on recording start and stop times)
missing.time.info <- read_xlsx('../all_motion_feature_data/MissingPercentTable.xlsx',.name_repair = "universal") %>% 
  mutate(Start.Timestamp = as.POSIXct(Start.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"), End.Timestamp = as.POSIXct(End.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York")) %>% 
  mutate(Recording.Duration = End.Timestamp - Start.Timestamp) %>%
  rename(AccelPatientNo_=Accel.Patient.No.) %>%
  mutate(Bed.miss.hours = Bed*Recording.Duration,LA.miss.hours = LA*Recording.Duration,LE.miss.hours = LE*Recording.Duration,LW.miss.hours = LW*Recording.Duration,RA.miss.hours = RA*Recording.Duration,RE.miss.hours = RE*Recording.Duration,RW.miss.hours = RW*Recording.Duration)

## Merge GCS and missing data info
unique.accel.nums <- sort(unique(gcs.data$AccelPatientNo_))
patient.ids <- paste('Patient No.',1:length(unique.accel.nums))
gcs.data <- left_join(gcs.data,data.frame(AccelPatientNo_ = unique.accel.nums, patient.ids), by = 'AccelPatientNo_')
missing.time.info <- left_join(missing.time.info,data.frame(AccelPatientNo_ = unique.accel.nums, patient.ids), by = 'AccelPatientNo_')

## Convert time-stamps to days in ICU stay
time.starting.info <- gcs.data %>%
  arrange(TakenInstant) %>%
  group_by(AccelPatientNo_, patient.ids) %>%
  summarise(first.timestamp = first(TakenInstant))
gcs.data$days <- NA
missing.time.info$start.days <- NA
missing.time.info$end.days <- NA
for (i in 1:nrow(time.starting.info)){
  curr.start.time <- time.starting.info$first.timestamp[i]
  curr.accel.no <- time.starting.info$AccelPatientNo_[i]
  curr.idx <- which(gcs.data$AccelPatientNo_ == curr.accel.no)
  gcs.data$days[curr.idx] <- as.numeric(difftime(gcs.data$TakenInstant[curr.idx], curr.start.time, units = "days"))
  curr.idx <- which(missing.time.info$AccelPatientNo_ == curr.accel.no)
  missing.time.info$start.days[curr.idx] <- as.numeric(difftime(missing.time.info$Start.Timestamp[curr.idx], curr.start.time, units = "days"))
  missing.time.info$end.days[curr.idx] <- as.numeric(difftime(missing.time.info$End.Timestamp[curr.idx], curr.start.time, units = "days"))
}
gcs.data <- left_join(gcs.data,
                      missing.time.info[,c('AccelPatientNo_','start.days','end.days')],
                      by = 'AccelPatientNo_')

## Define partitions for different plots
partition1 <- patient.ids[1:14]
partition2 <- patient.ids[1:14 + 14]
partition3 <- patient.ids[1:14 + 2*14]
partition4 <- patient.ids[1:14 + 3*14]
partition5 <- patient.ids[1:13 + 4*14]

## ggplot: produce GCSm trajectory and concurrent recording indicator plots for partitions of patients
ordering.ids <- patient.ids
gcs.plot.1 <- gcs.data %>%
  mutate(patient.ids = factor(patient.ids,levels = ordering.ids)) %>%
  filter(patient.ids %in% partition1) %>%
  drop_na(Best.Motor.Response) %>%
  ggplot(aes(x = days, y = Best.Motor.Response)) +
  geom_rect(inherit.aes = FALSE, 
            data = missing.time.info %>%
              mutate(patient.ids = factor(patient.ids,levels = ordering.ids)) %>%
              filter(patient.ids %in% partition1), 
            aes(xmin=start.days,xmax=end.days,ymin=0, ymax=7),
            fill='red',
            alpha = .2) +
  geom_line(size = .75,) +
  geom_point(size = .85) +
  ylab('GCSm')+
  scale_x_continuous(breaks= pretty_breaks(),expand = c(0,0)) +
  coord_cartesian(ylim=c(1, 6)) +
  facet_wrap(~patient.ids,
             scales = 'free_x',
             ncol = 2) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size=14, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

gcs.plot.2 <- gcs.data %>%
  mutate(patient.ids = factor(patient.ids,levels = ordering.ids)) %>%
  filter(patient.ids %in% partition2) %>%
  drop_na(Best.Motor.Response) %>%
  ggplot(aes(x = days, y = Best.Motor.Response)) +
  geom_rect(inherit.aes = FALSE, 
            data = missing.time.info %>%
              mutate(patient.ids = factor(patient.ids,levels = ordering.ids)) %>%
              filter(patient.ids %in% partition2), 
            aes(xmin=start.days,xmax=end.days,ymin=0, ymax=7),
            fill='red',
            alpha = .2) +
  geom_line(size = .75,) +
  geom_point(size = .85) +
  ylab('GCSm')+
  scale_x_continuous(breaks= pretty_breaks(),expand = c(0,0)) +
  coord_cartesian(ylim=c(1, 6)) +
  facet_wrap(~patient.ids,
             scales = 'free_x',
             ncol = 2) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size=14, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

gcs.plot.3 <- gcs.data %>%
  mutate(patient.ids = factor(patient.ids,levels = ordering.ids)) %>%
  filter(patient.ids %in% partition3) %>%
  drop_na(Best.Motor.Response) %>%
  ggplot(aes(x = days, y = Best.Motor.Response)) +
  geom_rect(inherit.aes = FALSE, 
            data = missing.time.info %>%
              mutate(patient.ids = factor(patient.ids,levels = ordering.ids)) %>%
              filter(patient.ids %in% partition3), 
            aes(xmin=start.days,xmax=end.days,ymin=0, ymax=7),
            fill='red',
            alpha = .2) +
  geom_line(size = .75,) +
  geom_point(size = .85) +
  ylab('GCSm')+
  scale_x_continuous(breaks= pretty_breaks(),expand = c(0,0)) +
  coord_cartesian(ylim=c(1, 6)) +
  facet_wrap(~patient.ids,
             scales = 'free_x',
             ncol = 2) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size=14, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

gcs.plot.4 <- gcs.data %>%
  mutate(patient.ids = factor(patient.ids,levels = ordering.ids)) %>%
  filter(patient.ids %in% partition4) %>%
  drop_na(Best.Motor.Response) %>%
  ggplot(aes(x = days, y = Best.Motor.Response)) +
  geom_rect(inherit.aes = FALSE, 
            data = missing.time.info %>%
              mutate(patient.ids = factor(patient.ids,levels = ordering.ids)) %>%
              filter(patient.ids %in% partition4), 
            aes(xmin=start.days,xmax=end.days,ymin=0, ymax=7),
            fill='red',
            alpha = .2) +
  geom_line(size = .75,) +
  geom_point(size = .85) +
  ylab('GCSm')+
  scale_x_continuous(breaks= pretty_breaks(),expand = c(0,0)) +
  coord_cartesian(ylim=c(1, 6)) +
  facet_wrap(~patient.ids,
             scales = 'free_x',
             ncol = 2) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size=14, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

gcs.plot.5 <- gcs.data %>%
  mutate(patient.ids = factor(patient.ids,levels = ordering.ids)) %>%
  filter(patient.ids %in% partition5) %>%
  drop_na(Best.Motor.Response) %>%
  ggplot(aes(x = days, y = Best.Motor.Response)) +
  geom_rect(inherit.aes = FALSE, 
            data = missing.time.info %>%
              mutate(patient.ids = factor(patient.ids,levels = ordering.ids)) %>%
              filter(patient.ids %in% partition5), 
            aes(xmin=start.days,xmax=end.days,ymin=0, ymax=7),
            fill='red',
            alpha = .2) +
  geom_line(size = .75,) +
  geom_point(size = .85) +
  ylab('GCSm')+
  scale_x_continuous(breaks= pretty_breaks(),expand = c(0,0)) +
  coord_cartesian(ylim=c(1, 6)) +
  facet_wrap(~patient.ids,
             scales = 'free_x',
             ncol = 2) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    strip.text = element_text(size=14, color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

### IX. Supplementary Figure 3: Plot histograms characterizing recording duration and start time of recording
## Load automatically extracted GCS labels
gcs.data <- read.csv('../clinical_data/02_clean_auto_GCS_table.csv') %>% 
  select(-X) %>% 
  mutate(TakenInstant = as.POSIXct(TakenInstant, tz = "America/New_York"))

## Load missing time info (includes information on recording start and stop times)
missing.time.info <- read_xlsx('../all_motion_feature_data/MissingPercentTable.xlsx',.name_repair = "universal") %>% 
  mutate(Start.Timestamp = as.POSIXct(Start.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"), End.Timestamp = as.POSIXct(End.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York")) %>% 
  mutate(Recording.Duration = End.Timestamp - Start.Timestamp) %>%
  rename(AccelPatientNo_=Accel.Patient.No.) %>%
  mutate(Bed.miss.hours = Bed*Recording.Duration,LA.miss.hours = LA*Recording.Duration,LE.miss.hours = LE*Recording.Duration,LW.miss.hours = LW*Recording.Duration,RA.miss.hours = RA*Recording.Duration,RE.miss.hours = RE*Recording.Duration,RW.miss.hours = RW*Recording.Duration)

## Combine GCS and missing time info by matching patient indices
unique.accel.nums <- sort(unique(gcs.data$AccelPatientNo_))
patient.ids <- paste('Patient No.',1:length(unique.accel.nums))
gcs.data <- left_join(gcs.data,data.frame(AccelPatientNo_ = unique.accel.nums, patient.ids), by = 'AccelPatientNo_')
missing.time.info <- left_join(missing.time.info,data.frame(AccelPatientNo_ = unique.accel.nums, patient.ids), by = 'AccelPatientNo_')

## Convert time stamps to days in ICU stay
time.starting.info <- gcs.data %>%
  arrange(TakenInstant) %>%
  group_by(AccelPatientNo_, patient.ids) %>%
  summarise(first.timestamp = first(TakenInstant))
gcs.data$days <- NA
missing.time.info$start.days <- NA
missing.time.info$end.days <- NA
for (i in 1:nrow(time.starting.info)){
  curr.start.time <- time.starting.info$first.timestamp[i]
  curr.accel.no <- time.starting.info$AccelPatientNo_[i]
  curr.idx <- which(gcs.data$AccelPatientNo_ == curr.accel.no)
  gcs.data$days[curr.idx] <- as.numeric(difftime(gcs.data$TakenInstant[curr.idx], curr.start.time, units = "days"))
  curr.idx <- which(missing.time.info$AccelPatientNo_ == curr.accel.no)
  missing.time.info$start.days[curr.idx] <- as.numeric(difftime(missing.time.info$Start.Timestamp[curr.idx], curr.start.time, units = "days"))
  missing.time.info$end.days[curr.idx] <- as.numeric(difftime(missing.time.info$End.Timestamp[curr.idx], curr.start.time, units = "days"))
}

## Arrange missing info dataframe by patient index
missing.time.info <- missing.time.info %>%
  arrange(AccelPatientNo_)

## Load patient clinical data to source recording start date 
source('./functions/load_patient_clinical_data.R')
patient.clinical.data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv') %>% 
  arrange(AccelPatientNo_) %>% 
  mutate(ptIdx = 1:nrow(.))

## Format variables into neat dataframe for ggplot
histogram.df <- data.frame(value = as.numeric(missing.time.info$Recording.Duration), 
                           type = 'Recording Duration (h)')
histogram.df <- rbind(histogram.df,
                      data.frame(value = as.numeric(patient.clinical.data$AccelRecordingStartDate - patient.clinical.data$NCCUAdmissionDate), 
                                 type = 'NCCU Admission to Recording (d)'),
                      data.frame(value = (missing.time.info$end.days-missing.time.info$start.days)/patient.clinical.data$DaysInNCCU, 
                                 type = 'Proportion of NCCU Stay Recorded'),
                      data.frame(value = as.numeric(patient.clinical.data$AccelRecordingStartDate - patient.clinical.data$NCCUAdmissionDate)/as.numeric(patient.clinical.data$NCCUDischargeDate - patient.clinical.data$NCCUAdmissionDate), 
                                 type = 'Proportion of NCCU Stay before Recording'))

## ggplot: produce characteristic histograms
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

### X. Supplementary Figure 5: Plot hyperparameter optimization tuning results
## Load compiled tuning results
compiled.tuning.df <- read.csv('../results/detection_results/regression_results/compiled_tuning_results.csv')

## ggplot: produce stratified line plots of tuning results
tuning.plot <- compiled.tuning.df %>%
  group_by(model,r,obs.window) %>%
  summarise(meanAUROCS = mean(aurocs),sdAUROCS = sd(aurocs),n = n()) %>%
  mutate(obs.window = factor(obs.window)) %>%
  mutate(obs.window = plyr::mapvalues(obs.window,from = c('0.5','1','3','6'), 
                                      to = c("0.5 hr window", "1 hr window","3 hr window","6 hr window")),
         Model = plyr::mapvalues(model,from = c("polr","polr.smote","mnlr","mnlr.smote"), 
                                 to = c("POLR","POLR.SMOTE","MNLR","MNLR.SMOTE"))) %>%
  ggplot(aes(x = r, y = meanAUROCS, color = Model)) +
  facet_wrap(~obs.window) +
  geom_line(size = .75) +
  geom_point(size = .75) + 
  scale_x_continuous(breaks = r,expand = c(0,0)) +
  geom_errorbar(aes(ymin = meanAUROCS - sdAUROCS,max = meanAUROCS + sdAUROCS),size = .75) +
  xlab('LOL Target Dimensionality (d)')+
  ylab('Macro-averaged AUROC')+
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 22, color = "black"),
    axis.title.y = element_text(size = 22, color = "black"),
    strip.text = element_text(size=22, color = "black"), 
    strip.placement = 'outside',
    strip.background = element_blank(),
    legend.position = 'bottom',
    legend.title = element_text(size = 18, color = "black"),
    legend.text=element_text(size=18),
    aspect.ratio = .75)

### XI. Supplementary Figure 6: Plot probability calibration curves
## Load calibration curve information from python
calib.curve.df <- read.csv('../results/detection_results/regression_results/compiled_calibration_curves.csv') %>%
  mutate(model.id = paste0(model,'_',obs_window))

## Loop through calibration curve results and interpolate onto a common axis
xq <- seq(from = 0, to = 1, length.out = 7)
calib.interpol <- data.frame(matrix(ncol = 5, nrow = 0))
uniq.model.ids <- unique(calib.curve.df$model.id)
for (i in 1:length(uniq.model.ids)){
  print(paste("Model no.",i,"out of",length(uniq.model.ids),"started."))
  uniq.imps <- unique(calib.curve.df$imp[calib.curve.df$model.id == uniq.model.ids[i]])
  for (j in 1:length(uniq.imps)){
    print(paste("Imputation no.",j,"out of",length(uniq.imps),"started."))
    unique.labels <- sort(unique(calib.curve.df$GCSm[calib.curve.df$model.id == uniq.model.ids[i] & calib.curve.df$imp == uniq.imps[j]]))
    for (k in 1:length(unique.labels)){
      print(paste("GCSm",unique.labels[k],"started."))
      curr.label <- unique.labels[k]
      curr.idx <- which(calib.curve.df$GCSm == curr.label & calib.curve.df$imp == uniq.imps[j] & calib.curve.df$model.id == uniq.model.ids[i])
      if (length(curr.idx) == 0){
        next
      } else if (length(curr.idx) == 1){
        calib.interpol <- rbind(calib.interpol, data.frame(prob_pred = xq[which.min(abs(calib.curve.df$prob_pred[curr.idx] - xq))], prob_true = calib.curve.df$prob_true[curr.idx], model.id = uniq.model.ids[i], imp =  uniq.imps[j], GCSm = curr.label))
        next
      }
      interpol.object <- approx(x = calib.curve.df$prob_pred[curr.idx],y = calib.curve.df$prob_true[curr.idx],xout = xq)
      calib.interpol <- rbind(calib.interpol, data.frame(prob_pred = interpol.object$x, prob_true = interpol.object$y, model.id = uniq.model.ids[i], imp = uniq.imps[j], GCSm = curr.label))
    }
  }
}

## Save interpolated calibration curve dataframe
write.csv(calib.interpol,'../results/detection_results/regression_results/interpolated_calibration_curves.csv',row.names = FALSE)

## Calculate mean and std of relevant calibration values in interpolated dataframe
calib.interpol.summ <- calib.interpol %>%
  drop_na(prob_true) %>%
  group_by(model.id,GCSm,prob_pred) %>%
  summarise(meanPROB_TRUE = mean(prob_true), upperPROB_TRUE = min(mean(prob_true)+sd(prob_true),1),lowerPROB_TRUE = max(mean(prob_true)-sd(prob_true),0))

## Filter out AUROC calibration curves (A)
auroc.calib.interpol <- as.data.frame(matrix(ncol=6,nrow=0)) 
for (i in 1:nrow(opt.group.aucs)){
  curr.model.id <- opt.group.aucs$model.id[i]
  curr.GCSm <- opt.group.aucs$GCSm[i]
  curr.frame <- filter(calib.interpol.summ, model.id == curr.model.id, GCSm == curr.GCSm)
  auroc.calib.interpol <- rbind(auroc.calib.interpol,curr.frame)
}
auroc.calib.interpol$obs.window <- sub(".*_", "", auroc.calib.interpol$model.id) 
auroc.calib.interpol$model.name <- str_sub(auroc.calib.interpol$model.id,1,4)

## ggplot: produce calibration curves of AUROC models
auroc.calib.plot <- auroc.calib.interpol %>%
  mutate(model.name = as.factor(toupper(model.name)),
         GCSm = paste('GCSm:',GCSm)) %>%
  ggplot(aes(x = prob_pred)) +
  facet_wrap(~ GCSm,
             ncol = 2,
             scales = 'free') +
  xlab("Mean Predicted Probability") +
  ylab("Fraction of Positives") +
  coord_cartesian(ylim = c(0,1),xlim = c(0,1)) +
  geom_abline(
    intercept = 0,
    slope = 1,
    alpha = 0.5,
    linetype = "dashed",
    size=.75
  ) +
  geom_line(aes(y = meanPROB_TRUE, color = factor(obs.window),linetype=model.name), alpha = 0.75, size=1.20) +
  geom_pointrange(aes(y = meanPROB_TRUE, ymin = lowerPROB_TRUE, ymax = upperPROB_TRUE, color = factor(obs.window)),position=position_dodge(0.05)) +
  theme_classic()+
  theme(
    strip.text = element_text(size=22), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size = 2),
    legend.position = "none",
    legend.title = element_text(size=22),
    legend.text=element_text(size=16),
    aspect.ratio = 1
  )

## Filter out AUPRC calibration curves (B)
auprc.calib.interpol <- as.data.frame(matrix(ncol=6,nrow=0)) 
for (i in 1:nrow(opt.group.aucprs)){
  curr.model.id <- opt.group.aucprs$model.id[i]
  curr.GCSm <- opt.group.aucprs$GCSm[i]
  curr.frame <- filter(calib.interpol.summ, model.id == curr.model.id, GCSm == curr.GCSm)
  auprc.calib.interpol <- rbind(auprc.calib.interpol,curr.frame)
}
auprc.calib.interpol$obs.window <- sub(".*_", "", auprc.calib.interpol$model.id) 
auprc.calib.interpol$model.name <- str_sub(auprc.calib.interpol$model.id,1,4)

## ggplot: produce calibration curves of AUPRC models
auprc.calib.plot <- auprc.calib.interpol %>%
  mutate(model.name = as.factor(toupper(model.name)),
         GCSm = paste('GCSm:',GCSm)) %>%
  ggplot(aes(x = prob_pred)) +
  facet_wrap(~ GCSm,
             ncol = 2,
             scales = 'free') +
  xlab("Mean Predicted Probability") +
  ylab("Fraction of Positives") +
  coord_cartesian(ylim = c(0,1),xlim = c(0,1)) +
  geom_abline(
    intercept = 0,
    slope = 1,
    alpha = 0.5,
    linetype = "dashed",
    size=.75
  ) +
  geom_line(aes(y = meanPROB_TRUE, color = factor(obs.window),linetype=model.name), alpha = 0.75, size=1.20) +
  geom_pointrange(aes(y = meanPROB_TRUE, ymin = lowerPROB_TRUE, ymax = upperPROB_TRUE, color = factor(obs.window)),position=position_dodge(0.05)) +
  theme_classic()+
  theme(
    strip.text = element_text(size=22), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size = 2),
    legend.position = "none",
    legend.title = element_text(size=22),
    legend.text=element_text(size=16),
    aspect.ratio = 1
  )
