# Load clinical patient data:
source('./functions/load_patient_clinical_data.R')
patient_clinical_data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv') %>% arrange(AccelPatientNo_) %>% mutate(ptIdx = 1:nrow(.))

# Load Motion Features (all)
if (!exists("all_motion_features")) {
  all_sensors <-
    readMat('../all_motion_feature_data/complete_sensor_data.mat')$sensors
  all_motion_features <- do.call(rbind, all_sensors)
  featureLabels <- read.csv('../all_motion_feature_data/feature_names.csv',header = FALSE)
  all_mf_times <- read.csv('../all_motion_feature_data/indexed_times.csv') %>% mutate(times = as.POSIXct(times,tz="America/New_York",format ="%d-%b-%Y %H:%M:%S"))
}

n <- length(all_motion_features)/length(featureLabels)

source('./functions/mf_to_dataframe.R')

# Recode all missing values to NA, find "totally missing" data-streams, and store all feature values in single DF:
out.MF2DF <- mf_to_dataframe(all_motion_features,n,verbose = TRUE) 
completeFeatureSet <- out.MF2DF[[1]]
totallyMissingSet <- out.MF2DF[[2]]
featureLabels <- unlist(featureLabels[1,])
complete_timestamps <- as.data.frame(matrix(ncol = 3, nrow = 0))
for (i in 1:n){
  curr_timestamps <- all_mf_times %>% filter(ptIdx == i)
  curr_timestamps$timeCount <- seq(nrow(curr_timestamps))
  complete_timestamps <- rbind(complete_timestamps,curr_timestamps)
}
complete_timestamps <- rename(complete_timestamps,timeStamps = times)
completeFeatureSet <- inner_join(completeFeatureSet,complete_timestamps,by=c("ptIdx","timeCount"))

myTMS <- totallyMissingSet %>% group_by(ptIdx, srIdx) %>% summarise()

# Load missing accelerometry information:
missing.time.info <- read_xlsx('../all_motion_feature_data/MissingPercentTable.xlsx',.name_repair = "universal") %>% 
  mutate(Start.Timestamp = as.POSIXct(Start.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"), End.Timestamp = as.POSIXct(End.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York")) %>% 
  mutate(Recording.Duration = End.Timestamp - Start.Timestamp) %>%
  rename(AccelPatientNo_=Accel.Patient.No.) %>%
  mutate(Bed.miss.hours = Bed*Recording.Duration,LA.miss.hours = LA*Recording.Duration,LE.miss.hours = LE*Recording.Duration,LW.miss.hours = LW*Recording.Duration,RA.miss.hours = RA*Recording.Duration,RE.miss.hours = RE*Recording.Duration,RW.miss.hours = RW*Recording.Duration)

# Load no motion accelerometry information:
noMotionTimeInfo <- read_xlsx('../all_motion_feature_data/NoMotionPercentTable.xlsx',.name_repair = "universal") %>% 
  mutate(Recording.Duration = missing.time.info$End.Timestamp - missing.time.info$Start.Timestamp) %>%
  rename(AccelPatientNo_=Accel.Patient.No.)
# Correct Totally Missing Sets in No Motion Info:
for (i in 1:nrow(myTMS)){
  currPtIdx <- myTMS$ptIdx[i]
  currSrIdx <- myTMS$srIdx[i]
  noMotionTimeInfo[currPtIdx,currSrIdx+1] <- NA
}
noMotionTimeInfo <- noMotionTimeInfo %>%
  mutate(Bed.nm.hours = Bed*Recording.Duration,LA.nm.hours = LA*Recording.Duration,LE.nm.hours = LE*Recording.Duration,LW.nm.hours = LW*Recording.Duration,RA.nm.hours = RA*Recording.Duration,RE.nm.hours = RE*Recording.Duration,RW.nm.hours = RW*Recording.Duration)

# LA percent nomotion
LA.perc.nomotion <- noMotionTimeInfo %>% drop_na(LA)
perc.LA.motion <- 100*as.numeric(sum(LA.perc.nomotion$LA.nm.hours))/as.numeric(sum(LA.perc.nomotion$Recording.Duration))

# LE percent nomotion
LE.perc.nomotion <- noMotionTimeInfo %>% drop_na(LE)
perc.LE.motion <- 100*as.numeric(sum(LE.perc.nomotion$LE.nm.hours))/as.numeric(sum(LE.perc.nomotion$Recording.Duration))

# LW percent nomotion
LW.perc.nomotion <- noMotionTimeInfo %>% drop_na(LW)
perc.LW.motion <- 100*as.numeric(sum(LW.perc.nomotion$LW.nm.hours))/as.numeric(sum(LW.perc.nomotion$Recording.Duration))

# RA percent nomotion
RA.perc.nomotion <- noMotionTimeInfo %>% drop_na(RA)
perc.RA.motion <- 100*as.numeric(sum(RA.perc.nomotion$RA.nm.hours))/as.numeric(sum(RA.perc.nomotion$Recording.Duration))

# RE percent nomotion
RE.perc.nomotion <- noMotionTimeInfo %>% drop_na(RE)
perc.RE.motion <- 100*as.numeric(sum(RE.perc.nomotion$RE.nm.hours))/as.numeric(sum(RE.perc.nomotion$Recording.Duration))

# RW percent nomotion
RW.perc.nomotion <- noMotionTimeInfo %>% drop_na(RW)
perc.RW.motion <- 100*as.numeric(sum(RW.perc.nomotion$RW.nm.hours))/as.numeric(sum(RW.perc.nomotion$Recording.Duration))

# Bed percent nomotion
Bed.perc.nomotion <- noMotionTimeInfo %>% drop_na(Bed) %>% filter(AccelPatientNo_ != 64)
perc.Bed.motion <- 100*as.numeric(sum(Bed.perc.nomotion$Bed.nm.hours))/as.numeric(sum(Bed.perc.nomotion$Recording.Duration))

# Total percent nomotion
Total.hours.motion <- as.numeric(sum(LA.perc.nomotion$LA.nm.hours)) + as.numeric(sum(LE.perc.nomotion$LE.nm.hours)) + as.numeric(sum(LW.perc.nomotion$LW.nm.hours)) + as.numeric(sum(RA.perc.nomotion$RA.nm.hours)) + as.numeric(sum(RE.perc.nomotion$RE.nm.hours)) + as.numeric(sum(RW.perc.nomotion$RW.nm.hours))
Total.hours.of.set <- as.numeric(sum(LA.perc.nomotion$Recording.Duration)) + as.numeric(sum(LE.perc.nomotion$Recording.Duration)) + as.numeric(sum(LW.perc.nomotion$Recording.Duration)) + as.numeric(sum(RA.perc.nomotion$Recording.Duration)) + as.numeric(sum(RE.perc.nomotion$Recording.Duration)) + as.numeric(sum(RW.perc.nomotion$Recording.Duration))
perc.total.motion <-  100*Total.hours.motion/Total.hours.of.set

Total.noMotionTimeInfo <- noMotionTimeInfo %>% mutate(num.col = (rowSums(!is.na(.))-2)/2) %>% mutate(total.hours = Recording.Duration*num.col)
Total.noMotionTimeInfo$Total.Hours.of.Motion <- rowSums(sapply(noMotionTimeInfo[,10:16], as.numeric),na.rm = TRUE)
Total.noMotionTimeInfo$Total.Percs.of.Motion <- 100*Total.noMotionTimeInfo$Total.Hours.of.Motion/as.numeric(Total.noMotionTimeInfo$total.hours)

complete.feature.set.long <- completeFeatureSet %>% dplyr::select(-timeCount) %>% pivot_longer(cols = c(Bed,LA,LE,LW,RA,RE,RW),names_to = "sensor") %>% mutate(AccelPatientNo_ = patient_clinical_data$AccelPatientNo_[ptIdx]) %>% dplyr::select(-ptIdx)
saveRDS(complete.feature.set.long,file='../all_motion_feature_data/long_feature_set.rds')

## SHUBS CHECKPOINT
# Load missing accelerometry information:
missing.time.info <- read_xlsx('../all_motion_feature_data/MissingPercentTable.xlsx',.name_repair = "universal") %>% 
  mutate(Start.Timestamp = as.POSIXct(Start.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"), End.Timestamp = as.POSIXct(End.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York")) %>% 
  mutate(Recording.Duration = End.Timestamp - Start.Timestamp) %>%
  rename(AccelPatientNo_=Accel.Patient.No.) %>%
  mutate(Bed.miss.hours = Bed*Recording.Duration,LA.miss.hours = LA*Recording.Duration,LE.miss.hours = LE*Recording.Duration,LW.miss.hours = LW*Recording.Duration,RA.miss.hours = RA*Recording.Duration,RE.miss.hours = RE*Recording.Duration,RW.miss.hours = RW*Recording.Duration)

# Load automatically extracted GCS labels:
gcs.data <- read.csv('../clinical_data/clean_auto_GCS_table.csv') %>% select(-X, -Coincide.with.Accel.Recording) %>% mutate(TakenInstant = as.POSIXct(TakenInstant, tz = "America/New_York")) %>% mutate(TakenDay = as.Date(TakenInstant, tz = "America/New_York")) %>% arrange(AccelPatientNo_,TakenInstant)
gcs.data$keep.detection <- FALSE
gcs.data$keep.prediction.lead.0 <- FALSE
gcs.data$keep.prediction.lead.1 <- FALSE
gcs.data$keep.prediction.lead.2 <- FALSE
gcs.data$keep.prediction.lead.6 <- FALSE

# Load long-form NA free complete feature set
na.free.complete.feature.set.long <- readRDS('../all_motion_feature_data/long_feature_set.rds') %>% drop_na(value)
na.free.complete.feature.set.long$keep.detection <- FALSE

# Find GCS data that corresponds to detection and prediction visualization exercises 
for (i in 1:nrow(missing.time.info)){
  curr.start.time <- missing.time.info$Start.Timestamp[i]
  curr.end.time <- missing.time.info$End.Timestamp[i]
  
  keep.detection.idx <- which(gcs.data$AccelPatientNo_ == missing.time.info$AccelPatientNo_[i] & (gcs.data$TakenInstant>=curr.start.time & gcs.data$TakenInstant<=curr.end.time+(6*60*60)))
  keep.prediction.lead.1.idx <- which(gcs.data$AccelPatientNo_ == missing.time.info$AccelPatientNo_[i] & (gcs.data$TakenInstant>=curr.start.time+(1*60*60) & gcs.data$TakenInstant<=curr.end.time+(7*60*60)))
  keep.prediction.lead.2.idx <- which(gcs.data$AccelPatientNo_ == missing.time.info$AccelPatientNo_[i] & (gcs.data$TakenInstant>=curr.start.time+(2*60*60) & gcs.data$TakenInstant<=curr.end.time+(8*60*60)))
  keep.prediction.lead.6.idx <- which(gcs.data$AccelPatientNo_ == missing.time.info$AccelPatientNo_[i] & (gcs.data$TakenInstant>=curr.start.time+(6*60*60) & gcs.data$TakenInstant<=curr.end.time+(12*60*60)))
  
  gcs.data$keep.detection[keep.detection.idx] <- TRUE
  gcs.data$keep.prediction.lead.0[keep.detection.idx] <- TRUE
  gcs.data$keep.prediction.lead.1[keep.prediction.lead.1.idx] <- TRUE
  gcs.data$keep.prediction.lead.2[keep.prediction.lead.2.idx] <- TRUE
  gcs.data$keep.prediction.lead.6[keep.prediction.lead.6.idx] <- TRUE
}

# Detection Plots:
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
  merged.filt.feat.set <- left_join(filt.feat.set,data.frame(time.seq,time.seq.idx), by = c("timeStamps" = "time.seq")) %>% mutate(GCST = detection.gcs.data$Glasgow.Coma.Scale.Score[i], GCSm = detection.gcs.data$Best.Motor.Response[i], GCSv = detection.gcs.data$Best.Verbal.Response[i], GCSe = detection.gcs.data$Eye.Opening[i])
  detection.plot.df <- rbind(detection.plot.df, merged.filt.feat.set)
  names(detection.plot.df) <- names(merged.filt.feat.set)
  print(paste("Detection row",i,"out of",nrow(detection.gcs.data),"completed"))
}
saveRDS(object = detection.plot.df, file = '../all_motion_feature_data/detection_plot_df.rds')

detection.plot.df <- readRDS('../all_motion_feature_data/detection_plot_df.rds') %>% dplyr::select(-keep.detection)

GCSm.detection.grouping <- detection.plot.df %>% drop_na(GCSm) %>% group_by(GCSm,featureType,sensor,time.seq.idx) %>% summarise(mean.feat.value = mean(value), sd.feat.value = sd(value), lower.bound = max((mean(value) - sd(value)),0), upper.bound = mean(value) + sd(value))
saveRDS(object = GCSm.detection.grouping, file = '../all_motion_feature_data/gcsm_detection_group_df.rds')

GCSm.detection.grouping <- readRDS('../all_motion_feature_data/gcsm_detection_group_df.rds') %>% 
  mutate(sensor = as.factor(sensor), GCSm = as.factor(GCSm))

source('./functions/movingAverage.R')

GCSm.detection.grouping$filt.mean.feat.values <- NA
GCSm.detection.grouping$filt.upper.bound <- NA
GCSm.detection.grouping$filt.lower.bound <- NA

for (curr.feat.type in unique(GCSm.detection.grouping$featureType)){
  for (curr.sensor in unique(GCSm.detection.grouping$sensor)){
    for (curr.GCSm in unique(GCSm.detection.grouping$GCSm)){
      
      row.idx <- GCSm.detection.grouping$featureType == curr.feat.type & GCSm.detection.grouping$sensor ==  curr.sensor & GCSm.detection.grouping$GCSm ==  curr.GCSm

      curr.time.seq.idx <- GCSm.detection.grouping$time.seq.idx[row.idx]
      curr.mean.feat.value <- GCSm.detection.grouping$mean.feat.value[row.idx]
      curr.upper.bound <- GCSm.detection.grouping$upper.bound[row.idx] 
      curr.lower.bound <- GCSm.detection.grouping$lower.bound[row.idx]
      
      curr.df <- data.frame(curr.time.seq.idx, curr.mean.feat.value, curr.upper.bound, curr.lower.bound) %>% arrange(curr.time.seq.idx)
      
      GCSm.detection.grouping$filt.mean.feat.values[row.idx] <- movingAverage(curr.df$curr.mean.feat.value, n=(10*12+1),fill.type = NA)
      GCSm.detection.grouping$filt.upper.bound[row.idx] <- movingAverage(curr.df$curr.upper.bound, n=(10*12+1),fill.type = NA)
      GCSm.detection.grouping$filt.lower.bound[row.idx] <- movingAverage(curr.df$curr.lower.bound, n=(10*12+1),fill.type = NA)
    }
  }
}

GCSm.plots <- GCSm.detection.grouping %>% filter(sensor != 'Bed') %>% ggplot(aes(x = time.seq.idx, y= filt.mean.feat.values)) +
  facet_grid(rows = vars(featureType), cols = vars(sensor), scales = "free") +
  geom_ribbon(aes(ymin = filt.lower.bound, ymax = filt.upper.bound, fill = GCSm), alpha=0.2) +
  geom_line(aes(color = GCSm),size=.75)
  
GCSe.detection.grouping <- detection.plot.df %>% drop_na(GCSe) %>% group_by(GCSe,featureType,sensor,time.seq.idx) %>% summarise(mean.feat.value = mean(value), sd.feat.value = sd(value), lower.bound = max((mean(value) - sd(value)),0), upper.bound = mean(value) + sd(value)) %>% 
  mutate(sensor = as.factor(sensor), GCSe = as.factor(GCSe))

# Prediction Lead 1 Plots:
prediction.lead.1.plot.df <- as.data.frame(matrix(nrow = 0,ncol = 10))
prediction.lead.1.gcs.data <- gcs.data %>% filter(keep.prediction.lead.1 == TRUE)

for (i in 1:nrow(prediction.lead.1.gcs.data)){
  print(paste("prediction.lead.1 row",i,"out of",nrow(prediction.lead.1.gcs.data),"started"))
  filt.feat.set <- na.free.complete.feature.set.long %>% filter(AccelPatientNo_ == prediction.lead.1.gcs.data$AccelPatientNo_[i] & timeStamps >= prediction.lead.1.gcs.data$TakenInstant[i] - (7*60*60) & timeStamps <= prediction.lead.1.gcs.data$TakenInstant[i] - (1*60*60))
  if (nrow(filt.feat.set) == 0){
    print(paste("prediction.lead.1 row",i,"out of",nrow(prediction.lead.1.gcs.data),"skipped"))
    next
  }
  time.seq <- seq(from = prediction.lead.1.gcs.data$TakenInstant[i] - (7*60*60),to = prediction.lead.1.gcs.data$TakenInstant[i] - (1*60*60),by= 5)
  time.seq.idx <- 1:length(time.sequence)
  merged.filt.feat.set <- left_join(filt.feat.set,data.frame(time.seq,time.seq.idx), by = c("timeStamps" = "time.seq")) %>% mutate(GCST = prediction.lead.1.gcs.data$Glasgow.Coma.Scale.Score[i], GCSm = prediction.lead.1.gcs.data$Best.Motor.Response[i], GCSv = prediction.lead.1.gcs.data$Best.Verbal.Response[i], GCSe = prediction.lead.1.gcs.data$Eye.Opening[i])
  prediction.lead.1.plot.df <- rbind(prediction.lead.1.plot.df, merged.filt.feat.set)
  names(prediction.lead.1.plot.df) <- names(merged.filt.feat.set)
  print(paste("prediction.lead.1 row",i,"out of",nrow(prediction.lead.1.gcs.data),"completed"))
}

saveRDS(object = prediction.lead.1.plot.df, file = '../all_motion_feature_data/prediction_lead_1_plot_df.rds')


# Examnine GCS scores of each patient 
# Load automatically extracted GCS labels:
gcs_data <- read.csv('../clinical_data/clean_auto_GCS_table.csv') %>% select(-X) %>% mutate(TakenInstant = as.POSIXct(TakenInstant, tz = "America/New_York"))

# Merge GCS.data and patient clinical data
merged.gcs.data <- left_join(gcs_data, patient_clinical_data, by = "AccelPatientNo_") %>% mutate(TakenDay = as.Date(TakenInstant, tz = "America/New_York")) %>% arrange(AccelPatientNo_,TakenInstant)

## Steps:
#-calculate number of evals per day per patient
summ.gcs.stats <- merged.gcs.data %>% filter(TakenDay >= NCCUAdmissionDate & TakenDay <= NCCUDischargeDate) %>% group_by(AccelPatientNo_, TakenDay) %>% summarise(no.evals = n())
#-filter out GCS scores that coincide with accelerometry recording time
accel.merged.gcs.data <- left_join(gcs_data, missingTimeInfo, by = "AccelPatientNo_") %>% arrange(AccelPatientNo_,TakenInstant) %>% filter(TakenInstant >= Start.Timestamp & TakenInstant <= End.Timestamp)
#-worst GCSm within 24 hours of admission
worst.gcsm.nccu.admission <- merged.gcs.data %>% filter(TakenDay >= NCCUAdmissionDate & TakenDay <= NCCUAdmissionDate+1) %>% group_by(AccelPatientNo_) %>% summarise(worst.GCSm = min(Best.Motor.Response,na.rm = TRUE))