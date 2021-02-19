#### Clinical Data Extraction and Examination ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# University of Cambridge
# Johns Hopkins University
# email address: sb2406@cam.ac.uk


library(tidyverse)
library(readxl)
library(R.matlab)
library(scales)

# Load patient clinical data
source('./functions/load_patient_clinical_data.R')
patient.clinical.data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv') %>% arrange(AccelPatientNo_) %>% mutate(ptIdx = 1:nrow(.))

# Load automatically extracted GCS labels with indicator labels:
gcs.data <- read.csv('../clinical_data/03_gcs_data_w_indicators.csv')

# Load missing time info (includes information on recording start and stop times)
missing.time.info <- read_xlsx('../all_motion_feature_data/MissingPercentTable.xlsx',.name_repair = "universal") %>% 
    mutate(Start.Timestamp = as.POSIXct(Start.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"), End.Timestamp = as.POSIXct(End.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York")) %>% 
    mutate(Recording.Duration = End.Timestamp - Start.Timestamp) %>%
    rename(AccelPatientNo_=Accel.Patient.No.) %>%
    mutate(Bed.miss.hours = Bed*Recording.Duration,LA.miss.hours = LA*Recording.Duration,LE.miss.hours = LE*Recording.Duration,LW.miss.hours = LW*Recording.Duration,RA.miss.hours = RA*Recording.Duration,RE.miss.hours = RE*Recording.Duration,RW.miss.hours = RW*Recording.Duration)

# Calculate gross static activity percentages by sensor
static.activity.percs <- miss.NM.TOD %>%
  group_by(sensor) %>%
  summarise(total.n = sum(n), total.NM = sum(n) - sum(num.NM) - sum(num.Missing),totalMissing  = sum(num.Missing)) %>%
  mutate(perc.active = sprintf('%.02f',100*total.NM/total.n), perc.Missing = sprintf('%.02f',100*totalMissing/total.n))

# Calculations of means and spreads for table 3
feature.type.mean.results <- data.frame(matrix(ncol = 9, nrow = 0))
kruskal.test.results <- data.frame(matrix(ncol = 6, nrow = 0))
imp.dirs <- list.files('../all_motion_feature_data/04_formatted_predictor_matrices/',pattern = 'imp*',include.dirs = TRUE)
for (i in imp.dirs){
  print(paste(i,'started.'))
  detection.folders <- list.files(file.path('../all_motion_feature_data/04_formatted_predictor_matrices',i),pattern='detection_window_*',include.dirs = TRUE)
  for (j in detection.folders){
    print(paste(j,'started.'))
    long.form.train.df <- readRDS(file.path('../all_motion_feature_data/04_formatted_predictor_matrices',i,j,'long_form_motor_train_values.rds')) %>% mutate(partition = 'train')
    long.form.test.df <- readRDS(file.path('../all_motion_feature_data/04_formatted_predictor_matrices',i,j,'long_form_motor_test_values.rds')) %>% mutate(partition = 'test')
    long.form.df <- rbind(long.form.train.df,long.form.test.df)
    rm(long.form.train.df,long.form.test.df)
    gc()
    curr.summary.df <- long.form.df %>% 
      group_by(featureType,GCSm) %>%
      summarise(meanValue = mean(value),sdValue = sd(value), medianValue = median(value), Q1Value =  quantile(value, 0.25), Q3Value =  quantile(value, 0.75)) %>%
      mutate(imputation = i, detection.window = j)
    feature.type.mean.results <- rbind(feature.type.mean.results, curr.summary.df)
    unique.feature.types <- unique(long.form.df$featureType)
    for (k in unique.feature.types){
      print(paste(k,'Kruskal Willis started.'))
      filt.long.form.df <- long.form.df %>% filter(featureType == k)
      curr.kw.test <- kruskal.test(value ~ GCSm, data = filt.long.form.df)
      kruskal.test.results <- rbind(kruskal.test.results,
                                    data.frame(imputation = i,
                                               detection.window = j,
                                               featureType = k,
                                               statistic = curr.kw.test$statistic,
                                               d.f = curr.kw.test$parameter,
                                               p.value = curr.kw.test$p.value))
    }
  }
}

# Save mean results dataframe
saveRDS(feature.type.mean.results,'../all_motion_feature_data/04_formatted_predictor_matrices/detection_feature_means.rds')

# Save Kruskal-Wallis test results
saveRDS(kruskal.test.results,'../all_motion_feature_data/04_formatted_predictor_matrices/detection_kruskal_wallis_test_results.rds')

# Calculate pairwise Wilcoxon mean differences
imp.dirs <- list.files('../all_motion_feature_data/04_formatted_predictor_matrices/',pattern = 'imp*',include.dirs = TRUE)
complete.compiled.wx.p.values <- as.data.frame(matrix(ncol = 8, nrow = 0))

for (i in imp.dirs){
  print(paste('imputation',i,'started.'))
  detection.folders <- list.files(file.path('../all_motion_feature_data/04_formatted_predictor_matrices',i),pattern='detection_window_*',include.dirs = TRUE)
  for (j in detection.folders){
    print(paste('detection folder',j,'started.'))
    long.form.train.df <- readRDS(file.path('../all_motion_feature_data/04_formatted_predictor_matrices',i,j,'long_form_motor_train_values.rds')) %>% mutate(partition = 'train')
    long.form.test.df <- readRDS(file.path('../all_motion_feature_data/04_formatted_predictor_matrices',i,j,'long_form_motor_test_values.rds')) %>% mutate(partition = 'test')
    long.form.df <- rbind(long.form.train.df,long.form.test.df)
    rm(long.form.train.df,long.form.test.df)
    gc()
    
    unique.feature.types <- unique(long.form.df$featureType)
    for (k in unique.feature.types){
      print(paste(k,'started.'))
      
      unique.sensors <- unique(long.form.df$sensor)
      
      for (m in unique.sensors){
        
        print(paste(m,'started.'))
        filt.long.form.df <- long.form.df %>% filter(featureType == k & sensor == m)
        
        greater.wx.test <- pairwise.wilcox.test(filt.long.form.df$value, filt.long.form.df$GCSm,p.adjust.method = "BH",alternative = 'greater')
        lesser.wx.test <- pairwise.wilcox.test(filt.long.form.df$value, filt.long.form.df$GCSm,p.adjust.method = "BH",alternative = 'less')
        two.sided.wx.test <- pairwise.wilcox.test(filt.long.form.df$value, filt.long.form.df$GCSm,p.adjust.method = "BH")
        
        
        temp.df.1 <- as.data.frame(greater.wx.test$p.value)
        temp.df.1$first.class <- rownames(temp.df.1)
        temp.df.1 <- temp.df.1 %>% pivot_longer(cols = -first.class, names_to = 'second.class',values_to = 'p.value') %>% 
          drop_na(p.value) %>% mutate(test = 'greater',
                                      imp = i,
                                      observation.window = j,
                                      feature.type = k,
                                      sensor = m)
        
        temp.df.2 <- as.data.frame(lesser.wx.test$p.value)
        temp.df.2$first.class <- rownames(temp.df.2)
        temp.df.2 <- temp.df.2 %>% pivot_longer(cols = -first.class, names_to = 'second.class',values_to = 'p.value') %>% 
          drop_na(p.value) %>% mutate(test = 'less',
                                      imp = i,
                                      observation.window = j,
                                      feature.type = k,
                                      sensor = m)
        
        temp.df.3 <- as.data.frame(two.sided.wx.test$p.value)
        temp.df.3$first.class <- rownames(temp.df.3)
        temp.df.3 <- temp.df.3 %>% pivot_longer(cols = -first.class, names_to = 'second.class',values_to = 'p.value') %>% 
          drop_na(p.value) %>% mutate(test = 'two.tailed',
                                      imp = i,
                                      observation.window = j,
                                      feature.type = k,
                                      sensor = m)
        
        curr.compiled.wx.p.values <- rbind(temp.df.1,temp.df.2,temp.df.3)
        
        complete.compiled.wx.p.values <- rbind(complete.compiled.wx.p.values,
                                               curr.compiled.wx.p.values)
      }
    }
  }
}

# Save Wilcoxon test results
saveRDS(complete.compiled.wx.p.values,'../all_motion_feature_data/04_formatted_predictor_matrices/detection_wilcoxon_test_results.rds')

# Calculate p-values across imputations (Wilcoxon)
grouped.complete.wx.p.values <- complete.compiled.wx.p.values %>%
  group_by(observation.window,feature.type,sensor,first.class,second.class,test) %>%
  mutate(z.norm = qnorm(p.value, lower.tail=FALSE)) %>%
  summarise(mean.p.value = mean(p.value),
            std.p.value = sd(p.value),
            mean.z.norm = mean(z.norm),
            var.z.norm = 1 + (1 + (1/n()))*var(z.norm),
            r.m = (1 + (1/n()))*var(z.norm),
            m = n()) %>%
  mutate(d.o.f. = (m-1)*((1 + (1/r.m))^2)) %>%
  mutate(p.m = pt(q = mean.z.norm, df = d.o.f., lower.tail = FALSE)) %>%
  mutate(p.value.imp = sprintf("%.2f",p.m))

# two-tailed Wilcoxon tests
two.tailed.non.sig.results <- grouped.complete.wx.p.values %>%
  filter(test == 'two.tailed',p.m >= 0.05)

combos.two.tailed.non.sig.results <- two.tailed.non.sig.results %>%
  group_by(first.class, second.class) %>%
  summarise(n = n())

feature.combos.2.tailed.non.sig <- two.tailed.non.sig.results %>%
  group_by(feature.type) %>%
  summarise(n = n())

sensors.combos.2.tailed.non.sig <- two.tailed.non.sig.results %>%
  group_by(sensor) %>%
  summarise(n = n())

# one-tailed Wilcoxon tests
one.tailed.non.sig.results <- grouped.complete.wx.p.values %>%
  filter(test != 'two.tailed',p.m >= 0.05)

one.tailed.sig.results <- grouped.complete.wx.p.values %>%
  filter(test != 'two.tailed',p.m <= 0.05 | mean.p.value == 0) %>%
  group_by(feature.type,test,first.class, second.class) %>%
  summarise(n = n())

combos.two.tailed.non.sig.results <- two.tailed.non.sig.results %>%
  group_by(first.class, second.class) %>%
  summarise(n = n())

feature.combos.2.tailed.non.sig <- two.tailed.non.sig.results %>%
  group_by(feature.type, first.class, second.class) %>%
  summarise(n = n())

sensors.combos.2.tailed.non.sig <- two.tailed.non.sig.results %>%
  group_by(sensor, first.class, second.class) %>%
  summarise(n = n())

# Load mean results dataframe
feature.type.mean.results <- readRDS('../all_motion_feature_data/04_formatted_predictor_matrices/detection_feature_means.rds')

# Load Kruskal-Wallis test results
kruskal.test.results <- readRDS('../all_motion_feature_data/04_formatted_predictor_matrices/detection_kruskal_wallis_test_results.rds')

# Calculate pooled (across imputations) means
totaled.means <- feature.type.mean.results %>% 
  group_by(featureType, detection.window, GCSm) %>% 
  summarise(pooled.mean = mean(meanValue), pooled.sd = sqrt( (1/n())*sum(sdValue^2)  + ((n()+1)/(n()))*(sd(meanValue)^2) )) 

sma.means.formatted <- totaled.means  %>% 
  filter(featureType == "sma") %>%
  mutate(formatted = sprintf("%.2f (%.2f)",pooled.mean,pooled.sd))

hlf.h.means.formatted <- totaled.means  %>% 
  filter(featureType == "freq_pairs1") %>%
  mutate(formatted = paste0(formatC(pooled.mean,digits = 2,format = 'E'),' (',formatC(pooled.sd,digits = 2,format = 'E'),')'))
           
hlf.l.means.formatted <- totaled.means  %>% 
  filter(featureType == "freq_pairs2") %>%
  mutate(formatted = paste0(formatC(pooled.mean,digits = 2,format = 'E'),' (',formatC(pooled.sd,digits = 2,format = 'E'),')'))

mfr.means.formatted <- totaled.means  %>% 
  filter(featureType == "med_freq") %>%
  mutate(formatted = sprintf("%.2f (%.2f)",pooled.mean,pooled.sd))

fde.means.formatted <- totaled.means  %>% 
  filter(featureType == "freq_entropy") %>%
  mutate(formatted = sprintf("%.2f (%.2f)",pooled.mean,pooled.sd))

bpw.means.formatted <- totaled.means  %>% 
  filter(featureType == "band_power") %>%
  mutate(formatted = paste0(formatC(pooled.mean,digits = 2,format = 'E'),' (',formatC(pooled.sd,digits = 2,format = 'E'),')'))

wvl.means.formatted <- totaled.means  %>% 
  filter(featureType == "wavelets") %>%
  mutate(formatted = paste0(formatC(pooled.mean,digits = 2,format = 'E'),' (',formatC(pooled.sd,digits = 2,format = 'E'),')'))

# Calculations of total (across GCS) means and spreads for table 3
total.means.feature.type.results <- data.frame(matrix(ncol = 5, nrow = 0))
for (i in imp.dirs){
  print(paste(i,'started.'))
  detection.folders <- list.files(file.path('../all_motion_feature_data/04_formatted_predictor_matrices',i),pattern='detection_window_*',include.dirs = TRUE)
  for (j in detection.folders){
    print(paste(j,'started.'))
    long.form.train.df <- readRDS(file.path('../all_motion_feature_data/04_formatted_predictor_matrices',i,j,'long_form_motor_train_values.rds')) %>% mutate(partition = 'train')
    long.form.test.df <- readRDS(file.path('../all_motion_feature_data/04_formatted_predictor_matrices',i,j,'long_form_motor_test_values.rds')) %>% mutate(partition = 'test')
    long.form.df <- rbind(long.form.train.df,long.form.test.df)
    rm(long.form.train.df,long.form.test.df)
    gc()
    curr.summary.df <- long.form.df %>% 
      group_by(featureType) %>%
      summarise(meanValue = mean(value),sdValue = sd(value)) %>%
      mutate(imputation = i, detection.window = j)
    total.means.feature.type.results <- rbind(total.means.feature.type.results, curr.summary.df)
  }
}

# Save total mean results dataframe
saveRDS(total.means.feature.type.results,'../all_motion_feature_data/04_formatted_predictor_matrices/detection_feature_means_across_GCSm.rds')

# Group by imputation and pool
pooled.total.means <- total.means.feature.type.results %>% 
  group_by(featureType, detection.window) %>% 
  summarise(pooled.mean = mean(meanValue), pooled.sd = sqrt( (1/n())*sum(sdValue^2)  + ((n()+1)/(n()))*(sd(meanValue)^2) ))

sma.means.formatted <- pooled.total.means  %>% 
  filter(featureType == "sma") %>%
  mutate(formatted = sprintf("%.2f (%.2f)",pooled.mean,pooled.sd))

hlf.h.means.formatted <- pooled.total.means  %>% 
  filter(featureType == "freq_pairs1") %>%
  mutate(formatted = paste0(formatC(pooled.mean,digits = 2,format = 'E'),' (',formatC(pooled.sd,digits = 2,format = 'E'),')'))

hlf.l.means.formatted <- pooled.total.means  %>% 
  filter(featureType == "freq_pairs2") %>%
  mutate(formatted = paste0(formatC(pooled.mean,digits = 2,format = 'E'),' (',formatC(pooled.sd,digits = 2,format = 'E'),')'))

mfr.means.formatted <- pooled.total.means  %>% 
  filter(featureType == "med_freq") %>%
  mutate(formatted = sprintf("%.2f (%.2f)",pooled.mean,pooled.sd))

fde.means.formatted <- pooled.total.means  %>% 
  filter(featureType == "freq_entropy") %>%
  mutate(formatted = sprintf("%.2f (%.2f)",pooled.mean,pooled.sd))

bpw.means.formatted <- pooled.total.means  %>% 
  filter(featureType == "band_power") %>%
  mutate(formatted = paste0(formatC(pooled.mean,digits = 2,format = 'E'),' (',formatC(pooled.sd,digits = 2,format = 'E'),')'))

wvl.means.formatted <- pooled.total.means  %>% 
  filter(featureType == "wavelets") %>%
  mutate(formatted = paste0(formatC(pooled.mean,digits = 2,format = 'E'),' (',formatC(pooled.sd,digits = 2,format = 'E'),')'))

## CARRYOVER FROM OTHER FILES
# Load missing accelerometry information:
missing.time.info <- read_xlsx('../all_motion_feature_data/MissingPercentTable.xlsx',.name_repair = "universal") %>% 
  mutate(Start.Timestamp = as.POSIXct(Start.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"), End.Timestamp = as.POSIXct(End.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York")) %>% 
  mutate(Recording.Duration = End.Timestamp - Start.Timestamp) %>%
  rename(AccelPatientNo_=Accel.Patient.No.) %>%
  mutate(
    Bed.miss.hours = Bed * Recording.Duration,
    LA.miss.hours = LA * Recording.Duration,
    LE.miss.hours = LE * Recording.Duration,
    LW.miss.hours = LW * Recording.Duration,
    RA.miss.hours = RA * Recording.Duration,
    RE.miss.hours = RE * Recording.Duration,
    RW.miss.hours = RW * Recording.Duration
  )

# Load patient clinical data
source('./functions/load_patient_clinical_data.R')
patient.clinical.data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv') %>% arrange(AccelPatientNo_) %>% mutate(ptIdx = 1:nrow(.))

# Load indices of totally missing series to exclude from analysis
complete.missing.idx <- read.csv('../all_motion_feature_data/complete_missing_idx.csv')

# Load no motion accelerometry information:
static.time.info <- read_xlsx('../all_motion_feature_data/NoMotionPercentTable.xlsx',.name_repair = "universal") %>% 
  mutate(Recording.Duration = missing.time.info$End.Timestamp - missing.time.info$Start.Timestamp) %>%
  rename(AccelPatientNo_=Accel.Patient.No.)
# Correct Totally Missing Sets in No Motion Info:
for (i in 1:nrow(complete.missing.idx)){
  currPtIdx <- complete.missing.idx$ptIdx[i]
  currSrIdx <- complete.missing.idx$srIdx[i]
  static.time.info[currPtIdx,currSrIdx+1] <- NA
}

static.time.info <- static.time.info %>%
  rowwise() %>%
  mutate(
    average.extremity.nm = mean(LA,LE,LW,RA,RE,RW, na.rm = TRUE),
    Bed.nm.hours = Bed * Recording.Duration,
    LA.nm.hours = LA * Recording.Duration,
    LE.nm.hours = LE * Recording.Duration,
    LW.nm.hours = LW * Recording.Duration,
    RA.nm.hours = RA * Recording.Duration,
    RE.nm.hours = RE * Recording.Duration,
    RW.nm.hours = RW * Recording.Duration
  )

gose.info <- patient.clinical.data %>%
  dplyr::select(AccelPatientNo_,GOS_EDischarge,GOS_E12Months)

# Bind GOSE info to no motion data
static.time.info <- left_join(static.time.info,gose.info,by = 'AccelPatientNo_')

# GOSE-specific static info
gose.discharge.static <- static.time.info %>%
  group_by(GOS_EDischarge) %>%
  summarise(n = sum(!is.na(average.extremity.nm)),
            mean.extremity.nm = 100*mean(average.extremity.nm,na.rm = TRUE),
            sd.extremity.nm = 100*sd(average.extremity.nm,na.rm = TRUE),
            median.extremity.nm = 100*median(average.extremity.nm,na.rm = TRUE),
            q1.extremity.nm = 100*quantile(average.extremity.nm,.25,na.rm = TRUE),
            q3.extremity.nm =100*quantile(average.extremity.nm,.75,na.rm = TRUE)) %>%
  mutate(formatted = sprintf('%0.02f (%0.02f - %0.02f)',median.extremity.nm,q1.extremity.nm,q3.extremity.nm))

gose.12mo.static <- static.time.info %>%
  group_by(GOS_E12Months) %>%
  summarise(n = sum(!is.na(average.extremity.nm)),
            mean.extremity.nm = 100*mean(average.extremity.nm,na.rm = TRUE),
            sd.extremity.nm = 100*sd(average.extremity.nm,na.rm = TRUE),
            median.extremity.nm = 100*median(average.extremity.nm,na.rm = TRUE),
            q1.extremity.nm = 100*quantile(average.extremity.nm,.25,na.rm = TRUE),
            q3.extremity.nm =100*quantile(average.extremity.nm,.75,na.rm = TRUE)) %>%
  mutate(formatted = sprintf('%0.02f (%0.02f - %0.02f)',median.extremity.nm,q1.extremity.nm,q3.extremity.nm))

# GCSm-specific static info:
gcs.data <- read.csv('../clinical_data/03_gcs_data_w_indicators.csv') %>% 
  mutate(TakenInstant = as.POSIXct(TakenInstant, tz = "America/New_York"))

grouped.gcsm.data <- gcs.data %>%
  select(AccelPatientNo_,Best.Motor.Response,det.indicator) %>%
  filter(det.indicator == 'TRUE') %>%
  group_by(AccelPatientNo_) %>%
  summarise(med.GCSm = floor(median(Best.Motor.Response,na.rm = TRUE)))

# Bind GCSm info to no motion data
static.time.info <- left_join(static.time.info,grouped.gcsm.data,by = 'AccelPatientNo_')

gcsm.specific.static <- static.time.info %>%
  group_by(med.GCSm) %>%
  summarise(n = sum(!is.na(average.extremity.nm)),
            mean.extremity.nm = 100*mean(average.extremity.nm,na.rm = TRUE),
            sd.extremity.nm = 100*sd(average.extremity.nm,na.rm = TRUE),
            median.extremity.nm = 100*median(average.extremity.nm,na.rm = TRUE),
            q1.extremity.nm = 100*quantile(average.extremity.nm,.25,na.rm = TRUE),
            q3.extremity.nm =100*quantile(average.extremity.nm,.75,na.rm = TRUE)) %>%
  mutate(formatted = sprintf('%0.02f (%0.02f - %0.02f)',median.extremity.nm,q1.extremity.nm,q3.extremity.nm))

# Correlations with median GCSm, GOSE at discharge, and GOSE at 12 months
gcsm.correlation <-cor.test(static.time.info$average.extremity.nm, static.time.info$med.GCSm,  method = "spearman", exact = FALSE)
gose.dis.correlation <-cor.test(static.time.info$average.extremity.nm, static.time.info$GOS_EDischarge,  method = "spearman", exact = FALSE)
gose.12m.correlation <-cor.test(static.time.info$average.extremity.nm, static.time.info$GOS_E12Months,  method = "spearman", exact = FALSE)

# LA percent nomotion
LA.perc.nomotion <- static.time.info %>% drop_na(LA)
perc.LA.motion <- 100*as.numeric(sum(LA.perc.nomotion$LA.nm.hours))/as.numeric(sum(LA.perc.nomotion$Recording.Duration))

# LE percent nomotion
LE.perc.nomotion <- static.time.info %>% drop_na(LE)
perc.LE.motion <- 100*as.numeric(sum(LE.perc.nomotion$LE.nm.hours))/as.numeric(sum(LE.perc.nomotion$Recording.Duration))

# LW percent nomotion
LW.perc.nomotion <- static.time.info %>% drop_na(LW)
perc.LW.motion <- 100*as.numeric(sum(LW.perc.nomotion$LW.nm.hours))/as.numeric(sum(LW.perc.nomotion$Recording.Duration))

# RA percent nomotion
RA.perc.nomotion <- static.time.info %>% drop_na(RA)
perc.RA.motion <- 100*as.numeric(sum(RA.perc.nomotion$RA.nm.hours))/as.numeric(sum(RA.perc.nomotion$Recording.Duration))

# RE percent nomotion
RE.perc.nomotion <- static.time.info %>% drop_na(RE)
perc.RE.motion <- 100*as.numeric(sum(RE.perc.nomotion$RE.nm.hours))/as.numeric(sum(RE.perc.nomotion$Recording.Duration))

# RW percent nomotion
RW.perc.nomotion <- static.time.info %>% drop_na(RW)
perc.RW.motion <- 100*as.numeric(sum(RW.perc.nomotion$RW.nm.hours))/as.numeric(sum(RW.perc.nomotion$Recording.Duration))

# Bed percent nomotion
Bed.perc.nomotion <- static.time.info %>% drop_na(Bed) %>% filter(AccelPatientNo_ != 64)
perc.Bed.motion <- 100*as.numeric(sum(Bed.perc.nomotion$Bed.nm.hours))/as.numeric(sum(Bed.perc.nomotion$Recording.Duration))

# Total percent nomotion
Total.hours.motion <- as.numeric(sum(LA.perc.nomotion$LA.nm.hours)) + as.numeric(sum(LE.perc.nomotion$LE.nm.hours)) + as.numeric(sum(LW.perc.nomotion$LW.nm.hours)) + as.numeric(sum(RA.perc.nomotion$RA.nm.hours)) + as.numeric(sum(RE.perc.nomotion$RE.nm.hours)) + as.numeric(sum(RW.perc.nomotion$RW.nm.hours))
Total.hours.of.set <- as.numeric(sum(LA.perc.nomotion$Recording.Duration)) + as.numeric(sum(LE.perc.nomotion$Recording.Duration)) + as.numeric(sum(LW.perc.nomotion$Recording.Duration)) + as.numeric(sum(RA.perc.nomotion$Recording.Duration)) + as.numeric(sum(RE.perc.nomotion$Recording.Duration)) + as.numeric(sum(RW.perc.nomotion$Recording.Duration))
perc.total.motion <-  100*Total.hours.motion/Total.hours.of.set

Total.static.time.info <- static.time.info %>% mutate(num.col = (rowSums(!is.na(.))-2)/2) %>% mutate(total.hours = Recording.Duration*num.col)
Total.static.time.info$Total.Hours.of.Motion <- rowSums(sapply(static.time.info[,10:16], as.numeric),na.rm = TRUE)
Total.static.time.info$Total.Percs.of.Motion <- 100*Total.static.time.info$Total.Hours.of.Motion/as.numeric(Total.static.time.info$total.hours)

complete.feature.set.long <- completeFeatureSet %>% dplyr::select(-timeCount) %>% pivot_longer(cols = c(Bed,LA,LE,LW,RA,RE,RW),names_to = "sensor") %>% mutate(AccelPatientNo_ = patient.clinical.data$AccelPatientNo_[ptIdx]) %>% dplyr::select(-ptIdx)
saveRDS(complete.feature.set.long,file='../all_motion_feature_data/long_feature_set.rds')

## SHUBS CHECKPOINT
# Load patient clinical data
source('./functions/load_patient_clinical_data.R')
patient.clinical.data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv') %>% arrange(AccelPatientNo_) %>% mutate(ptIdx = 1:nrow(.))

# Load missing accelerometry information:
missing.time.info <- read_xlsx('../all_motion_feature_data/MissingPercentTable.xlsx',.name_repair = "universal") %>% 
  mutate(Start.Timestamp = as.POSIXct(Start.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York"), End.Timestamp = as.POSIXct(End.Timestamp,format = '%d-%b-%Y %H:%M:%S',tz = "America/New_York")) %>% 
  mutate(Recording.Duration = End.Timestamp - Start.Timestamp) %>%
  rename(AccelPatientNo_=Accel.Patient.No.) %>%
  mutate(Bed.miss.hours = Bed*Recording.Duration,LA.miss.hours = LA*Recording.Duration,LE.miss.hours = LE*Recording.Duration,LW.miss.hours = LW*Recording.Duration,RA.miss.hours = RA*Recording.Duration,RE.miss.hours = RE*Recording.Duration,RW.miss.hours = RW*Recording.Duration)

# Examnine GCS scores of each patient 
# Load automatically extracted GCS labels:
gcs_data <- read.csv('../clinical_data/02_clean_auto_GCS_table.csv') %>% select(-X) %>% mutate(TakenInstant = as.POSIXct(TakenInstant, tz = "America/New_York"))

# Merge GCS.data and patient clinical data
merged.gcs.data <- left_join(gcs_data, patient.clinical.data, by = "AccelPatientNo_") %>% mutate(TakenDay = as.Date(TakenInstant, tz = "America/New_York")) %>% arrange(AccelPatientNo_,TakenInstant)

## Steps:
#-calculate number of evals per day per patient
summ.gcs.stats <- merged.gcs.data %>% filter(TakenDay >= NCCUAdmissionDate & TakenDay <= NCCUDischargeDate) %>% group_by(AccelPatientNo_, TakenDay) %>% summarise(no.evals = n())
#-filter out GCS scores that coincide with accelerometry recording time
accel.merged.gcs.data <- left_join(gcs_data, missingTimeInfo, by = "AccelPatientNo_") %>% arrange(AccelPatientNo_,TakenInstant) %>% filter(TakenInstant >= Start.Timestamp & TakenInstant <= End.Timestamp)
#-worst GCSm within 24 hours of admission
worst.gcsm.nccu.admission <- merged.gcs.data %>% filter(TakenDay >= NCCUAdmissionDate & TakenDay <= NCCUAdmissionDate+1) %>% group_by(AccelPatientNo_) %>% summarise(worst.GCSm = min(Best.Motor.Response,na.rm = TRUE))
