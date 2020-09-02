library(tidyverse)
library(readxl)

source('./functions/load_patient_clinical_data.R')

gcs_table <- read_xlsx('../clinical_data/automatic_GCS_table.xlsx') %>% drop_na(Value) %>% pivot_wider(id_cols = c(TakenInstant,PrimaryMrn,BirthDate,LastName,FirstName), names_from = DisplayName, values_from = Value)

patient_clinical_data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv')

accel_no_match <- patient_clinical_data %>% select(StudyPatientNo_,AccelPatientNo_,DOB)

id_gcs_table <- left_join(gcs_table,accel_no_match,by = c("BirthDate" = "DOB")) %>% drop_na(StudyPatientNo_)

missingTable <- read_xlsx('../all_motion_feature_data/MissingPercentTable.xlsx',.name_repair = "universal") %>% mutate(Start.Timestamp = as.POSIXct(Start.Timestamp,format = '%d-%b-%Y %H:%M:%S'), End.Timestamp = as.POSIXct(End.Timestamp,format = '%d-%b-%Y %H:%M:%S'))

# Add indicator of whether the data point coincides with accelerometry recording:

id_gcs_table$Coincide.with.Accel.Recording <- NA

for (i in 1:nrow(missingTable)){
  currAccelPatientNo <- missingTable$Study.Patient.No.[i]
  accelPatientIdx <- which(id_gcs_table$AccelPatientNo_ == currAccelPatientNo)
  
  filtGCS <- id_gcs_table[accelPatientIdx,]
  
  
  id_gcs_table$Coincide.with.Accel.Recording[accelPatientIdx] = filtGCS$TakenInstant <= missingTable$End.Timestamp[i] & filtGCS$TakenInstant >= missingTable$Start.Timestamp[i]
}

finalGCSTable <- id_gcs_table[,c(10,11,1,6,7,8,9,12)]
names(finalGCSTable) <- make.names(names(finalGCSTable))

write.csv(finalGCSTable,'../clinical_data/clean_auto_GCS_table.csv')

finalGCSTable <- read.csv('../clinical_data/clean_auto_GCS_table.csv') %>% select(-X)

man_table <- read_xlsx('../clinical_data/manual_GCS_table.xlsx') %>% mutate(Coincide.with.Accel.Recording = `Coincide with accel` == 1,Time...3 = sprintf("%04d", Time...3)) %>% filter(`Study Patient Number` != 1) %>% mutate(combTime = paste(as.character(Date),Time...3)) %>% mutate(TakenInstant = as.POSIXct(combTime,format = '%Y-%m-%d %H%M'))
man_table <- man_table[,c(1,4:7,10,12)]

names(man_table)[1:5] <- c("StudyPatientNo_","Eye.Opening","Best.Verbal.Response","Best.Motor.Response","Glasgow.Coma.Scale.Score")
accel_no_match <- patient_clinical_data %>% select(StudyPatientNo_,AccelPatientNo_)

comb_man_gcs_table <- left_join(man_table,accel_no_match,by = "StudyPatientNo_") %>% drop_na(AccelPatientNo_)
final_man_gcs_table <- comb_man_gcs_table[,c(1,8,7,5,4,3,2,6)]
final_man_gcs_table <- final_man_gcs_table[order(final_man_gcs_table$TakenInstant),]

write.csv(final_man_gcs_table,'../clinical_data/clean_manu_GCS_table.csv')