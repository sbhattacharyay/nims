#### Clinical Data Extraction and Examination ####
#
# Eshan Joshi, Shubhayu Bhattacharyay
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

library(tidyverse)

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