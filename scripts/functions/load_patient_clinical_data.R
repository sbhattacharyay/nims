load_patient_clinical_data <- function(directory,gose_thresh=5,mrs_thresh=3) {
  na_strings <-
    c("NA","N A","N / A","N/A","N/ A","Not Available","Not available","NaN","","NaT")
  patient_clinical_data <-
    read.csv(directory, na.strings = na_strings) %>%
    arrange(AccelPatientNo_) %>%
    mutate(
      Sex = as.factor(Sex),
      ConsentDate = as.Date(ConsentDate, '%d-%b-%Y'),
      AccelRecordingStartDate = as.Date(AccelRecordingStartDate, '%d-%b-%Y'),
      AccelDataMountDate = as.Date(AccelDataMountDate, '%d-%b-%Y'),
      Hospital_EDAdmissionDate = as.Date(Hospital_EDAdmissionDate, '%d-%b-%Y'),
      NCCUAdmissionDate = as.Date(NCCUAdmissionDate, '%d-%b-%Y'),
      DOB = as.Date(DOB, '%d-%b-%Y'),
      NCCUDischargeDate = as.Date(NCCUDischargeDate, '%d-%b-%Y'),
      HospitalDischargeDate = as.Date(HospitalDischargeDate, '%d-%b-%Y'),
      DeathDate = as.Date(DeathDate, '%d-%b-%Y'),
      ActiveProblem1 = as.character(ActiveProblem1),
      ActiveProblem2 = as.character(ActiveProblem2),
      ActiveProblem3 = as.character(ActiveProblem3),
      ActiveProblem4 = as.character(ActiveProblem4),
      ActiveProblem5 = as.character(ActiveProblem5),
      CVA = as.factor(CVA),
      ICH = as.factor(ICH),
      SAH = as.factor(SAH),
      BrainTumorOrLesion = as.factor(BrainTumorOrLesion),
      SDHOrEDH = as.factor(SDHOrEDH),
      TBI = as.factor(TBI),
      DiedDuringThisHospitalStay_ = as.factor(DiedDuringThisHospitalStay_),
      Destination = as.factor(Destination),
      Death12Months = as.factor(Death12Months)
    ) %>%
    mutate(
      fav_mort_dis = factor(DiedDuringThisHospitalStay_ == "0", labels = c("Unfav", "Fav")),
      fav_mort_12m = factor(Death12Months == "0", labels = c("Unfav", "Fav")),
      fav_GOSE_dis = factor(GOS_EDischarge >= gose_thresh, labels = c("Unfav", "Fav")),
      fav_GOSE_12m = factor(GOS_E12Months >= gose_thresh, labels = c("Unfav", "Fav")),
      fav_mRS_dis = factor(mRSDischarge <= mrs_thresh, labels = c("Unfav", "Fav")),
      fav_mRS_12m = factor(mRS12Months <= mrs_thresh, labels = c("Unfav", "Fav"))
    )
  
  clinical_variable_list <-
    data.frame(
      colnames(patient_clinical_data),
      sapply(patient_clinical_data, class),
      sapply(patient_clinical_data, typeof)    
    ) %>%
    rename(
      variable = 1,
      class = 2,
      type = 3
    )
  rownames(clinical_variable_list) <- c()
  write.csv(clinical_variable_list,
            '../clinical_data/clinical_variable_list.csv')
  
  return(patient_clinical_data)
}