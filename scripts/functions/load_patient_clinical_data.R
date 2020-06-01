load_patient_clinical_data <- function(directory) {
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
      Death12Months = as.factor(Death12Months),
      GOSE_fav_12mo = as.factor(GOSE_fav_12mo),
      mRS_fav_12mo = as.factor(mRS_fav_12mo)
    )
  return(patient_clinical_data)
}