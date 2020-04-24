load_patient_clinical_data <- function() {
  na_strings <-
    c("NA",
      "N A",
      "N / A",
      "N/A",
      "N/ A",
      "Not Available",
      "Not available",
      "NaN",
      "",
      "NaT")
  patient_clinical_data <-
    read.csv('../clinical_data/patient_clinical_data.csv', na.strings = na_strings) %>%
    arrange(PY_pNum) %>%
    mutate(
      gender = as.factor(gender),
      date_admission = as.Date(date_admission, '%d-%b-%Y'),
      date_discharge = as.Date(date_discharge, '%d-%b-%Y')
    ) %>%
    mutate(
      diagnosis_1 = as.character(diagnosis_1),
      diagnosis_2 = as.character(diagnosis_2),
      diagnosis_3 = as.character(diagnosis_3)
    ) %>%
    mutate(
      stroke = as.factor(stroke),
      ich = as.factor(ich),
      sah = as.factor(sah),
      bt = as.factor(bt),
      sdh = as.factor(sdh),
      tbi = as.factor(tbi)
    ) %>%
    mutate(
      date_consented = as.Date(date_consented, '%d-%b-%Y'),
      date_eeg = as.Date(date_eeg, '%d-%b-%Y'),
      date_jfk_en = as.Date(date_jfk_en, '%d-%b-%Y')
    ) %>%
    mutate(
      date_jfk_dis = as.Date(date_jfk_dis, '%d-%b-%Y'),
      death = as.factor(death),
      code_destination = as.factor(code_destination)
    ) %>%
    mutate(
      destination = as.character(destination),
      outcome_12mo = as.factor(outcome_12mo),
      reason_for_loss = as.character(reason_for_loss),
      notes = as.character(notes)
    ) %>%
    mutate(
      favorable=as.factor(gose >= 4),
      favorable_12mo = as.factor(gose_12mo >=4),
      death_12mo = as.factor(outcome_12mo ==1)
    )
  patient_clinical_data[patient_clinical_data$outcome_12mo == 2,'death_12mo']<- rep(NA,sum(patient_clinical_data$outcome_12mo == 2))
  return(patient_clinical_data)
}