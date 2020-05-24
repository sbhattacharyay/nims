update_clinicalVariableList <- function(){
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
  clinicalVariableList <-
    read.csv('../clinical_data/clinicalVariableList.csv', na.strings = na_strings)
  clinicalVariableList <-
    data.frame(
      colnames(patient_clinical_data),
      sapply(patient_clinical_data, class),
      sapply(patient_clinical_data, typeof),
      colnames(patient_clinical_data)
    ) %>%
    rename(
      variable = 1,
      type = 3,
      class = 2,
      variableSubtitle = 4
    )
  rownames(clinicalVariableList) <- c()
  write.csv(clinicalVariableList,
            '../clinical_data/clinicalVariableList.csv')
  return(clinicalVariableList)
}