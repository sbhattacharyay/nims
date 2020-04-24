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
      1:length(colnames(patient_clinical_data)),
      colnames(patient_clinical_data),
      sapply(patient_clinical_data, class),
      sapply(patient_clinical_data, typeof),
      colnames(patient_clinical_data),
      clinicalVariableList$variableTitle,
      clinicalVariableList$figName,
      clinicalVariableList$figTitle,
      clinicalVariableList$units
    ) %>%
    rename(
      order = 1,
      variable = 2,
      type = 3,
      class = 4,
      variableSubtitle = 5,
      variableTitle = 6,
      figName = 7,
      figTitle = 8,
      units = 9
    )
  rownames(clinicalVariableList) <- c()
  write.csv(clinicalVariableList,
            '../clinical_data/clinicalVariableList.csv')
  return(clinicalVariableList)
}