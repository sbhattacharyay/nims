update_clinicalVariableList <- function(directory){
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
    read.csv(directory, na.strings = na_strings)
  clinicalVariableList <-
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
  rownames(clinicalVariableList) <- c()
  write.csv(clinicalVariableList,
            '../clinical_data/clinicalVariableList.csv')
  return(clinicalVariableList)
}