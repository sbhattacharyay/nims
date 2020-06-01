load_tf_patient_covariates <- function(directory) {
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
  
  tf_patient_covariates <-
    read.csv(directory, na.strings = na_strings) %>%
    arrange(PY_pNum) %>%
    mutate_at(4:10, as.factor) %>%
    mutate(SDH.TBI = factor(as.integer(TBI == 1 | SDH == 1))) %>%
    mutate_at(14:17, as.factor)
  return(tf_patient_covariates)
}