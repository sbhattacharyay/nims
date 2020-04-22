#### Clean Clinical Table Data  ####

library('tidyverse')
library('naniar')
library('writexl')

na_strings <- c("NA", "N A", "N / A", "N/A", "N/ A", "Not Available", "NOt available","","NaT","NaN")

patientData<-read.csv('../clinical_data/original_patient_data.csv')%>%
  replace_with_na_all(condition = ~.x %in% na_strings)
clinicalVariableList<-read_excel('../clinical_data/clinicalVariableList.xlsx')
colnames(patientData)<-clinicalVariableList$variable
clinicalVariableList<-mutate(clinicalVariableList,class=sapply(patientData, class))

write.csv(patientData,'../clinical_data/clean_patient_data.csv')