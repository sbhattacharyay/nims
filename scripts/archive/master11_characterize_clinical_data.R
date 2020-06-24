#### Master Script 10: Characterize Clinical Data ####
# Decoding Quantitative Motor Features for Classification and Prediction
# in Severe Acquired Brain Injury
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

library(devtools)
if (!require(lolR)) install_github('neurodata/lol', build_vignettes=TRUE, force=TRUE)
library(lolR)
library(R.matlab)
library(tidyverse)
library(ggplot2)
library(plotly)
library(naniar)
library(MASS)

source('./functions/load_patient_clinical_data.R')
source('./functions/update_clinicalVariableList.R')
source('./functions/get_motion_features.R')
source('./functions/lol_project_motion_features.R')
source('./functions/viz_lol_2D.R')
source("./functions/generateRootDir.R")
source("./functions/plot_descriptive_figs.R")

# Load patient clinical data (sorts by PY numbering and corrects variable types)
patient_clinical_data<-load_patient_clinical_data()

# Load and update clinical variable list
clinicalVariableList<-update_clinicalVariableList()

plotDescriptiveFigs(patient_clinical_data,clinicalVariableList)