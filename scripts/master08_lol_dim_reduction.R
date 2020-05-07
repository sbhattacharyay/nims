#### Master Script 8: Linear Optimal Low Rank (LOL) Dimensionality Reduction ####
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

# Load patient clinical data (sorts by PY numbering and corrects variable types)
patient_clinical_data<-load_patient_clinical_data()

# Load and update clinical variable list
clinicalVariableList<-update_clinicalVariableList()

# Load Motion Features Organized by Time of Day (TOD)
if(!exists("tod_sensors")) {
  tod_sensors<-readMat('../motion_feature_data/bed_corrected_imputed_complete_sensor_data.mat')$bed.corrected.sensors
}

# Load Motion Features Organized by Time from Recording (TFR)
if(!exists("tfr_sensors")) {
  tfr_sensors<-readMat('../tfr_motion_feature_data/bed_corrected_imputed_complete_sensor_data.mat')$bed.corrected.sensors
}

# Organizing Features for both motion feature categorizations
tod_temp_var<-do.call(rbind,tod_sensors)
tfr_temp_var<-do.call(rbind,tfr_sensors)
# ranked_temp_var<-lapply(temp_var, function(x) t(apply(x,1,sort,decreasing = TRUE)))
# fft_temp_var<-lapply(temp_var, function(x) abs(t(mvfft(t(x)))))

tod_motion_features <- get_motion_features(tod_temp_var)
tfr_motion_features <- get_motion_features(tfr_temp_var)

# LOL Projections:
r<-3

Y<-patient_clinical_data$gose
tod_GOSE_LOL<-lol_project_motion_features(tod_motion_features,Y,r)
tfr_GOSE_LOL<-lol_project_motion_features(tfr_motion_features,Y,r)

Y<-patient_clinical_data$gose_12mo
tod_GOSE12mo_LOL<-lol_project_motion_features(tod_motion_features,Y,r)
tfr_GOSE12mo_LOL<-lol_project_motion_features(tfr_motion_features,Y,r)

Y<-patient_clinical_data$favorable
tod_fav_LOL<-lol_project_motion_features(tod_motion_features,Y,r)
tfr_fav_LOL<-lol_project_motion_features(tfr_motion_features,Y,r)

Y<-patient_clinical_data$favorable_12mo
tod_fav12mo_LOL<-lol_project_motion_features(tod_motion_features,Y,r)
tfr_fav12mo_LOL<-lol_project_motion_features(tfr_motion_features,Y,r)

Y<-patient_clinical_data$death
tod_death_LOL<-lol_project_motion_features(tod_motion_features,Y,r)
tfr_death_LOL<-lol_project_motion_features(tfr_motion_features,Y,r)

Y<-patient_clinical_data$death_12mo
tod_death12mo_LOL<-lol_project_motion_features(tod_motion_features,Y,r)
tfr_death12mo_LOL<-lol_project_motion_features(tfr_motion_features,Y,r)

# LOL Projection Visualization:
Y<-patient_clinical_data$gose
viz_lol_2D(tod_GOSE_LOL,Y,"gose_tod",height = 1500, width = 1800, ptSize = 20)
viz_lol_2D(tfr_GOSE_LOL,Y,"gose_tfr",height = 1500, width = 1800, ptSize = 20)

Y<-patient_clinical_data$gose_12mo
viz_lol_2D(tod_GOSE12mo_LOL,Y,"gose_12mo_tod",height = 1500, width = 1800, ptSize = 20)
viz_lol_2D(tfr_GOSE12mo_LOL,Y,"gose_12mo_tfr",height = 1500, width = 1800, ptSize = 20)

Y<-patient_clinical_data$favorable
viz_lol_2D(tod_fav_LOL,Y,"fav_tod",height = 1500, width = 1800, ptSize = 20)
viz_lol_2D(tfr_fav_LOL,Y,"fav_tfr",height = 1500, width = 1800, ptSize = 20)

Y<-patient_clinical_data$favorable_12mo
viz_lol_2D(tod_fav12mo_LOL,Y,"fav_12mo_tod",height = 1500, width = 1800, ptSize = 20)
viz_lol_2D(tfr_fav12mo_LOL,Y,"fav_12mo_tfr",height = 1500, width = 1800, ptSize = 20)

Y<-patient_clinical_data$death
viz_lol_2D(tod_death_LOL,Y,"death_tod",height = 1500, width = 1800, ptSize = 20)
viz_lol_2D(tfr_death_LOL,Y,"death_tfr",height = 1500, width = 1800, ptSize = 20)

Y<-patient_clinical_data$death_12mo
viz_lol_2D(tod_death12mo_LOL,Y,"death_12mo_tod",height = 1500, width = 1800, ptSize = 20)
viz_lol_2D(tfr_death12mo_LOL,Y,"death_12mo_tfr",height = 1500, width = 1800, ptSize = 20)

# Coefficient Analysis:
viz_lol_coeff()