#### Master Script 9: Classification of Mortality Outcomes ####
# Decoding Quantitative Motor Features for Classification and Prediction
# in Severe Acquired Brain Injury
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

library(devtools)
if (!require(lolR))
  install_github('neurodata/lol',
                 build_vignettes = TRUE,
                 force = TRUE)
library(lolR)
library(R.matlab)
library(tidyverse)
library(ggplot2)
library(plotly)
library(naniar)
library(MASS)
library(glmnet)
library(caret)
library(kernlab)
library(rlist)
library(cvAUC)
library(gridExtra)
library(ggplotify)
library(grid)
library(cowplot)

source('./functions/load_patient_clinical_data.R')
source('./functions/update_clinicalVariableList.R')
source('./functions/get_motion_features.R')
source('./functions/lol_project_motion_features.R')
source("./functions/load_tf_patient_covariates.R")
source('./functions/cv_lol_project_motion_features.R')
source('./functions/prepare_training_covariates.R')
source('./functions/prepare_testing_covariates.R')
source('./functions/train_caret_models.R')
source('./functions/predict_caret_models.R')
source('./functions/classification_function_shiny_dis.R')
source('./functions/classification_function_shiny_12mo.R')
source('./functions/get_auc_info.R')
source('./functions/get_auc_plots.R')
source('./functions/get_auc_info_ci.R')
source('./functions/generateRootDir.R')

time_choice<-'tfr'
time_slide<-c(0,8)
classifier_choice<-c("svmRadialWeights","knn","lda","glmnet")
r<-6
mf_choice<-c("band_powerFeats","freq_entropyFeats","freq_pairsFeats","med_freqFeats","smaFeats","waveletsFeats")
sensor_loc<-c("left_ank","left_el","left_wr","right_ank","right_el","right_wr")

preds_dis<-classification_function_shiny_dis(time_choice,time_slide,classifier_choice,r,mf_choice,sensor_loc)
preds_12m<-classification_function_shiny_12mo(time_choice,time_slide,classifier_choice,r,mf_choice,sensor_loc)

auc_dis <- get_auc_info(preds_dis)
auc_dis_ci <- get_auc_info_ci(preds_dis)
auc_12m <- get_auc_info(preds_12m)
auc_12m_ci <- get_auc_info_ci(preds_12m)

auc_dis_plots<-get_auc_plots(auc_dis,auc_dis_ci)
auc_12m_plots<-get_auc_plots(auc_12m,auc_dis_ci)

setwd(generateRootDir())
save_plot(do.call(plot_grid, c(unlist(auc_dis_plots, recursive = F), ncol=2,align='hv')), file="auc_dis.pdf", 
          ncol=2,base_asp = 1.1,base_height=12, base_width=5)
save_plot(do.call(plot_grid, c(unlist(auc_12m_plots, recursive = F), ncol=2,align='hv')), file="auc_12m.pdf", 
          ncol=2,base_asp = 1.1,base_height=12, base_width=5)
setwd('../../scripts')