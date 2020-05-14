library(devtools)
if (!require(lolR))
  install_github('neurodata/lol',
                 build_vignettes = TRUE,
                 force = TRUE)
library(lolR)
library(R.matlab)
library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(tibble)
library(stringr)
library(forcats)
library(readxl)
library(plotly)
library(naniar)
library(MASS)
library(glmnet)
library(caret)
library(kernlab)
library(rlist)
library(nnet)
library(e1071)
library(randomForest)
library(foreach)
library(gridExtra)
library(ggplotify)
library(precrec)
library(cvAUC)
library(grid)
library(cowplot)
library(shinycssloaders)
source('./functions/load_patient_clinical_data.R')
source('./functions/update_clinicalVariableList.R')
source('./functions/get_motion_features.R')
source('./functions/lol_project_motion_features.R')
source("./functions/strat_cross_val_splits.R")
source("./functions/load_tf_patient_covariates.R")
source('./functions/cv_lol_project_motion_features.R')
source('./functions/prepare_training_covariates.R')
source('./functions/prepare_testing_covariates.R')
source('./functions/classification_function_shiny_dis.R')
source('./functions/classification_function_shiny_12mo.R')
source('./functions/train_caret_models.R')
source('./functions/predict_caret_models.R')
source('./functions/getAUC.R')

optim_mf <- function(classification_function){
  use_fun <- match.fun(classification_function)
  motion_feats <- c("band_powerFeats","freq_entropyFeats","freq_pairsFeats", 
                    "med_freqFeats","smaFeats","waveletsFeats")
  motion_feats <- foreach(n=1:length(motion_feats),.combine = c) %do% combn(motion_feats,n,simplify=FALSE)
  
  svm_preds<-lapply(motion_feats, function(x) use_fun("tfr",c(0,8),'svmRadialWeights',
                                                       3,x,c('APACHE'),c("left_ank","left_el","left_wr", 
                                                                          "right_ank","right_el","right_wr")))
  names(svm_preds)<-sapply(motion_feats, function(x) paste(x, sep=" ", collapse="."))
  svm_auc<-lapply(svm_preds, function(x)getAUC(x))
  svm_auc_mean<-svm_auc %>% map(1) %>% as.data.frame()
  svm_auc_sd<-svm_auc %>% map(2) %>% as.data.frame()
  
  knn_preds<-lapply(motion_feats, function(x) use_fun("tfr",c(0,8),'knn',
                                                     3,x,c('APACHE'),c("left_ank","left_el","left_wr", 
                                                                      "right_ank","right_el","right_wr")))
  names(knn_preds)<-sapply(motion_feats, function(x) paste(x, sep=" ", collapse="."))
  knn_auc<-lapply(knn_preds, function(x)getAUC(x))
  knn_auc_mean<-knn_auc %>% map(1) %>% as.data.frame()
  knn_auc_sd<-knn_auc %>% map(2) %>% as.data.frame()
  
  glmnet_preds<-lapply(motion_feats, function(x) use_fun("tfr",c(0,8),'glmnet',
                                                         3,x,c('APACHE'),c("left_ank","left_el","left_wr", 
                                                                            "right_ank","right_el","right_wr")))
  names(glmnet_preds)<-sapply(motion_feats, function(x) paste(x, sep=" ", collapse="."))
  glmnet_auc<-sapply(glmnet_preds, function(x)getAUC(x))
  glmnet_auc_mean<-glmnet_auc %>% map(1) %>% as.data.frame()
  glmnet_auc_sd<-glmnet_auc %>% map(2) %>% as.data.frame()
  
  lda_preds<-lapply(motion_feats, function(x) use_fun("tfr",c(0,8),'lda',
                                                       3,x,c('APACHE'),c("left_ank","left_el","left_wr", 
                                                                          "right_ank","right_el","right_wr")))
  names(lda_preds)<-sapply(motion_feats, function(x) paste(x, sep=" ", collapse="."))
  lda_auc<-lapply(lda_preds, function(x)getAUC(x))
  lda_auc_mean<-lda_auc %>% map(1) %>% as.data.frame()
  lda_auc_sd<-lda_auc %>% map(2) %>% as.data.frame()  
  
  avNNet_preds<-lapply(motion_feats, function(x) use_fun("tfr",c(0,8),'avNNet',
                                                          3,x,c('APACHE'),c("left_ank","left_el","left_wr", 
                                                                            "right_ank","right_el","right_wr")))
  names(avNNet_preds)<-sapply(motion_feats, function(x) paste(x, sep=" ", collapse="."))
  avNNet_auc<-lapply(avNNet_preds, function(x)getAUC(x))
  avNNet_auc_mean<-lda_auc %>% map(1) %>% as.data.frame()
  avNNet_auc_sd<-lda_auc %>% map(2) %>% as.data.frame()  
  
  mf_auc_means <- list(svm=svm_auc_mean,knn=knn_auc_mean,glmnet=glmnet_auc_mean,
                    lda=lda_auc_mean,avNNet=avNNet_auc_mean)
  mf_auc_stds <- list(svm=svm_auc_sd,knn=knn_auc_sd,glmnet=glmnet_auc_sd,
                       lda=lda_auc_sd,avNNet=avNNet_auc_sd)  
  return(list(means=mf_auc_means, stds=mf_auc_stds))
}

mf_dis<-optim_mf("classification_function_shiny_dis")
save(mf_dis,file="./vignettes/motion_feat_gridsearch_dis.rdata")
mf_12mo<-optim_mf("classification_function_shiny_12mo")
save(mf_12mo,file="./vignettes/motion_feat_gridsearch_12mo.rdata")