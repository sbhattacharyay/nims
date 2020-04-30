library(ggplot)
library(tidyverse)
source("./functions/generateRootDir.R")

plot_AUC_Barplots <- function(DirCSV,newDirName,height = 500, width = 650) {
  
  rootDir<-generateRootDir(file.path("AUC_barplots",newDirName))
  
  curr_data <- read.csv(DirCSV) %>% rename(Model.Covariates = 1)
  
  X<-as.character(curr_data$Model.Covariates[1])
  n_last <- 3                                
  ModelType<-toupper(substr(X, nchar(X) - n_last + 1, nchar(X)))
  
  OutcomeType <- substr(DirCSV, nchar(DirCSV) - 8 + 1, nchar(DirCSV) - 4)
  
  curr_data<-curr_data %>% mutate(Model.Covariates = substr(as.character(Model.Covariates),1,nchar(as.character(Model.Covariates)) - (n_last+1)))
  
  viz<-ggplot(curr_data) +
    geom_bar( aes(x=reorder(Model.Covariates,-tod_GOSE_test_val.mean), y=tod_GOSE_test_val.mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=reorder(Model.Covariates,-tod_GOSE_test_val.mean), ymin=tod_GOSE_test_val.mean-tod_GOSE_test_val.stddev, ymax=tod_GOSE_test_val.mean+tod_GOSE_test_val.stddev), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(x="Predictors used in model",y="AUC (std.dev)",title = paste(ModelType,": (Time of Day MF, GOSE-Trained, Predicting",OutcomeType,"Fav)"))
  fileName <- file.path(rootDir,"tod_GOSE_Fav.png")
  png(file = fileName,
      width = width,
      height = height)
  print(viz)
  dev.off()
  
  viz<-ggplot(curr_data) +
    geom_bar( aes(x=reorder(Model.Covariates,-tod_fav_test_val.mean), y=tod_fav_test_val.mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=reorder(Model.Covariates,-tod_fav_test_val.mean), ymin=tod_fav_test_val.mean-tod_fav_test_val.stddev, ymax=tod_fav_test_val.mean+tod_fav_test_val.stddev), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(x="Predictors used in model",y="AUC (std.dev)",title = paste(ModelType,": (Time of Day MF, Fav-Trained, Predicting",OutcomeType,"Fav)"))
  fileName <- file.path(rootDir,"tod_Fav_Fav.png")
  png(file = fileName,
      width = width,
      height = height)
  print(viz)
  dev.off()
  
  viz<-ggplot(curr_data) +
    geom_bar( aes(x=reorder(Model.Covariates,-tod_death_test_val.mean), y=tod_death_test_val.mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=reorder(Model.Covariates,-tod_death_test_val.mean), ymin=tod_death_test_val.mean-tod_death_test_val.stddev, ymax=tod_death_test_val.mean+tod_death_test_val.stddev), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(x="Predictors used in model",y="AUC (std.dev)",title = paste(ModelType,": (Time of Day MF, Death-Trained, Predicting",OutcomeType,"Death)"))
  fileName <- file.path(rootDir,"tod_Death_Death.png")
  png(file = fileName,
      width = width,
      height = height)
  print(viz)
  dev.off()
  
  viz<-ggplot(curr_data) +
    geom_bar( aes(x=reorder(Model.Covariates,-tfr_GOSE_test_val.mean), y=tfr_GOSE_test_val.mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=reorder(Model.Covariates,-tfr_GOSE_test_val.mean), ymin=tfr_GOSE_test_val.mean-tfr_GOSE_test_val.stddev, ymax=tfr_GOSE_test_val.mean+tfr_GOSE_test_val.stddev), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(x="Predictors used in model",y="AUC (std.dev)",title = paste(ModelType,": (Time from Recording MF, GOSE-Trained, Predicting",OutcomeType,"Fav)"))
  fileName <- file.path(rootDir,"tfr_GOSE_Fav.png")
  png(file = fileName,
      width = width,
      height = height)
  print(viz)
  dev.off()
  
  viz<-ggplot(curr_data) +
    geom_bar( aes(x=reorder(Model.Covariates,-tfr_fav_test_val.mean), y=tfr_fav_test_val.mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=reorder(Model.Covariates,-tfr_fav_test_val.mean), ymin=tfr_fav_test_val.mean-tfr_fav_test_val.stddev, ymax=tfr_fav_test_val.mean+tfr_fav_test_val.stddev), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(x="Predictors used in model",y="AUC (std.dev)",title = paste(ModelType,": (Time from Recording MF, Fav-Trained, Predicting",OutcomeType,"Fav)"))
  fileName <- file.path(rootDir,"tfr_Fav_Fav.png")
  png(file = fileName,
      width = width,
      height = height)
  print(viz)
  dev.off()
  
  viz<-ggplot(curr_data) +
    geom_bar( aes(x=reorder(Model.Covariates,-tfr_death_test_val.mean), y=tfr_death_test_val.mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=reorder(Model.Covariates,-tfr_death_test_val.mean), ymin=tfr_death_test_val.mean-tfr_death_test_val.stddev, ymax=tfr_death_test_val.mean+tfr_death_test_val.stddev), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(x="Predictors used in model",y="AUC (std.dev)",title = paste(ModelType,": (Time from Recording MF, Death-Trained, Predicting",OutcomeType,"Death)"))
  fileName <- file.path(rootDir,"tfr_Death_Death.png")
  png(file = fileName,
      width = width,
      height = height)
  print(viz)
  dev.off()
}

for (i in list.files(path="../plots/2020-04-30/",pattern="^test.*csv")){
  curr_Path <- file.path("../plots/2020-04-30/",i)
  currName <- substr(i, 1, nchar(i) - 4)
  plot_AUC_Barplots(curr_Path,currName)
}