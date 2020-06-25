library(tidyverse)
source('./functions/list_do.R')
source('./functions/list_cbind.R')

buildPredResDF <- function(modelsDF,patient_clinical_data,outerFolds,impNo = 1) {
  dir.create(file.path('..','all_motion_feature_data','ML_predictions',paste0('iter',unique(modelsDF$iter))),recursive=TRUE,showWarnings=FALSE)
  uniqSegWindows <- unique(modelsDF$seg_window)
  for(currSegWindow in uniqSegWindows) {
    sw_modelsDF <- modelsDF %>% filter(seg_window == currSegWindow)
    uniqLabels <- unique(sw_modelsDF$label)
    for (currLabel in uniqLabels){
      avail_Ptidx <- which(!is.na(patient_clinical_data[,currLabel]))
      l_sw_modelsDF <- sw_modelsDF %>% filter(label == currLabel)
      currFeatSet <- readRDS(file.path('../all_motion_feature_data/complete_LOL_set',paste0(currSegWindow,'_min'),paste0(currLabel,'.rds')))
      uniqFormType <- unique(l_sw_modelsDF$formType)
      predsList <- vector(mode = "list", length = length(uniqFormType))
      names(predsList) <- as.character(uniqFormType)
      for (currFormType in uniqFormType){
        ft_l_sw_modelsDF <- l_sw_modelsDF %>% filter(formType == currFormType)
        uniqFolds <- unique(ft_l_sw_modelsDF$fold)
        foldsList <- vector(mode = "list", length = length(uniqFolds))
        names(foldsList) <- paste0(rep("Fold",length(uniqFolds)),1:length(uniqFolds))
        clinOnlyFoldsList <- vector(mode = "list", length = length(uniqFolds))
        names(clinOnlyFoldsList) <- paste0(rep("Fold",length(uniqFolds)),1:length(uniqFolds))
        for (currFold in uniqFolds){
          ground_truth <- currFeatSet[[currFold]]$test[,as.character(currLabel)]
          f_ft_l_sw_modelsDF <- ft_l_sw_modelsDF %>% filter(fold == currFold)
          currTestPredsDF <- list.cbind(lapply(f_ft_l_sw_modelsDF$predPath,FUN = function(x){readRDS(x) %>% dplyr::select(Fav)}))
          names(currTestPredsDF) <- f_ft_l_sw_modelsDF$method
          foldsList[[currFold]] <- cbind(ground_truth,currTestPredsDF)
          clinOnlyFoldsList[[currFold]] <- patient_clinical_data %>% filter(ptIdx %in% avail_Ptidx[outerFolds[[currLabel]][[currFold]]]) %>% dplyr::select(all_of(c('APACHEIIMortalityRisk',as.character(currLabel)))) %>% mutate(probFav = 1 - APACHEIIMortalityRisk) %>% dplyr::select(-APACHEIIMortalityRisk) %>% rename(ground_truth = 1)
        }
        predsList[[currFormType]] <- foldsList
        predsList$CLINonly <- clinOnlyFoldsList
      }
     saveRDS(predsList,file.path('..','all_motion_feature_data','ML_predictions',paste0('iter',unique(modelsDF$iter)),paste0(paste(currLabel,paste0(currSegWindow,'_min'),paste0('imp',impNo),sep = '.'),'.rds'))) 
    }
  }
}