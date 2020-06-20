#### Master Script 6: Temporal Segmentation and Predictor Preparation ####
# Decoding Quantitative Motor Features for Classification and Prediction
# in Severe Acquired Brain Injury
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

library(tidyverse)

if (!exists("compiledImputations")) {
  compiledImputations <- readRDS('../all_motion_feature_data/bc_i_c_dataset.rds')
}

for (seg_window in c(5,10,30,60,180)) {
  print(paste("Segment Window Size",seg_window,"min Started"))
  
  new_dir_name <- file.path("..","all_motion_feature_data",paste0("final_processed_features_",seg_window,"min"))
  
  n <- length(unique(compiledImputations[[1]]$ptIdx))
  
  # Algorithm strategy: assign new label to data that corresponds to chunk Index and provide a timeIdx within chunk
  dir.create(new_dir_name,showWarnings = FALSE)
  for (l in 1:length(compiledImputations)){
    print(paste("Imputation No.",l,"Started"))
    currImp <- compiledImputations[[l]]
    newImp <- as.data.frame(matrix(nrow = 0,ncol = 12))
    for (i in 1:n) {
      curr_patient <- currImp %>% filter(ptIdx == i)
      currPtTimePoints <- unique(curr_patient$timeCount)
      no_splits <- floor(length(currPtTimePoints)/(seg_window*12))
      truncTime <- currPtTimePoints[1:((seg_window*12)*no_splits)]
      indexSplits <- split(truncTime, ceiling(seq_along(truncTime)/(seg_window*12)))
      indexDF <- as.data.frame(indexSplits) %>% gather() %>% rename(segIdx=1,timeCount=2) %>% mutate(segIdx = gsub("[^0-9.-]", "", segIdx)) %>% mutate(segIdx=as.numeric(segIdx))
      indexDF$segTimePoint <- rep(1:length(indexSplits[[1]]),length(indexSplits))
      curr_patient <- curr_patient %>% inner_join(.,indexDF,by="timeCount") %>% dplyr::select(-timeCount)
      newImp <- rbind(newImp,curr_patient)
      print(paste("Patient No.",i,"Complete"))
    }
    imp_file_path <- file.path(new_dir_name,paste0("imp",l,".csv"))
    write.csv(newImp,imp_file_path)
    print(paste("Imputation No.",l,"Finished and Saved"))
  }
  #rm(compiledImputations)
  print(paste("Segment Window Size",seg_window,"min Complete"))
}