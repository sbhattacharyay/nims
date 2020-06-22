motion_feature_preparation <-function(seg_window,mf_choice,sr_choice,labels.temp){
  
  
  seg_dir <- list.dirs('../all_motion_feature_data/')[str_detect(list.dirs('../all_motion_feature_data/'),pattern=paste0('final_processed_features_',as.character(seg_window)))]
  seg_files <- list.files(seg_dir,pattern = "*.csv",full.names = TRUE)
  first_imp <-
    read_csv(seg_files[1]) %>% dplyr::select(-1, -2) %>% 
    filter(featureType %in% mf_choice) %>% 
    dplyr::select(all_of(c(sr_choice, "ptIdx", "segTimePoint", "featureType", "segIdx")))
  
  
  print('Casting data frame ...')
  curr_data_table <- dcast(melt(imp,id.vars = c("ptIdx","featureType","segIdx","segTimePoint")),ptIdx+segIdx~featureType+variable+segTimePoint)
  print('Casting complete')
  
  labelIdx <- which(colnames(curr_data_table)==outcome)
  amrIIIdx <- which(colnames(curr_data_table)=="APACHEIIMortalityRisk")
  
  avail_Ptidx <- !is.na(patient_clinical_data[,outcome])
  avail_seq <- seq(length(avail_Ptidx))[avail_Ptidx]
  
  # Filter out missing patients:
  filt_curr_data_table <- curr_data_table %>% filter(ptIdx %in% avail_seq)
  set.seed(123)
  outcomes_avail<-patient_clinical_data[avail_seq,outcome]
  
  # Nested stratified cross-validation: Create 5 folds for outer cross-validation (training vs. testing)
  outer_cvIdx <- createFolds(outcomes_avail,outer_fold_count,list = TRUE);
  
  for(i in 1:outer_fold_count){
    urr_test_data <- filt_curr_data_table %>% filter(ptIdx %in% outer_cvIdx[[i]])
    curr_train_data <- filt_curr_data_table %>% filter(!ptIdx %in% outer_cvIdx[[i]])
    
    set.seed(1876)
    #down_train <- SMOTE(as.formula(paste(outcome,"~.")), data  = curr_train_data[,c(mf_col_idx,amrIIIdx,labelIdx)])
    down_train <- downSample(x = curr_train_data[,!(names(curr_train_data) == outcome)],
                             y = curr_train_data[,outcome], yname = outcome)
    
    n_train <- nrow(down_train)
    n_mfeat <- length(mf_col_idx)
    
    if (n_mfeat > n_train/5){
      print("LOL initatied ...")
      r <- 125
      LOL <-  lol.project.lol(as.matrix(down_train[,(mf_col_idx-2)]), down_train[,outcome], 125)
      reduced_MF_DF <- cbind(as.data.frame(LOL$Xr[,1:r]),down_train[,outcome])
      print("LOL complete")
    } else {
      r <- n_mfeat
      reduced_MF_DF <- cbind(as.data.frame(as.matrix(down_train[,(mf_col_idx-2)])),down_train[,outcome])
    }
    names(reduced_MF_DF)<-c(paste0(rep("MF.",r),1:r),outcome)
    APACHEIIMortalityRisk <- down_train[,c("APACHEIIMortalityRisk")]
    combined_DF <- cbind(reduced_MF_DF,APACHEIIMortalityRisk)
    
  }
  
}