motion_feature_preparation <-function(patient_clinical_data,seg_window,mf_choice,sr_choice,outerFolds,impNo = 1,saveDir,DSseed = 2020){
  
  if (seg_window == 5){
    r <- 1000
  } else if (seg_window == 10) {
    r <- 500
  } else if (seg_window == 30) {
    r <- 175
  } else if (seg_window == 60) {
    r <- 80
  } else if (seg_window == 180) {
    r <- 25
  } else {
    stop("Incompatible segment window. Segment window must be 5, 10, 30, 60, or 180 min")
  }
  dir.create(saveDir,recursive = T,showWarnings = F)
  seg_files <-
    list.files(
      list.dirs('../all_motion_feature_data/')[str_detect(
        list.dirs('../all_motion_feature_data/'),
        pattern = paste0('final_processed_features_', as.character(seg_window))
      )],
      pattern = "*.csv",
      full.names = TRUE
    )
  imp <-
    read_csv(seg_files[impNo]) %>% dplyr::select(-1, -2) %>% 
    filter(featureType %in% mf_choice) %>% 
    dplyr::select(all_of(c(sr_choice, "ptIdx", "segTimePoint", "featureType", "segIdx")))
  mf_col_idx <- 3:(2+seg_window*12*length(sr_choice)*length(mf_choice))
  
  print(paste("Casting imputation no.",impNo,"data into predictor form ..."))
  cast_imp <- dcast(melt(imp,id.vars = c("ptIdx","featureType","segIdx","segTimePoint")),ptIdx+segIdx~featureType+variable+segTimePoint)
  print(paste("Casting complete"))
  
  for (l in 1:length(outerFolds)){
    curr_label <- names(outerFolds)[l]
    curr_outerFold <- outerFolds[[l]]
    avail_Ptidx <- which(!is.na(patient_clinical_data[,curr_label]))
    filt_cast_imp <- cast_imp %>% filter(ptIdx %in% avail_Ptidx) %>% inner_join(.,patient_clinical_data,by="ptIdx")
    labelIdx <- which(names(filt_cast_imp) == curr_label)
    apacheIIriskIdx <- which(names(filt_cast_imp) == "APACHEIIMortalityRisk")
    feature_list <- vector(mode = "list", length = length(curr_outerFold))
    for (j in 1:length(curr_outerFold)){
      curr_test <- filt_cast_imp %>% filter(ptIdx %in% avail_Ptidx[curr_outerFold[[j]]])
      curr_train <- filt_cast_imp  %>% filter(!ptIdx %in% avail_Ptidx[curr_outerFold[[j]]])
      set.seed(seed = DSseed)
      down_train <- downSample(x = curr_train[,!(names(curr_train) == curr_label)],
                               y = curr_train[,curr_label], yname = curr_label)
      print(paste("LOL for",curr_label,"fold:",j,"initatied ..."))
      curr_LOL <-  lol.project.lol(as.matrix(down_train[,mf_col_idx]), down_train[,curr_label], r)
      train_DF <- cbind(down_train$ptIdx,down_train$segIdx,as.data.frame(curr_LOL$Xr[,1:r]),down_train[,c(apacheIIriskIdx,labelIdx)])
      names(train_DF)<-c("ptIdx","segIdx",paste0(rep("MF.",r),1:r),"APACHEIIMortalityRisk",curr_label)
      test_DF <- cbind(curr_test$ptIdx,curr_test$segIdx,as.data.frame(as.matrix(curr_test[,mf_col_idx])%*%curr_LOL$A[,1:r]),curr_test[,c(apacheIIriskIdx,labelIdx)])
      names(test_DF)<-c("ptIdx","segIdx",paste0(rep("MF.",r),1:r),"APACHEIIMortalityRisk",curr_label)
      print(paste("LOL for",curr_label,"fold:",j,"complete."))
      tempList <- list(train_DF,test_DF)
      names(tempList) <- c("train","test")
      feature_list[[j]] <- tempList
    }
    print(paste('Saving',curr_label,'initiated ...'))
    saveRDS(feature_list, file = file.path(saveDir,paste0(curr_label,'.rds')))
    print(paste('Saving',curr_label,'complete.'))
  }
}