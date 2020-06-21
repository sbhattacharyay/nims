source('./functions/list_do.R')
source('./functions/list_cbind.R')

seg_classification_function <-function(curr_data_table,outcome,classifier_choice,patient_clinical_data,mf_col_idx,outer_fold_count,inner_fold_count){
  
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
  
  combined_train_predictions <- vector(mode = "list", length = outer_fold_count)
  clin_only_train_predictions  <- vector(mode = "list", length = outer_fold_count)
  mf_only_train_predictions  <- vector(mode = "list", length = outer_fold_count)
  
  combined_test_predictions <- vector(mode = "list", length = outer_fold_count)
  clin_only_test_predictions  <- vector(mode = "list", length = outer_fold_count)
  mf_only_test_predictions  <- vector(mode = "list", length = outer_fold_count)
  
  for(i in 1:outer_fold_count){
    curr_test_data <- filt_curr_data_table %>% filter(ptIdx %in% outer_cvIdx[[i]])
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
    
    combined_formula <- as.formula(paste(outcome,"~."))
    mf_only_formula <- as.formula(paste(outcome,"~",paste(paste0(rep("MF.",r),1:r),collapse = " + ")))
    clin_only_formula <- as.formula(paste(outcome,"~ APACHEIIMortalityRisk"))
    
    # label_balance <- prop.table(table(combined_DF[,outcome]))
    # lower_class <- names(label_balance)[which.min(label_balance)]
    # model_weights <- ifelse(down_train[,outcome] == lower_class,
    #                         label_balance[[1]],
    #                         label_balance[[2]])
    tc <- trainControl(method="repeatedcv",
                       number=inner_fold_count,
                       repeats=1,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       allowParallel = TRUE)
                       #verboseIter = FALSE,
                       #savePredictions = "all",
                       #returnResamp = "all")
    print(paste("Training combined models started for fold",i))
    set.seed(2020)
    combined_models <- lapply(classifier_choice,function(x) train(
      combined_formula, 
      data = combined_DF, 
      method = x,
      metric = "ROC",
      trControl = tc,
      # weights = model_weights,
      preProcess = c("YeoJohnson","center","scale"),
      tuneLength = 3
    ))
    names(combined_models)<-classifier_choice
    print(paste("Training combined models complete for fold",i))
    
    print(paste("Training motion feature only models started for fold",i))
    set.seed(2020)
    mf_only_models <- lapply(classifier_choice,function(x) train(
      mf_only_formula, 
      data = combined_DF, 
      method = x,
      metric = "ROC",
      trControl = tc,
      # weights = model_weights,
      preProcess = c("YeoJohnson","center","scale"),
      tuneLength = 3
    ))
    names(mf_only_models)<-classifier_choice
    print(paste("Training motion feature only complete for fold",i))
    
    print(paste("Training clinical feature only models started for fold",i))
    train_patient_table <- patient_clinical_data[-outer_cvIdx[[i]],]
    train_patient_table[,outcome] <- (train_patient_table[,outcome])
    
    label_balance <- prop.table(table(train_patient_table[,outcome]))
    lower_class <- names(label_balance)[which.min(label_balance)]
    model_weights <- ifelse(train_patient_table[,outcome] == lower_class,
                            label_balance[[1]],
                            label_balance[[2]])
    tc <- trainControl(method="repeatedcv",
                       number=4,
                       repeats=1,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary,
                       allowParallel = TRUE)
                       #verboseIter = FALSE,
                       #savePredictions = "all",
                       #returnResamp = "all")

    set.seed(2020)
    if ('glmnet_h2o' %in% classifier_choice){temp_classifier_choice<-replace(classifier_choice, classifier_choice == 'glmnet_h2o', 'glm')}
    clin_only_models <- lapply(temp_classifier_choice,function(x) train(
      clin_only_formula, 
      data = train_patient_table, 
      method = x,
      metric = "ROC",
      trControl = tc,
      weights = model_weights,
      preProcess = c("YeoJohnson","center","scale"),
      tuneLength = 3
    ))
    names(clin_only_models)<-temp_classifier_choice
    print(paste("Training clinical feature only models complete for fold",i))
    
    # Self-validation predictions:
    combined_train_preds <- cbind(down_train[,outcome],as.data.frame(list.cbind(lapply(combined_models, function(x) x%>%predict(combined_DF,type = "prob")%>%dplyr::select(2)))))    
    names(combined_train_preds)<-c('ground_truth',names(combined_models))
    
    mf_only_train_preds <- cbind(down_train[,outcome],as.data.frame(list.cbind(lapply(mf_only_models, function(x) x%>%predict(combined_DF,type = "prob")%>%dplyr::select(2)))))    
    names(mf_only_train_preds)<-c('ground_truth',names(mf_only_models))
    
    clin_only_train_preds <- cbind(train_patient_table[,outcome],as.data.frame(list.cbind(lapply(clin_only_models, function(x) x%>%predict(train_patient_table,type = "prob")%>%dplyr::select(2)))))    
    names(clin_only_train_preds)<-c('ground_truth',names(clin_only_models))
    
    combined_train_predictions[[i]] <- combined_train_preds
    clin_only_train_predictions[[i]] <- clin_only_train_preds
    mf_only_train_predictions[[i]] <- mf_only_train_preds
    
    # Test set predictions:
    if (n_mfeat > n_train/2){
      test_MF_DF <- cbind(as.data.frame(as.matrix(curr_test_data[,mf_col_idx])%*%LOL$A[,1:r]),curr_test_data[,outcome])
    } else {
      test_MF_DF <- cbind(as.data.frame(curr_test_data[,mf_col_idx]),curr_test_data[,outcome])
    }
    names(test_MF_DF)<-c(paste0(rep("MF.",r),1:r),outcome)
    APACHEIIMortalityRisk <- curr_test_data[,c("APACHEIIMortalityRisk")]
    combined_test_DF <- cbind(test_MF_DF,APACHEIIMortalityRisk)
    test_patient_table <- patient_clinical_data[outer_cvIdx[[i]],]

    combined_test_preds <- cbind(curr_test_data[,outcome],as.data.frame(list.cbind(lapply(combined_models, function(x) x%>%predict(combined_test_DF,type = "prob")%>%dplyr::select(2)))))    
    names(combined_test_preds)<-c('ground_truth',names(combined_models))
    
    mf_only_test_preds <- cbind(curr_test_data[,outcome],as.data.frame(list.cbind(lapply(mf_only_models, function(x) x%>%predict(combined_test_DF,type = "prob")%>%dplyr::select(2)))))    
    names(mf_only_test_preds)<-c('ground_truth',names(mf_only_models))
    
    clin_only_test_preds <- cbind(test_patient_table[,outcome],as.data.frame(list.cbind(lapply(clin_only_models, function(x) x%>%predict(test_patient_table,type = "prob")%>%dplyr::select(2)))))    
    names(clin_only_test_preds)<-c('ground_truth',names(clin_only_models))

    combined_test_predictions[[i]] <- combined_test_preds
    clin_only_test_predictions[[i]] <- clin_only_test_preds
    mf_only_test_predictions[[i]] <- mf_only_test_preds
    print(paste("Fold",i,"Complete"))
  }
  
  outList1 <-
    list(
      combined_test_predictions,
      clin_only_test_predictions,
      mf_only_test_predictions
    )
  names(outList1) <-
    c(
      "combined_test_predictions",
      "clin_only_test_predictions",
      "mf_only_test_predictions"
    )
  outList2 <-
    list(
      combined_train_predictions,
      clin_only_train_predictions,
      mf_only_train_predictions
    )
  names(outList2) <-
    c(
      "combined_train_predictions",
      "clin_only_train_predictions",
      "mf_only_train_predictions"
    )
  
  finalOutList <- list(outList1,outList2)
  return(finalOutList)
}