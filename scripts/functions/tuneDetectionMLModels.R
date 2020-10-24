tuneDetectionMLModels <- function(curr_SMOTE, inner_fold_count, classifier_choice, save.path, seed.list, ordinalLabels) {

  # Rename variables of dataframes to fit formula
  names(curr_SMOTE)[1:(ncol(curr_SMOTE)-1)] <- paste(paste0(rep("MF.",(ncol(curr_SMOTE)-1)),1:(ncol(curr_SMOTE)-1)))
  curr_SMOTE <- curr_SMOTE[1:(ncol(curr_SMOTE)-1)]
  names(ordinalLabels) <- paste(paste0(rep("label.",ncol(ordinalLabels)),1:ncol(ordinalLabels)))
  curr_SMOTE <- cbind(curr_SMOTE,ordinalLabels)

  labelCols <- which(grepl("label",names(curr_SMOTE)))
  
  for (labelCol in labelCols){
    levels(curr_SMOTE[,labelCol])<-c("less","greater")
  }
  
  for (MLmethod in classifier_choice) {
    print(paste("Training of",MLmethod,"initiated ..."))
    
    currMethodMdls <- vector(mode = "list", length = ncol(ordinalLabels))
    
    for (labelIdx in 1:ncol(ordinalLabels)) {
      gc()
      print(paste("Ordinal label:",labelIdx,"initiated ..."))
      curr_label <- names(ordinalLabels)[labelIdx]
      
      ML_formula <- as.formula(paste(curr_label,"~",paste(paste0(rep("MF.",ncol(curr_SMOTE)-ncol(ordinalLabels)),1:(ncol(curr_SMOTE)-ncol(ordinalLabels))),collapse = " + ")))
      
      label_balance <- prop.table(table(ordinalLabels[,labelIdx]))
      lower_class <- names(label_balance)[which.min(label_balance)]
      model_weights <- ifelse(ordinalLabels[,labelIdx] == lower_class,
                              label_balance[[1]],
                              label_balance[[2]])
      
      train.control <- caret::trainControl(method="repeatedcv",
                                    number=inner_fold_count,
                                    repeats=1,
                                    summaryFunction = twoClassSummary, 
                                    classProbs = TRUE,
                                    seeds = seed.list,
                                    allowParallel = TRUE)
      
      if (MLmethod != "lda") {tune.grid <- get(paste("tune.grid.", MLmethod, sep = ""))}
      
      if(MLmethod == "lda") {
        currMethodMdls[[labelIdx]] <- train(ML_formula, 
                                            curr_SMOTE, 
                                            method = MLmethod,
                                            metric = "ROC", 
                                            trControl = train.control, 
                                            weights = model_weights,
                                            preProcess = c("YeoJohnson", "center", "scale"))
      } else {
        currMethodMdls[[labelIdx]] <- train(ML_formula, 
                                            curr_SMOTE, 
                                            method = MLmethod,
                                            metric = "ROC", 
                                            trControl = train.control, 
                                            weights = model_weights,
                                            preProcess = c("YeoJohnson", "center", "scale"),
                                            tuneGrid = tune.grid)
      }
      print(paste("Ordinal label:",labelIdx,"completed."))
    }
    saveRDS(currMethodMdls,file = file.path(save.path,paste0(MLmethod,".rds")))
    gc()
    print(paste("Training of",MLmethod,"complete."))
  }
}