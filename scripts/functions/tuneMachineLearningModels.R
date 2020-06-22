tuneMachineLearningModels <- function(Iter, DeepIter, classifier_choice, path.D, seed.list, labelsList, fullData.imput.in) {
  
  curr_data_table,outcome,classifier_choice,patient_clinical_data,mf_col_idx,outer_fold_count,inner_fold_count
  
  for (MLmethod in classifier_choice) {
    
    if (MLmethod != "lda") {
      tune.grid <- get(paste("tune.grid.", MLmethod, sep = ""))
    }
    
    for(i in 1:nrow(labelsList)) {
      
      label <- labelsList$label[i]
      temp <- labelsList$temp[i]
      
      set.seed(42)
      if(MLmethod == "glm") {
        classifier <- train(formula.custom, 
                            trainingData, 
                            method = MLmethod,
                            metric = "ROC", 
                            trControl = train.control, 
                            preProcess = c("YeoJohnson", "center", "scale"))
        
      } else if (MLmethod == "APACHE") {
        classifier <- train(formula.APACHE, 
                            trainingData, 
                            method = "glm",
                            metric = "ROC", 
                            trControl = train.control, 
                            preProcess = c("YeoJohnson", "center", "scale"))
        
      } else if (MLmethod == "DeepNN") {
        tuneDeepLearningModel(tune.grid = tune.grid, 
                              dataset_input = trainingData, 
                              day = day, 
                              cumulative = cumul, 
                              seed.list = seed.list, 
                              file.path = path.D,
                              iter = DeepIter, 
                              verbose = TRUE,
                              no.parallel.cores = no.parallel.cores,
                              tranches.per.core = 2,
                              imputVarTest = FALSE,
                              augmentDeaths = augmentDeaths)
        
      } else {
        classifier <- train(formula.simple, 
                            trainingData, 
                            method = MLmethod,
                            metric = "ROC", 
                            trControl = train.control, 
                            preProcess = c("YeoJohnson", "center", "scale"),
                            tuneGrid = tune.grid)
      }
      
      if(MLmethod != "DeepNN") {
        saveRDSFiles(classifier.x = classifier, MLmethod = MLmethod, day = day, iter = Iter, cumulative = cumul, path = path.D)
        rm(classifier)
      }
      gc()
    }
    
  }
  
  
}