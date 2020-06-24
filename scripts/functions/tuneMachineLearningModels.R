tuneMachineLearningModels <- function(seg_window,Iter, DeepIter, classifier_choice, path.D, seed.list, labelsList) {
  
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
  
  for (labelIdx in 1:nrow(labelsList)){
    curr_label <- labelsList$label.name[labelIdx]
    print(paste("Label:",curr_label,"initiated ..."))
    currFeatSet <- readRDS(file.path('../all_motion_feature_data/complete_LOL_set',paste0(seg_window,'_min'),paste0(curr_label,'.rds')))
    mf.only.formula <- as.formula(paste(curr_label,"~",paste(paste0(rep("MF.",r),1:r),collapse = " + ")))
    combined.formula <- as.formula(paste(curr_label,"~",paste(paste0(rep("MF.",r),1:r),collapse = " + "),"+ APACHEIIMortalityRisk"))
    for (foldIdx in 1:length(currFeatSet)){
      print(paste("Fold no.",foldIdx,"initiated ..."))
      curr_train <- currFeatSet[[foldIdx]]$train %>% drop_na(all_of(curr_label))
      curr_test <- currFeatSet[[foldIdx]]$test %>% drop_na(all_of(curr_label))
      for (MLmethod in classifier_choice) {
        print(paste("Training of",MLmethod,"initiated ..."))
        if (MLmethod != "lda") {tune.grid <- get(paste("tune.grid.", MLmethod, sep = ""))}
        set.seed(1876)
        if(MLmethod == "lda") {
          combined.classifier <- train(combined.formula, 
                                       curr_train, 
                                       method = MLmethod,
                                       metric = "ROC", 
                                       trControl = train.control, 
                                       preProcess = c("YeoJohnson", "center", "scale"))
          mf.only.classifier <- train(mf.only.formula, 
                                      curr_train, 
                                      method = MLmethod,
                                      metric = "ROC", 
                                      trControl = train.control, 
                                      preProcess = c("YeoJohnson", "center", "scale"))          
        } else if (MLmethod == "DeepNN") {
          tuneDeepLearningModel(tune.grid = tune.grid,
                                dataset_input = curr_train,
                                seg_window = seg_window,
                                labelName = curr_label,
                                fold = foldIdx,
                                combined = TRUE,
                                seed.list = seed.list,
                                file.path = path.D,
                                iter = DeepIter,
                                verbose = TRUE,
                                no.parallel.cores = no.parallel.cores,
                                tranches.per.core = 2,
                                imputVarTest = FALSE)
          tuneDeepLearningModel(tune.grid = tune.grid,
                                dataset_input = curr_train,
                                seg_window = seg_window,
                                labelName = curr_label,
                                fold = foldIdx,
                                combined = FALSE,
                                seed.list = seed.list,
                                file.path = path.D,
                                iter = DeepIter,
                                verbose = TRUE,
                                no.parallel.cores = no.parallel.cores,
                                tranches.per.core = 2,
                                imputVarTest = FALSE)
        } else {
          combined.classifier <- train(combined.formula, 
                                       curr_train, 
                                       method = MLmethod,
                                       metric = "ROC", 
                                       trControl = train.control, 
                                       preProcess = c("YeoJohnson", "center", "scale"),
                                       tuneGrid = tune.grid)
          mf.only.classifier <- train(mf.only.formula, 
                                      curr_train, 
                                      method = MLmethod,
                                      metric = "ROC", 
                                      trControl = train.control, 
                                      preProcess = c("YeoJohnson", "center", "scale"),
                                      tuneGrid = tune.grid)          
        }
        if(MLmethod != "DeepNN") {
          saveRDSFiles(
            classifier.x = combined.classifier,
            newData = curr_test,
            MLmethod = MLmethod,
            seg_window = seg_window,
            labelName = curr_label,
            fold = foldIdx,
            combined = TRUE,
            iter = Iter,
            path = path.D
          )
          saveRDSFiles(
            classifier.x = mf.only.classifier,
            newData = curr_test,
            MLmethod = MLmethod,
            seg_window = seg_window,
            labelName = curr_label,
            fold = foldIdx,
            combined = FALSE,
            iter = Iter,
            path = path.D
          )
          rm(combined.classifier)
          rm(mf.only.classifier)
        }
        gc()
        print(paste("Training of",MLmethod,"complete."))
      }
      print(paste("Fold no.",foldIdx,"complete."))
    }
    print(paste("Label:",curr_label,"complete."))
  }
}