# Function to save the model, results, and predictions for the re-samples for caret and deep learning models

saveRDSFiles <- function(classifier.x, newData, MLmethod, seg_window, labelName, fold, combined, iter, path, saveClassifier = TRUE) {
  
  curr_test_preds <- classifier.x %>% predict(newData,type = "prob")
  
  if (combined == TRUE) {
    formType <- "combined"
  } else {
    formType <- "MFonly"
  }
  
  for (folder in c("pred", "tuning", "results")) {
    if(!dir.exists(paste(path, folder, MLmethod, sep = "/"))) {
      dir.create(paste(path, folder, MLmethod, sep = "/"), recursive = TRUE)
    }
  }
  
  if(saveClassifier == TRUE) saveRDS(classifier.x, paste(path,"/tuning/",MLmethod,"/ROC.fold",fold,".",MLmethod,".",labelName,".",formType,".iter",iter,".rds", sep = ""))
  saveRDS(classifier.x$pred,   paste(path,"/pred/",MLmethod,"/ROC.fold",fold,".",MLmethod,".",labelName,".",formType,".iter",iter,".rds", sep = ""))
  saveRDS(classifier.x$results,paste(path,"/results/",MLmethod,"/ROC.fold",fold,".",MLmethod,".",labelName,".",formType,".iter",iter,".rds", sep = ""))
}
