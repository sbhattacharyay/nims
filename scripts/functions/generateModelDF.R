# Function to retrieve the predictions and results from machine learning model output files
generateModelDF <- function(path.D = "../all_motion_feature_data/ML_results",outer_fold_count,iter=1) {
  # generate list of possible classifiers
  combinations <- expand.grid(MLmethod = c("adaboost", "avNNet","DeepNN","glmnet", "parRF", "svmRadialWeights","lda"),
                              seg_window = c(5,10,30,60,180),
                              label.name = c("fav_mort_dis","fav_GOSE_dis","fav_mRS_dis","fav_mort_12m","fav_GOSE_12m","fav_mRS_12m"),
                              combined = c("combined", "MFonly"),
                              fold = 1:outer_fold_count)
  
  # output model df is currently empty
  out.model.df.length <- 0
  out.model.df <- data.frame()
  
  # try each classifier and if it exists, add it to a list
  for (i in 1:nrow(combinations)) {
    # select classifier to try from combinations
    combTry <- combinations[i,]
    
    #get the right path
    path <- path.D
    
    # try to load the results (smallest of the three options --> fastest)
    fileExists <- file.exists(file.path(path,paste0(combTry$seg_window,"_min"),"results",combTry$MLmethod,paste("ROC",paste0(combTry$seg_window,"_min"),paste0("fold",combTry$fold),combTry$MLmethod,combTry$label.name,combTry$combined,paste0("iter",iter),"rds",sep = ".")))
    
    out.model.df.row <- data.frame(method = combTry$MLmethod,
                                   seg_window = combTry$seg_window,
                                   label = combTry$label.name,
                                   fold = combTry$fold,
                                   classifierPath = file.path(path,paste0(combTry$seg_window,"_min"),"tuning",combTry$MLmethod,paste("ROC",paste0(combTry$seg_window,"_min"),paste0("fold",combTry$fold),combTry$MLmethod,combTry$label.name,combTry$combined,paste0("iter",iter),"rds",sep = ".")), 
                                   predPath = file.path(path,paste0(combTry$seg_window,"_min"),"pred",combTry$MLmethod,paste("ROC",paste0(combTry$seg_window,"_min"),paste0("fold",combTry$fold),combTry$MLmethod,combTry$label.name,combTry$combined,paste0("iter",iter),"rds",sep = ".")), 
                                   resultsPath = file.path(path,paste0(combTry$seg_window,"_min"),"results",combTry$MLmethod,paste("ROC",paste0(combTry$seg_window,"_min"),paste0("fold",combTry$fold),combTry$MLmethod,combTry$label.name,combTry$combined,paste0("iter",iter),"rds",sep = ".")), 
                                   formType = combTry$combined)

    if (fileExists) {
      if (out.model.df.length == 0) {
        out.model.df <- out.model.df.row
        out.model.df.length <- out.model.df.length + 1
      } else {
        out.model.df <- rbind(out.model.df, out.model.df.row)
        out.model.df.length <- out.model.df.length + 1
      }
    }
  }
  
  out.model.df$iter <- iter
  
  return(out.model.df)
}