# Function to retrieve the predictions and results from machine learning model output files

generateModelDF <- function(path.D = "../all_motion_feature_data/ML_results") {
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
    fileExists <- file.exists(paste(path, "/results/", combTry$MLmethod, "/ROC.tuningML.", combTry$MLmethod,".day", combTry$day, combTry$cumulative, ".iter", combTry$iter, ".results.RDS", sep = ""))
    
    out.model.df.row <- data.frame(method = combTry$MLmethod, 
                                   day = combTry$day, 
                                   iter = combTry$iter, 
                                   classifierPath = paste(path, "/tuning/", combTry$MLmethod, "/ROC.tuningML.", combTry$MLmethod,".day", combTry$day, combTry$cumulative, ".iter", combTry$iter, ".RDS", sep = ""), 
                                   predPath = paste(path, "/pred/", combTry$MLmethod, "/ROC.tuningML.", combTry$MLmethod,".day", combTry$day, combTry$cumulative, ".iter", combTry$iter, ".pred.RDS", sep = ""), 
                                   resultsPath = paste(path, "/results/", combTry$MLmethod, "/ROC.tuningML.", combTry$MLmethod,".day", combTry$day, combTry$cumulative, ".iter", combTry$iter, ".results.RDS", sep = ""), 
                                   stringsAsFactors = FALSE,
                                   cumulative = ifelse(combTry$cumulative == "", FALSE, TRUE))
    
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
  
  return(out.model.df)
}

