get_sens_info <- function(cps_set){
  outputList <- list()
  for (i in 1:length(cps_set)){
    curr_pred_space <- cps_set[[i]]
    no_classifiers <- length(curr_pred_space$Fold1)
    classifier_names <- names(curr_pred_space$Fold1)
    value_DF <- data.frame()
    for (j in 1:no_classifiers){
      vec_of_values <- c()
      for (m in 1:length(curr_pred_space)){
        curr_cp <- curr_pred_space[[m]][[j]]
        vec_of_values <- c(vec_of_values,curr_cp$sensitivity)
      }
      newRow<-data.frame(classifier = classifier_names[j], 
                         mean_sensitivity = mean(vec_of_values),
                         std_sensitivity = sd(vec_of_values))
      value_DF <- rbind(value_DF,newRow)
    }
    outputList[[i]] <- value_DF
  }
  names(outputList) <- names(cps_set)
  return(outputList)
}