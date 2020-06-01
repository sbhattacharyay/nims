get_NPV_info <- function(cps_set){
  outputDF <- data.frame()
  for (i in 1:length(cps_set)){
    curr_pred_space <- cps_set[[i]]
    curr_pred_space_name <- names(cps_set)[i]
    no_classifiers <- length(curr_pred_space$Fold1)
    classifier_names <- names(curr_pred_space$Fold1)
    value_DF <- data.frame()
    for (j in 1:no_classifiers){
      values <- c()
      for (m in 1:length(curr_pred_space)){
        curr_cp <- curr_pred_space[[m]][[j]]
        curr_sens <- curr_cp$sensitivity
        curr_prev <- curr_cp$prevalence
        curr_spec <- curr_cp$specificity
        curr_NPV <- (curr_spec*(1-curr_prev))/(((1-curr_sens)*curr_prev)+((1-curr_prev)*curr_spec))
        values <- c(values,curr_NPV)
      }
      newRow<-data.frame(classifier = classifier_names[j], 
                         mean_value = mean(values),
                         std_value = sd(values),
                         values)
      value_DF <- rbind(value_DF,newRow)
    }
    outputDF <- rbind(outputDF,data.frame(value_DF,pred_space=curr_pred_space_name))
  }
  return(outputDF)
}