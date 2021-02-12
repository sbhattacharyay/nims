library(pROC)

macro.average.AUROC <- function(results.df){
  unique.true.labels <- unique(results.df$true.labels)
  aurocs <- c()
  for (curr.label in unique.true.labels){
    temp.df <- data.frame(binary.true.labels = as.numeric(results.df$true.labels == curr.label)) %>%
      mutate(pred.prob = results.df[[paste0('GCSm=',as.character(curr.label))]])
    roc.obj <- roc(temp.df$binary.true.labels, temp.df$pred.prob, direction = '<',quiet = TRUE)
    aurocs <- c(aurocs, auc(roc.obj))
  }
  
  return(mean(aurocs))
}