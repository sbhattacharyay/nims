get_auc_info <- function(list_of_preds){
  
  outList <- (lapply(seq_along(list_of_preds), function(x) {
    lapply(seq_along(list_of_preds[[1]][[1]])[-1], function(y) {
      cvAUC(list(list_of_preds[[x]][[1]][[y]],
                 list_of_preds[[x]][[2]][[y]],
                 list_of_preds[[x]][[3]][[y]],
                 list_of_preds[[x]][[4]][[y]],
                 list_of_preds[[x]][[5]][[y]]),
            list(list_of_preds[[x]][[1]]$ground_truth,
                 list_of_preds[[x]][[2]]$ground_truth,
                 list_of_preds[[x]][[3]]$ground_truth,
                 list_of_preds[[x]][[4]]$ground_truth,
                 list_of_preds[[x]][[5]]$ground_truth))})
  }))
  names(outList) <- names(list_of_preds)
  lapply(seq_along(outList), function(x){names(outList[[x]]) <<- names(list_of_preds[[x]][[1]])[-1]})
  return(outList)
}