getAUC <- function(list_of_preds){
  outList_mean <- (sapply(seq_along(list_of_preds), function(x) {
    sapply(seq_along(list_of_preds[[1]][[1]])[-1], function(y) {
      mean(AUC(list(list_of_preds[[x]][[1]][[y]],
                    list_of_preds[[x]][[2]][[y]],
                    list_of_preds[[x]][[3]][[y]],
                    list_of_preds[[x]][[4]][[y]],
                    list_of_preds[[x]][[5]][[y]]),
               list(list_of_preds[[x]][[1]]$ground_truth,
                    list_of_preds[[x]][[2]]$ground_truth,
                    list_of_preds[[x]][[3]]$ground_truth,
                    list_of_preds[[x]][[4]]$ground_truth,
                    list_of_preds[[x]][[5]]$ground_truth)))})
  }))
  outList_sd <- (sapply(seq_along(list_of_preds), function(x) {
    sapply(seq_along(list_of_preds[[1]][[1]])[-1], function(y) {
      sd(AUC(list(list_of_preds[[x]][[1]][[y]],
                  list_of_preds[[x]][[2]][[y]],
                  list_of_preds[[x]][[3]][[y]],
                  list_of_preds[[x]][[4]][[y]],
                  list_of_preds[[x]][[5]][[y]]),
             list(list_of_preds[[x]][[1]]$ground_truth,
                  list_of_preds[[x]][[2]]$ground_truth,
                  list_of_preds[[x]][[3]]$ground_truth,
                  list_of_preds[[x]][[4]]$ground_truth,
                  list_of_preds[[x]][[5]]$ground_truth)))})
  }))
  names(outList_mean) <- names(list_of_preds)
  names(outList_sd) <- names(list_of_preds)
  return(list(auc_means=outList_mean, auc_stdevs=outList_sd))
}